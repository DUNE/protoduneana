#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/TrackHitMeta.h"

#include <TTree.h>
#include <TBranch.h>

#include <cstdio>
#include <algorithm>


#define LOG(x)         (fLog ? printf("\tLOG " #x ": " "\033[1;9%dm" "%s" "\033[0m\n", x?2:1, x?"true":"false") : 0, x)
#define ASSERT(x) if (!(fLog ? printf("\tAST " #x ": " "\033[1;9%dm" "%s" "\033[0m\n", x?2:1, x?"true":"false") : 0, x)) continue; 
#define DEBUG(x)  if  ((fLog ? printf("\tDBG " #x ": " "\033[1;9%dm" "%s" "\033[0m\n", x?1:2, x?"true":"false") : 0, x)) exit(EXIT_FAILURE);

using PtrHit    = art::Ptr<recob::Hit>;
using VecPtrHit = std::vector<art::Ptr<recob::Hit>>;
using PtrTrk    = art::Ptr<recob::Track>;
using VecPtrTrk = std::vector<art::Ptr<recob::Track>>;
using PtrShw    = art::Ptr<recob::Shower>;
using VecPtrShw = std::vector<art::Ptr<recob::Shower>>;

namespace ana {
    /*
    PDVD: (beam direction along Z, inside side1)
    ┌─────────┬─────────┬─────────┬─────────┐  X
    │ side1   │ side1   │ side1   │ side1   │  🡩
    │ sec0    │ sec1    │ sec2    │ sec3    │  │  TOP
    │ tpc8,10 │ tpc9,11 │ tpc12,14│ tpc13,15│  │
    │         │         │         │         │  │
    ├─────────┼─────────┼─────────┼─────────┤  │
    │ side0   │ side0   │ side0   │ side0   │  │
    │ sec4    │ sec5    │ sec6    │ sec7    │  │  BOT
    │ tpc0,2  │ tpc1,3  │ tpc4,6  │ tpc5,7  │  │ 
    │         │         │         │         │  │
    └─────────┴─────────┴─────────┴─────────┘
     ─────────────────────────────────────> Y

    PDHD: (beam direction along Z, inside side0)
    ┌─────────┬─────────┐ Y
    │ side0   │ side1   │ 🡩
    │ sec0    │ sec1    │ │
    │ tpc1,5  │ tpc2,6  │ │
    │         │         │ │
    └─────────┴─────────┘  
     ─────────────────> X
       BOT       TOP
    */

    enum EnumDet { kPDVD, kPDHD };
    std::vector<unsigned> const n_sec = {
        8, // PDVD
        2  // PDHD
    };
    std::vector<std::map<unsigned, int>> const tpc2sec = {
        { // PDVD
            {0, 4}, {2, 4},
            {1, 5}, {3, 5},
            {4, 6}, {6, 6},
            {5, 7}, {7, 7},
            {8, 0}, {10, 0},
            {9, 1}, {11, 1},
            {12, 2}, {14, 2},
            {13, 3}, {15, 3}
        }, { // PDHD
            {4, -1}, {0, -1},
            {5, 0}, {1, 0},
            {6, 1}, {2, 1},
            {7, -1}, {3, -1}
        }
    };
    std::vector<std::map<int, std::pair<unsigned, unsigned>>> const sec2tpc = {
        { // PDVD
            {0, {8, 10} },
            {1, {9, 11} },
            {2, {12, 14} },
            {3, {13, 15} },
            {4, {0, 2} },
            {5, {1, 3} },
            {6, {4, 6} },
            {7, {5, 7} }
        }, { // PDHD
            {0, {1, 5} },
            {1, {2, 6} }
        }
    };
    std::vector<std::map<unsigned, int>> const tpc2side = {
        { // PDVD
            { 0, 0 }, { 1, 0 }, { 2, 0 }, { 3, 0 },
            { 4, 0 }, { 5, 0 }, { 6, 0 }, { 7, 0 },
            { 8, 1 }, { 9, 1 }, {10, 1 }, {11, 1 },
            {12, 1 }, {13, 1 }, {14, 1 }, {15, 1 }
        }, { // PDHD
            { 0, -1 }, { 1, 0 }, { 2, 1 }, { 3, -1 },
            { 4, -1 }, { 5, 0 }, { 6, 1 }, { 7, -1 }
        }
    };
    std::vector<std::map<int, std::vector<int>>> const side2secs = {
        { // PDVD
            { 0, {4, 5, 6, 7} }, { 1, {0, 1, 2, 3} }
        }, { // PDHD
            { 0, {0} }, { 1, {1} }
        }
    };
    std::vector<std::map<int, int>> const sec2side = {
        { // PDVD
            {0, 1}, {1, 1}, {2, 1}, {3, 1},
            {4, 0}, {5, 0}, {6, 0}, {7, 0},
        }, { // PDHD
            {0, 0}, {1, 1} 
        }
    };

    // y = m*x + p
    struct LinearRegression {
        static constexpr unsigned nmin = 4;
        unsigned n=0;
        double mx=0, my=0, mx2=0, my2=0, mxy=0;
        double varx=0, vary=0, cov=0;
        double m=0, p=0, r2=0;
        double lp=0, corr=0;
        void add(double x, double y) {
            mx+=x; my+=y; mx2+=x*x; my2+=y*y; mxy+=x*y; n++;
        }
        void compute() {
            if (n < nmin) return;
            mx/=n; my/=n; mx2/=n; my2/=n; mxy/=n;
            varx=mx2-mx*mx;
            vary=my2-my*my;
            cov=mxy-mx*my;
            m= n<nmin ? 0 : cov/varx;
            p= my - m*mx;
            r2= n<nmin ? 0 : cov*cov / (varx*vary);
            lp= .5*(varx+vary + sqrt(pow(varx-vary, 2) + 4*cov*cov));
            corr= 1 - 4 * (varx*vary - cov*cov) / pow(varx+vary, 2);
        }
        double projection(double x, double y) const {
            return (x + m*(y-p)) / sqrt(1 + m*m);
        }
        double distance(double x, double y) const {
            return abs(y - (m*x + p)) / sqrt(1 + m*m);
        }
        double theta(int dirx) const {
            return atan2(dirx*cov, dirx*(lp-vary));
        }
        void SetBranches(TTree* t, const char* pre="") {
            t->Branch(Form("%sRegM", pre), &m);
            t->Branch(Form("%sRegP", pre), &p);
            t->Branch(Form("%sRegR2", pre), &r2);
        }
    };

    template<typename T>
    struct Bounds {
        T min, max;
        Bounds() : min(std::numeric_limits<T>::max()), max(std::numeric_limits<T>::lowest()) {}
        Bounds(T m, T M) : min(m), max(M) {}
        bool isInside(T x, float r=0) const { 
            return min+r <= x && x <= max-r;
        }

        friend std::ostream& operator<<(std::ostream& os, const Bounds& b) {
            return os << "[" << b.min << ", " << b.max << "]";
        }
    };
    template<typename T>
    struct Bounds3D {
        Bounds<T> x, y, z;
        Bounds3D() : x(), y(), z() {}
        Bounds3D(geo::Point_t const& min, geo::Point_t const& max) :
            x{T(min.x()), T(max.x())}, y{T(min.y()), T(max.y())}, z{T(min.z()), T(max.z())} {}
        // Bounds3D(geo::BoxBoundedGeo const& bb) :
        //     x{bb.MinX(), bb.MaxX()}, y{bb.MinY(), bb.MaxY()}, z{bb.MinZ(), bb.MaxZ()} {}
        bool isInside(geo::Point_t const& p, float r=0) const {
            return x.isInside(p.x(), r) && y.isInside(p.y(), r) && z.isInside(p.z(), r);
        }
        bool isInside(TVector3 const& p, float r=0) const {
            return x.isInside(p.x(), r) && y.isInside(p.y(), r) && z.isInside(p.z(), r);
        }
        bool isInsideYZ(geo::Point_t const& p, float r=0) const {
            return y.isInside(p.y(), r) && z.isInside(p.z(), r);
        }
        geo::Point_t min() const {
            return geo::Point_t{x.min, y.min, z.min};
        }
        geo::Point_t max() const {
            return geo::Point_t{x.max, y.max, z.max};
        }

        friend std::ostream& operator<<(std::ostream& os, const Bounds3D& b) {
            return os << b.min() << " -> " << b.max();
        }
    };

    // axis of equation ay*Y - az*Z = 0
    struct Axis {
        double ay, az;
        double space(geo::WireGeo wiregeo) {
            return ay * wiregeo.GetCenter().Y() + az * wiregeo.GetCenter().Z();
        }
    };

    struct Vec2 {
        float space, drift;
        Vec2() : space(0), drift(0) {}
        Vec2(float s, float d) : space(s), drift(d) {}
        float norm() const { return sqrt(space*space + drift*drift); }
        float angle() const { return atan2(drift, space); }    
        float dot(Vec2 const& v) const { return space*v.space + drift*v.drift; }
        Vec2 operator-(Vec2 const& v) const { return Vec2{space - v.space, drift - v.drift}; }
        void SetBranches(TTree* t, const char* pre="") {
            t->Branch(Form("%sSpace", pre), &space);
            t->Branch(Form("%sDrift", pre), &drift);
        }
    };
    struct Hit {
        unsigned tpc;
        int section;
        float space;
        unsigned channel;
        float tick;
        float adc;
        Hit() : tpc(0), section(0), space(0), channel(0), tick(0), adc(0) {}
        Hit(unsigned T, int S, float s, unsigned c, float t, float a) :
            tpc(T), section(S), space(s), channel(c), tick(t), adc(a) {}

        Vec2 vec(float tick2cm=1) const { return Vec2{space, tick * tick2cm}; }
        // double operator-(Hit const& h) const {
        //     if (section != h.section) return std::numeric_limits<double>::max();
        //     return sqrt(pow(space - h.space, 2) + pow(tick - h.tick, 2));
        // }
        friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
            return os << "tpc:" << hit.tpc << " space:" << hit.space << " ch:" << hit.channel << " tick:" << hit.tick << " ADC:" << hit.adc;
        }
        void SetBranches(TTree* t, const char* pre="") {
            t->Branch(Form("%sHitTPC", pre), &tpc);
            t->Branch(Form("%sHitSection", pre), &section);
            t->Branch(Form("%sHitSpace", pre), &space);
            t->Branch(Form("%sHitChannel", pre), &channel);
            t->Branch(Form("%sHitTick", pre), &tick);
            t->Branch(Form("%sHitADC", pre), &adc);
        }
    };
    struct Hits {
        unsigned N;
        std::vector<unsigned> tpc;
        std::vector<int> section;
        std::vector<float> space;
        std::vector<unsigned> channel;
        std::vector<float> tick;
        std::vector<float> adc;
        Hits() : N(0), tpc(), section(), space(), channel(), tick(), adc() {}
        void push_back(Hit const& hit) {
            N++;
            tpc.push_back(hit.tpc);
            section.push_back(hit.section);
            space.push_back(hit.space);
            channel.push_back(hit.channel);
            tick.push_back(hit.tick);
            adc.push_back(hit.adc);
        }
        void clear() {
            N = 0;
            tpc.clear();
            section.clear();
            space.clear();
            channel.clear();
            tick.clear();
            adc.clear();
        }
        unsigned size() const { return N; }
        bool empty() const { return !N; }
        float energy() const {
            float e = 0;
            for (float a : adc) e+=a;
            return e;
        }
        Vec2 barycenter(float tick2cm=1) const {
            float bs=0, bt=0;
            unsigned bn=0;
            for (unsigned i=0; i<N; i++) {
                bs += this->space[i];
                bt += this->tick[i] * tick2cm;
                bn++;
            }
            return bn ? Vec2{bs/bn, bt/bn} : Vec2{};
        }

        void SetBranches(TTree *t, const char* pre="") {
            t->Branch(Form("%sNHit", pre), &N);
            t->Branch(Form("%sHitTPC", pre), &tpc);
            t->Branch(Form("%sHitSection", pre), &section);
            t->Branch(Form("%sHitSpace", pre), &space);
            t->Branch(Form("%sHitChannel", pre), &channel);
            t->Branch(Form("%sHitTick", pre), &tick);
            t->Branch(Form("%sHitADC", pre), &adc);
        }

        Hit at(unsigned i) const { return Hit{tpc[i], section[i], space[i], channel[i], tick[i], adc[i]}; }
        struct iterator {
            const Hits* hits;
            unsigned i;
            iterator(Hits const* h, unsigned i=0) : hits(h), i(i) {}
            iterator& operator++() { i++; return *this; }
            bool operator!=(iterator const& h) const { return i != h.i; }
            Hit operator*() const { return hits->at(i); }
        };
        iterator begin() const { return iterator(this); }
        iterator end() const { return iterator(this, N); }
        
        bool find(recob::Hit hit) {
            auto it_chan = std::find(channel.begin(), channel.end(), hit.Channel());
            auto it_tick = std::find(tick.begin(), tick.end(), hit.PeakTime());
            unsigned i = std::distance(channel.begin(), it_chan);
            if (i == std::distance(tick.begin(), it_tick)) return i;
            return N;
        }
    };
    struct Point {
        float x, y, z;
        Point() : x(0), y(0), z(0) {}
        Point(float x, float y, float z) : x(x), y(y), z(z) {}
        Point(double x, double y, double z) : x(float(x)), y(float(y)), z(float(z)) {}
        Point(geo::Point_t const& p) : x(float(p.x())), y(float(p.y())), z(float(p.z())) {}
        Point(TVector3 const& v) : x(float(v.X())), y(float(v.Y())), z(float(v.Z())) {}
        float r2() const { return x*x + y*y + z*z; }

        Point operator+(Point const& p) const { return Point{x+p.x, y+p.y, z+p.z}; }
        Point operator+(geo::Point_t const& p) const { return Point{x+(float)p.x(), y+(float)p.y(), z+(float)p.z()}; }
        Point operator+=(Point const& p) { x += p.x; y += p.y; z += p.z; return *this; }
        Point operator-(Point const& p) const { return Point{x-p.x, y-p.y, z-p.z}; }
        Point operator-(geo::Point_t const& p) const { return Point{x-(float)p.x(), y-(float)p.y(), z-(float)p.z()}; }
        Point operator-=(Point const& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
        Point operator*(float f) const { return Point{f*x, f*y, f*z}; }

        operator bool() const { return !(x==0 && y==0 && z==0); }

        friend std::ostream& operator<<(std::ostream& os, const Point& p) {
            return os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
        }

        void SetBranches(TTree* t, const char* pre="") {
            t->Branch(Form("%sPointX", pre), &x);
            t->Branch(Form("%sPointY", pre), &y);
            t->Branch(Form("%sPointZ", pre), &z);
        }
    };
    struct Points {
        unsigned N;
        std::vector<float> x, y, z;
        Points() : N(0), x(), y(), z() {}
        void push_back(Point const& p) {
            N++;
            x.push_back(p.x);
            y.push_back(p.y);
            z.push_back(p.z);
        }
        void push_back(geo::Point_t const& p) { push_back(Point{p}); }
        void clear() {
            N = 0;
            x.clear();
            y.clear();
            z.clear();
        }
        unsigned size() const { return N; }
        bool empty() const { return !N; }
        Point barycenter() const {
            Point b;
            for (unsigned i=0; i<N; i++) b += at(i);
            return b * (1.F / N);
        }


        void SetBranches(TTree* t, const char* pre="") {
            t->Branch(Form("%sNPoint", pre), &N);
            t->Branch(Form("%sPointX", pre), &x);
            t->Branch(Form("%sPointY", pre), &y);
            t->Branch(Form("%sPointZ", pre), &z);
        }

        Point at(unsigned i) const { return Point{x[i], y[i], z[i]}; }
        struct iterator {
            const Points* points;
            unsigned i;
            iterator(Points const* p, unsigned i=0) : points(p), i(i) {}
            iterator& operator++() { i++; return *this; }
            bool operator!=(iterator const& p) { return i != p.i; }
            Point operator*() { return points->at(i); }
        };
        iterator begin() const { return iterator(this); }
        iterator end() const { return iterator(this, N); }
    };

    template<typename AnaStruct>
    void SetBranches(TTree* t, const char* pre, AnaStruct* x) {   x->SetBranches(t, pre); }

    // reco to truth functions
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    // Computationally intensive!!
    template<typename Data> simb::MCParticle const* trk2mcp(
        PtrTrk const& pt,
        detinfo::DetectorClocksData const& clockData,
        art::FindManyP<recob::Hit, Data> const& fmp_trk2hit
    ) {
        std::unordered_map<int, float> map_tid_ene;
        for (PtrHit const& p_hit : fmp_trk2hit.at(pt.key()))
            for (sim::TrackIDE ide : bt_serv->HitToTrackIDEs(clockData, p_hit))
                map_tid_ene[ide.trackID] += ide.energy;

        float max_ene = -1;
        int tid_max = 0;
        for (std::pair<int, float> p : map_tid_ene)
            if (p.second > max_ene)
                max_ene = p.second, tid_max = p.first;
        return max_ene == -1 ? nullptr : pi_serv->TrackIdToParticle_P(tid_max);
    }

    // Computationally intensive!!
    template<typename Data> PtrTrk mcp2trk(
        simb::MCParticle const* mcp, 
        VecPtrTrk const& vp_trk,
        detinfo::DetectorClocksData const& clockData,
        art::FindManyP<recob::Hit, Data> const& fmp_trk2hit
    ) {
        std::unordered_map<PtrTrk, unsigned> map_trk_nhit;
        for (PtrTrk p_trk : vp_trk)
            map_trk_nhit[p_trk] += bt_serv->TrackIdToHits_Ps(clockData, mcp->TrackId(), fmp_trk2hit.at(p_trk.key())).size();

        unsigned max = 0;
        PtrTrk p_trk_from_mcp;
        for (std::pair<PtrTrk, unsigned> p : map_trk_nhit)
            if (p.second > max)
                max = p.second, p_trk_from_mcp = p.first;
        return p_trk_from_mcp;
    }
    template<typename Data> VecPtrTrk mcp2trks(
        simb::MCParticle const* mcp,
        VecPtrTrk const& vpt_ev,
        detinfo::DetectorClocksData const& clockData,
        art::FindManyP<recob::Hit, Data> const& fmp_trk2hit
    ) {
        if (!mcp) return VecPtrTrk{};
        VecPtrTrk vpt_mcp;
        for (PtrTrk pt_ev : vpt_ev) {
            simb::MCParticle const* mcp_trk = trk2mcp(pt_ev, clockData, fmp_trk2hit);
            if (mcp_trk && mcp_trk->TrackId() == mcp->TrackId())
                vpt_mcp.push_back(pt_ev);
        }
        return vpt_mcp;
    }
    VecPtrHit mcp2hits(
        simb::MCParticle const* mcp,
        VecPtrHit const& vph_ev,
        detinfo::DetectorClocksData const& clockData,
        bool use_eve,
        std::vector<float> *energyFracs = nullptr
    ) {
        if (!mcp) return VecPtrHit{};
        VecPtrHit vph_mcp;
        for (PtrHit ph_ev : vph_ev)
            for (sim::TrackIDE ide : (use_eve
                ? bt_serv->HitToEveTrackIDEs(clockData, ph_ev)
                : bt_serv->HitToTrackIDEs(clockData, ph_ev)
            )) {
                if (ide.trackID == mcp->TrackId()) {
                    vph_mcp.push_back(ph_ev);
                    if (energyFracs) energyFracs->push_back(ide.energyFrac);
                }
            }
        return vph_mcp;
    }

    simb::MCParticle const* shw2mcp(
        PtrShw const& ps,
        detinfo::DetectorClocksData const& clockData,
        art::FindManyP<recob::Hit> const& fmp_shw2hit
    ) {
        std::unordered_map<int, float> map_tid_ene;
        for (PtrHit const& p_hit : fmp_shw2hit.at(ps.key()))
            for (sim::TrackIDE ide : bt_serv->HitToTrackIDEs(clockData, p_hit))
                map_tid_ene[ide.trackID] += ide.energy;

        float max_ene = -1;
        int tid_max = 0;
        for (std::pair<int, float> p : map_tid_ene)
            if (p.second > max_ene)
                max_ene = p.second, tid_max = p.first;
        return max_ene == -1 ? nullptr : pi_serv->TrackIdToParticle_P(tid_max);
    }
    PtrShw mcp2shw(
        simb::MCParticle const* mcp, 
        VecPtrShw const& vps_ev,
        detinfo::DetectorClocksData const& clockData,
        art::FindManyP<recob::Hit> const& fmp_shw2hit
    ) {
        std::unordered_map<PtrShw, unsigned> map_shw_nhit;
        for (PtrShw ps : vps_ev)
            map_shw_nhit[ps] += bt_serv->TrackIdToHits_Ps(clockData, mcp->TrackId(), fmp_shw2hit.at(ps.key())).size();

        unsigned max = 0;
        PtrShw p_shw_from_mcp;
        for (std::pair<PtrShw, unsigned> p : map_shw_nhit)
            if (p.second > max)
                max = p.second, p_shw_from_mcp = p.first;
        return p_shw_from_mcp;
    }
    VecPtrShw mcp2shws(
        simb::MCParticle const* mcp,
        VecPtrShw const& vps_ev,
        detinfo::DetectorClocksData const& clockData,
        art::FindManyP<recob::Hit> const& fmp_shw2hit
    ) {
        if (!mcp) return VecPtrShw{};
        VecPtrShw vps_mcp;
        for (PtrShw ps_ev : vps_ev) {
            simb::MCParticle const* mcp_shw = shw2mcp(ps_ev, clockData, fmp_shw2hit);
            if (mcp_shw && mcp_shw->TrackId() == mcp->TrackId())
                vps_mcp.push_back(ps_ev);
        }
        return vps_mcp;
    }




    // Analysis Module Subclass
    struct SortedHits {
        // std::vector<VecPtrHit> vph_sec; // sorted hits per section
        VecPtrHit vph; // sorted hits
        std::vector<int> secs; // sorted sections
        std::vector<LinearRegression> regs; // regressions per side

        std::vector<size_t> secs_first;
        size_t side_first;

        operator bool() const { return !vph.empty(); }
        PtrHit start() const { return vph.front(); }
        PtrHit end() const { return vph.back(); }
        bool is_cc() const { return side_first != 0; }
        PtrHit cc_first() const { return is_cc() ? vph[side_first - 1] : PtrHit(); }
        PtrHit cc_second() const { return is_cc() ? vph[side_first] : PtrHit(); }
        bool is_sc() const { return !secs_first.empty(); }
        VecPtrHit scs() const {
            VecPtrHit res(0);
            for (size_t i : secs_first) {
                if (i != 0 and i != side_first) {
                    res.push_back(vph[i-1]);
                    res.push_back(vph[i]);
                }
            }
            return res;
        }
        int end_sec() const { return secs.back(); }
        LinearRegression end_reg(int det) const {
            return regs[ana::sec2side.at(det).at(end_sec())];
        }
        VecPtrHit::iterator after_cathode_it() { return vph.begin() + side_first; }
        VecPtrHit::iterator endsec_it() { return vph.begin() + secs_first.back(); }

        SortedHits() : vph(0), secs(0), regs(2), secs_first(0), side_first(0) {}

        // PtrHit start, end;
        // std::pair<PtrHit, PtrHit> cc; // cathode crossing
        // VecPtrHit sc; // section crossing

        // SortedHits() : vph(0), secs(0), regs(2) {}
        // operator bool() const { return !vph.empty(); }
        // bool is_cc() const { return cc.first && cc.second; }
        // bool is_sc() const { return !sc.empty(); }
        // int end_sec() const { return secs.back(); }
        // LinearRegression end_reg(int det) const {
        //     return regs[ana::sec2side.at(det).at(end_sec())];
        // }

        // size_t bot_index=0;
        // size_t endsec_index=0;
        // VecPtrHit::iterator bot_it() { return vph.begin() + bot_index; }
        // VecPtrHit::iterator endsec_it() { return vph.begin() + endsec_index; }
    };

    class MichelAnalyzer {
    public:
        MichelAnalyzer(fhicl::ParameterSet const& p);

        // Utilities
        art::ServiceHandle<art::TFileService> asFile;
        // art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        // art::ServiceHandle<cheat::BackTrackerService> bt_serv;

        const geo::GeometryCore* asGeo;
        const geo::WireReadoutGeom* asWire;
        const detinfo::DetectorPropertiesService* asDetProp;
        const detinfo::DetectorClocksService* asDetClocks;
        // error: taking address of rvalue [-fpermissive]
        // const detinfo::DetectorClocksData* clockData;

        int geoDet;
        enum EnumDet { kPDVD, kPDHD };
        float fTick2cm;

        art::InputTag   tag_mcp, tag_sed, tag_wir,
                        tag_hit, tag_clu, tag_trk,
                        tag_shw, tag_spt, tag_pfp;

        // VecPtrHit vph_ev;
        // VecPtrTrk vpt_ev;
        // VecPtrShw vps_ev;
        // error: taking address of rvalue [-fpermissive]
        // std::unique_ptr<art::FindManyP<recob::Hit>> fmp_trk2hit;
        // std::unique_ptr<art::FindOneP<recob::Track>> fop_hit2trk;
        // art::FindManyP<recob::Hit> *fmp_trk2hit;
        // art::FindOneP<recob::Track> *fop_hit2trk;
        // art::FindManyP<recob::Hit> *fmp_shw2hit;
        // art::FindOneP<recob::Shower> *fop_hit2shw;

        bool IsUpright(recob::Track const& T);
        std::string GetParticleName(int pdg);
        Axis GetAxis(geo::PlaneID) const;
        double GetSpace(geo::WireID) const;
        Hit GetHit(PtrHit const&) const;
        double GetDistance(PtrHit const&, PtrHit const&, bool) const;
        double GetDistance(ana::Hit const&, ana::Hit const&, bool) const;
        simb::MCParticle const* GetMichelMCP(simb::MCParticle const*) const;
        std::vector<float> GetdQds(
            VecPtrHit const& vph,
            unsigned smoothing_length,
            unsigned *i_max = nullptr
        ) const;
        std::vector<float> GetdQds(
            VecPtrHit::const_iterator first,
            VecPtrHit::const_iterator last,
            unsigned smoothing_length,
            unsigned *i_max = nullptr
        ) const;

        // SortedHits GetSortedHits(
        //     VecPtrHit const& vph_unsorted,
        //     int dirz = 1,
        //     geo::View_t view = geo::kW
        // ) const;
        SortedHits GetSortedHits_dirX(
            VecPtrHit const& vph_unsorted,
            int dirX = -1, // decreasing X
            geo::View_t view = geo::kW
        ) const;
        // SortedHits GetSortedHits_PDHD(
        //     VecPtrHit const& vph_unsorted,
        //     int dirX,
        //     geo::View_t view = geo::kW
        // ) const;
    };
}

ana::MichelAnalyzer::MichelAnalyzer(fhicl::ParameterSet const& p) :
    asFile()
{
    asGeo = &*art::ServiceHandle<geo::Geometry>{};
    asWire = &art::ServiceHandle<geo::WireReadout>{}->Get();
    asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>{};    
    asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>{};

    for (std::vector<std::string> prod : p.get<std::vector<std::vector<std::string>>>("Products")) {
        std::string 
            process  = prod[0],
            label    = prod[1],
            instance = prod[2],
            type     = prod[3];

        art::InputTag tag(label,instance);

        if      (type == "simb::MCParticle")        tag_mcp = tag;
        else if (type == "sim::SimEnergyDeposit")   tag_sed = tag;
        else if (type == "recob::Wire")             tag_wir = tag;
        else if (type == "recob::Hit")              tag_hit = tag;
        else if (type == "recob::Cluster")          tag_clu = tag;
        else if (type == "recob::Track")            tag_trk = tag;
        else if (type == "recob::Shower")           tag_shw = tag;
        else if (type == "recob::SpacePoint")       tag_spt = tag;
        else if (type == "recob::PFParticle")       tag_pfp = tag;
    }

    if (asGeo->DetectorName().find("vd") != std::string::npos)
        geoDet = kPDVD;
    else if (asGeo->DetectorName().find("hd") != std::string::npos)
        geoDet = kPDHD;
    else {
        std::cout << "\033[1;91m" "unknown geometry: "
            << asGeo->DetectorName() << "\033[0m" << std::endl;
        exit(1);
    }

    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);
    // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();
}    

bool ana::MichelAnalyzer::IsUpright(recob::Track const& T) {
    if (geoDet == kPDVD)
        return T.Start().X() > T.End().X();
    if (geoDet == kPDHD)
        return T.Start().Y() > T.End().Y();
    return false;
}
std::string ana::MichelAnalyzer::GetParticleName(int pdg) {
    std::vector<std::string> static periodic_table = { "",
        "H",                                                                                                  "He", 
        "Li", "Be",                                                             "B",  "C",  "N",  "O",  "F",  "Ne",
        "Na", "Mg",                                                             "Al", "Si", "P",  "S",  "Cl", "Ar",
        "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe"
    };

    switch (pdg) {
        case 11: return "e-";
        case -11: return "e+";
        case 12: return "ve";
        case -12: return "-ve";
        case 13: return "µ-";
        case -13: return "µ+";
        case 14: return "vµ";
        case -14: return "-vµ";
        case 22: return "γ";
        case 2212: return "p";
        case 2112: return "n";
        case 111: return "π0";
        case 211: return "π+";
        case -211: return "π-";
    }

    if (pdg > 1000000000) {
        unsigned ex = pdg % 10;
        unsigned A = (pdg / 10) % 1000;
        unsigned Z = (pdg / 10000) % 1000;
        unsigned L = (pdg / 10000000);
        if (L==100 && Z && Z < periodic_table.size())
            return Form("%u%s%s", A, periodic_table[Z].c_str(), ex ? "*" : "");
    }

    return Form("%d", pdg);
}

ana::Axis ana::MichelAnalyzer::GetAxis(geo::PlaneID pid) const {
    geo::WireGeo w0 = asWire->Wire(geo::WireID{pid, 0});
    geo::WireGeo w1 = asWire->Wire(geo::WireID{pid, 1});
    int dy = w1.GetCenter().Y() > w0.GetCenter().Y() ? 1 : -1;
    int dz = w1.GetCenter().Z() > w0.GetCenter().Z() ? 1 : -1;
    return { dy * w0.CosThetaZ(), dz * w0.SinThetaZ() };
}
double ana::MichelAnalyzer::GetSpace(geo::WireID wid) const {
    return GetAxis(wid).space(asWire->Wire(wid));
}
ana::Hit ana::MichelAnalyzer::GetHit(PtrHit const& ph) const {
    geo::WireID wid = ph->WireID();
    return {
        wid.TPC,
        ana::tpc2sec.at(geoDet).at(wid.TPC),
        float(GetSpace(wid)),
        ph->Channel(),
        ph->PeakTime(),
        ph->ROISummedADC()
    };
}

// 2D projected distance in cm
double ana::MichelAnalyzer::GetDistance(PtrHit const& ph1, PtrHit const& ph2, bool allow_different_side=false) const {
    return GetDistance(GetHit(ph1), GetHit(ph2), allow_different_side);
}
double ana::MichelAnalyzer::GetDistance(ana::Hit const& h1, ana::Hit const& h2, bool allow_different_side=false) const {
    int side1 = ana::sec2side.at(geoDet).at(h1.section);
    int side2 = ana::sec2side.at(geoDet).at(h2.section);
    if (!allow_different_side && (side1 != side2))
        return std::numeric_limits<double>::max();
    return sqrt(pow(h2.space - h1.space, 2) + pow((h2.tick - h1.tick) * fTick2cm, 2));
}

simb::MCParticle const* ana::MichelAnalyzer::GetMichelMCP(
    simb::MCParticle const* muon
) const {
    if (!muon) return nullptr;
    if (muon->NumberDaughters() < 3) return nullptr;

    simb::MCParticle const* michel = nullptr;
    bool has_numu = false, has_nue = false;
    for (
        int i_dau=muon->NumberDaughters()-3;
        i_dau<muon->NumberDaughters();
        i_dau++
    ) {
        simb::MCParticle const* dau = pi_serv->TrackIdToParticle_P(muon->Daughter(i_dau));    
        if (!dau) continue;
        switch (abs(dau->PdgCode())) {
            case 14: has_numu = true; break;
            case 12: has_nue = true; break;
            case 11: michel = dau; break;
            default: break;
        }
    }
    if (has_numu && has_nue && michel) return michel;
    return nullptr;
}

// sort hits, get first and last point,
// points at cathode crossing if any,
// points at section crossing if any
// (a section is a set of two adjacent TPCs)
//
// cause of failure:
// - no section with at least 4 (ana::LinearRegression::nmin) hits 
// ana::SortedHits ana::MichelAnalyzer::GetSortedHits(
//     VecPtrHit const& vph_unsorted,
//     int dirz,
//     geo::View_t view
// ) const {
//     ana::SortedHits sh;
//     // sh.vph_sec.resize(ana::n_sec[geoDet]);
//     std::vector<VecPtrHit> vph_sec(ana::n_sec[geoDet]);
//     std::vector<float> sec_mz(ana::n_sec[geoDet], 0);
//     for (PtrHit const& ph : vph_unsorted) {
//         if (ph->View() != view) continue;
//         int side = ana::tpc2side.at(geoDet).at(ph->WireID().TPC);
//         if (side == -1) continue; // skip hits on the other side of the anodes
//         double z = GetSpace(ph->WireID());
//         double t = ph->PeakTime() * fTick2cm;
//         sh.regs[side].add(z, t);
//         int sec = ana::tpc2sec.at(geoDet).at(ph->WireID().TPC);
//         vph_sec[sec].push_back(ph);
//         sec_mz[sec] += z;
//     }
//     sh.regs[0].compute();
//     sh.regs[1].compute();

//     // sec order according to the direction in Z
//     for (unsigned sec=0; sec<ana::n_sec[geoDet]; sec++) {
//         if (vph_sec[sec].size() < ana::LinearRegression::nmin) {
//             vph_sec[sec].clear();
//             sec_mz[sec] = 0;
//         } else {
//             sec_mz[sec] /= vph_sec[sec].size();
//             sh.secs.push_back(sec);
//         }
//     }
//     if (sh.secs.empty()) return ana::SortedHits{};

//     std::sort(
//         sh.secs.begin(), sh.secs.end(),
//         [&sec_mz, dirz](int sec1, int sec2) {
//             return (sec_mz[sec2] - sec_mz[sec1]) * dirz > 0;
//         }
//     );

//     for (unsigned i=0; i<sh.secs.size(); i++) {
//         int sec = sh.secs[i];
//         ana::LinearRegression const& reg = sh.regs[ana::sec2side.at(geoDet).at(sec)];
//         std::sort(
//             vph_sec[sec].begin(), vph_sec[sec].end(),
//             [&](PtrHit const& ph1, PtrHit const& ph2) {
//                 double s1 = reg.projection(GetSpace(ph1->WireID()), ph1->PeakTime() * fTick2cm);
//                 double s2 = reg.projection(GetSpace(ph2->WireID()), ph2->PeakTime() * fTick2cm);
//                 return (s2 - s1) * dirz > 0;
//             }
//         );

//         if (sec==sh.secs.front())
//             sh.start = vph_sec[sec].front();
//         else {
//             if (ana::sec2side.at(geoDet).at(sec) != ana::sec2side.at(geoDet).at(sh.secs[i-1])) {
//                 sh.cc.second = vph_sec[sec].front(); 
//             } else
//                 sh.sc.push_back(vph_sec[sec].front());
//         } 
//         if (sec==sh.secs.back())
//             sh.end = vph_sec[sec].back();
//         else {
//             if (ana::sec2side.at(geoDet).at(sec) != ana::sec2side.at(geoDet).at(sh.secs[i+1])) {
//                 sh.cc.first = vph_sec[sec].back(); 
//             } else
//                 sh.sc.push_back(vph_sec[sec].back());
//         }
//     }

//     sh.bot_index = 0;
//     for (int sec : sh.secs) {
//         if (ana::sec2side.at(geoDet).at(sec) == 1) 
//             sh.bot_index += vph_sec[sec].size();
//         if (sec == sh.secs.back())
//             sh.endsec_index = sh.vph.size();

//         for (PtrHit const& ph : vph_sec[sec])
//             sh.vph.push_back(ph);
//     }

//     return sh;
// }

ana::SortedHits ana::MichelAnalyzer::GetSortedHits_dirX(
    VecPtrHit const& vph_unsorted,
    int dirX, 
    geo::View_t view
) const {
    ana::SortedHits sh;
    // sh.vph_sec.resize(ana::n_sec[geoDet]);
    std::vector<VecPtrHit> vph_sec(ana::n_sec[geoDet]);
    std::vector<float> sec_mt(ana::n_sec[geoDet], 0);
    for (PtrHit const& ph : vph_unsorted) {
        if (ph->View() != view) continue;
        int side = ana::tpc2side.at(geoDet).at(ph->WireID().TPC);
        if (side == -1) continue; // skip hits on the other side of the anodes in PDHD
        double z = GetSpace(ph->WireID());
        double t = ph->PeakTime() * fTick2cm;
        sh.regs[side].add(z, t);
        int sec = ana::tpc2sec.at(geoDet).at(ph->WireID().TPC);
        vph_sec[sec].push_back(ph);
        sec_mt[sec] += t;
    }
    sh.regs[0].compute();
    sh.regs[1].compute();

    for (unsigned sec=0; sec<ana::n_sec[geoDet]; sec++) {
        if (vph_sec[sec].size() < ana::LinearRegression::nmin) {
            vph_sec[sec].clear();
            sec_mt[sec] = 0;
        } else {
            sh.secs.push_back(sec);
            sec_mt[sec] /= vph_sec[sec].size();
        }
    }
    if (sh.secs.empty()) return ana::SortedHits{};
    
    // first section in the bot volume
    std::vector<int>::iterator bot_it = std::find_if(
        sh.secs.begin(), sh.secs.end(),
        [&](int sec) { return ana::sec2side.at(geoDet).at(sec) == 0; }
    );

    std::sort(
        sh.secs.begin(), bot_it,
        [&sec_mt, &dirX](int sec1, int sec2) {
            // in top volume, increasing X <-> decreasing tick
            return dirX * (sec_mt[sec1]-sec_mt[sec2]) > 0;
        }
    );
    std::sort(
        bot_it, sh.secs.end(),
        [&sec_mt, &dirX](int sec1, int sec2) {
            // in bot volume, increasing X <-> increasing tick
            return dirX * (sec_mt[sec1]-sec_mt[sec2]) < 0;
        }
    );

    for (unsigned i=0; i<sh.secs.size(); i++) {
        int sec = sh.secs[i];
        int side = ana::sec2side.at(geoDet).at(sec);
        ana::LinearRegression const& reg = sh.regs[side];

        // top volume (side==1): increasing X <-> decreasing tick <-> decreasing s for m>0
        // bot volume (side==0): increasing X <-> increasing tick <-> increasing s for m>0
        int sign = dirX * (side == 0 ? 1 : -1) * (reg.m > 0 ? 1 : -1);

        VecPtrHit& vph = vph_sec[sec];
        std::sort(
            vph.begin(), vph.end(),
            [&](PtrHit const& ph1, PtrHit const& ph2) {
                double s1 = reg.projection(GetSpace(ph1->WireID()), ph1->PeakTime() * fTick2cm);
                double s2 = reg.projection(GetSpace(ph2->WireID()), ph2->PeakTime() * fTick2cm);
                return (s2 - s1) * sign > 0;
            }
        );

        // if (sec==sh.secs.front())
        //     sh.start = vph.front();
        // else {
        //     int prev_side = ana::sec2side.at(geoDet).at(sh.secs[i-1]);
        //     if (prev_side != side) {
        //         sh.cc.second = vph.front(); 
        //     } else
        //         sh.sc.push_back(vph.front());
        // } 
        // if (sec==sh.secs.back())
        //     sh.end = vph.back();
        // else {
        //     int next_side = ana::sec2side.at(geoDet).at(sh.secs[i+1]);
        //     if (side != next_side) {
        //         sh.cc.first = vph.back(); 
        //     } else
        //         sh.sc.push_back(vph.back());
        // }
    }

    sh.vph.clear();
    int prev_side=-1;
    for (int sec : sh.secs) {
        sh.secs_first.push_back(sh.vph.size());

        int side = ana::sec2side.at(geoDet).at(sec);
        if (prev_side != -1 && side != prev_side)
            sh.side_first = sh.vph.size();
        prev_side = side;

        for (PtrHit const& ph : vph_sec[sec])
            sh.vph.push_back(ph);
    }

    // sh.bot_index = 0;
    // sh.vph.clear();
    // for (int sec : sh.secs) {
    //     int side = ana::sec2side.at(geoDet).at(sec);
    //     if (side == 1) 
    //         sh.bot_index += vph_sec[sec].size();
    //     if (sec == sh.secs.back())
    //         sh.endsec_index = sh.vph.size();

    //     for (PtrHit const& ph : vph_sec[sec])
    //         sh.vph.push_back(ph);
    // }

    return sh;
}

// ana::SortedHits ana::MichelAnalyzer::GetSortedHits_PDHD(
//     VecPtrHit const& vph_unsorted,
//     int dirX,
//     geo::View_t view
// ) const {
//     ana::SortedHits sh;
//     // sh.vph_sec.resize(ana::n_sec[geoDet]);
//     std::vector<VecPtrHit> vph_sec(ana::n_sec[geoDet]);
//     std::vector<float> sec_mt(ana::n_sec[geoDet], 0);
//     for (PtrHit const& ph : vph_unsorted) {
//         if (ph->View() != view) continue;
//         int side = ana::tpc2side.at(geoDet).at(ph->WireID().TPC);
//         if (side == -1) continue; // skip hits on the other side of the anodes
//         double z = GetSpace(ph->WireID());
//         double t = ph->PeakTime() * fTick2cm;
//         sh.regs[side].add(z, t);
//         int sec = ana::tpc2sec.at(geoDet).at(ph->WireID().TPC);
//         vph_sec[sec].push_back(ph);
//         sec_mt[sec] += t;
//     }
//     sh.regs[0].compute();
//     sh.regs[1].compute();

//     for (unsigned sec=0; sec<ana::n_sec[geoDet]; sec++) {
//         if (vph_sec[sec].size() < ana::LinearRegression::nmin) {
//             vph_sec[sec].clear();
//             sec_mt[sec] = 0;
//         } else {
//             sh.secs.push_back(sec);
//             sec_mt[sec] /= vph_sec[sec].size();
//         }
//     }
//     if (sh.secs.empty()) return ana::SortedHits{};

//     // if the track crosses the cathode, there are hits in both sections
//     // sort them according to the provided X direction
//     if (sh.secs.size() == 2 || dirX < 0) {
//         std::swap(sh.secs[0], sh.secs[1]);
//     }

//     for (unsigned i=0; i<sh.secs.size(); i++) {
//         int sec = sh.secs[i];
//         int side = ana::sec2side.at(geoDet).at(sec);
//         ana::LinearRegression const& reg = sh.regs[side];

//         // dirX > 0: left to right <-> increasing tick <-> increasing s for m>0
//         // dirX < 0: right to left <-> decreasing tick <-> decreasing s for m>0
//         int sign = (dirX > 0 ? 1 : -1) * (reg.m > 0 ? 1 : -1);

//         VecPtrHit& vph = vph_sec[sec];
//         std::sort(
//             vph.begin(), vph.end(),
//             [&](PtrHit const& ph1, PtrHit const& ph2) {
//                 double s1 = reg.projection(GetSpace(ph1->WireID()), ph1->PeakTime() * fTick2cm);
//                 double s2 = reg.projection(GetSpace(ph2->WireID()), ph2->PeakTime() * fTick2cm);
//                 return (s2 - s1) * sign > 0;
//             }
//         );

//         if (sec==sh.secs.front())
//             sh.start = vph.front();
//         else {
//             int prev_side = ana::sec2side.at(geoDet).at(sh.secs[i-1]);
//             if (prev_side == 1 && side == 0) {
//                 sh.cc.second = vph.front(); 
//             } else
//                 sh.sc.push_back(vph.front());
//         } 
//         if (sec==sh.secs.back())
//             sh.end = vph.back();
//         else {
//             int next_side = ana::sec2side.at(geoDet).at(sh.secs[i+1]);
//             if (side == 1 && next_side == 0) {
//                 sh.cc.first = vph.back(); 
//             } else
//                 sh.sc.push_back(vph.back());
//         }
//     }

//     sh.bot_index = 0; // MEANINGLESS IN PDHD
//     sh.vph.clear();
//     for (int sec : sh.secs) {
//         int side = ana::sec2side.at(geoDet).at(sec);
//         if (sec == sh.secs.back())
//             sh.endsec_index = sh.vph.size();

//         for (PtrHit const& ph : vph_sec[sec])
//             sh.vph.push_back(ph);
//     }

//     return sh;
// }

std::vector<float> ana::MichelAnalyzer::GetdQds(
    VecPtrHit const& vph,
    unsigned smoothing_length,
    unsigned *i_max
) const {
    size_t const length = vph.size();
    if (length <= smoothing_length) return std::vector<float>(length, 0.F);

    float max = -1;
    std::vector<float> v_dQds(smoothing_length, 0.F);
    for (auto iph=vph.begin()+smoothing_length; iph!=vph.end(); iph++) {
        VecPtrHit::const_iterator jph = iph-smoothing_length;

        float dQ = std::accumulate(
            jph, iph, 0.,
            [](float sum, PtrHit const& ph) {
                return sum+ph->ROISummedADC();
            }
        ) / smoothing_length;

        float dx = 0;
        for (auto kph=jph; kph!=iph; kph++) {
            if (kph == vph.begin())
                dx += GetDistance(*kph, *(kph+1));
            else if (kph == vph.end()-1)
                dx += GetDistance(*(kph-1), *kph);
            else
                dx += .5 * (
                    GetDistance(*(kph-1), *kph)
                    + GetDistance(*kph, *(kph+1))
                );
        }
        dx /= smoothing_length;

        float dQds = dQ / dx;

        if (i_max && dQds > max) {
            max = dQds;
            *i_max = std::distance(vph.begin(), iph);
        }

        v_dQds.push_back(dQ / dx);
    }
    for (unsigned i=0; i<smoothing_length; i++)
        v_dQds[i] = v_dQds[smoothing_length];

    return v_dQds;
}

std::vector<float> ana::MichelAnalyzer::GetdQds(
    VecPtrHit::const_iterator first,
    VecPtrHit::const_iterator last,
    unsigned smoothing_length,
    unsigned *i_max
) const {
    ptrdiff_t const length = std::distance(first,last);
    if (length <= smoothing_length) return std::vector<float>(length, 0.F);

    float max = -1;
    std::vector<float> v_dQds(smoothing_length, 0.F);
    for (auto iph=first+smoothing_length; iph!=last; iph++) {
        VecPtrHit::const_iterator jph = iph-smoothing_length;

        float dQ = std::accumulate(
            jph, iph, 0.,
            [](float sum, PtrHit const& ph) {
                return sum+ph->ROISummedADC();
            }
        ) / smoothing_length;

        float dx = 0;
        for (auto kph=jph; kph!=iph; kph++) {
            if (kph == first)
                dx += GetDistance(*kph, *(kph+1));
            else if (kph == last-1)
                dx += GetDistance(*(kph-1), *kph);
            else
                dx += .5 * (
                    GetDistance(*(kph-1), *kph)
                    + GetDistance(*kph, *(kph+1))
                );
        }
        dx /= smoothing_length;

        float dQds = dQ / dx;

        if (i_max && dQds > max) {
            max = dQds;
            *i_max = std::distance(first, iph);
        }

        v_dQds.push_back(dQ / dx);
    }
    for (unsigned i=0; i<smoothing_length; i++)
        v_dQds[i] = v_dQds[smoothing_length];

    return v_dQds;
}