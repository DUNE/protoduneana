#include "utils.h"

/* 
Naming conventions,
some variables are prefixed by their type followed by an underscore:
 - `ph_`  for `art::Ptr<recob::Hit>`
 - `vph_` for `std::vector<art::Ptr<recob::Hit>>`

 - `pt_`  for `art::Ptr<recob::Track>`
 - `vpt_` for `std::vector<art::Ptr<recob::Track>>`

 - `ps_`  for `art::Ptr<recob::Shower>`
 - `vps_` for `std::vector<art::Ptr<recob::Shower>>`

 - `tag_` for `art::InputTag`
 - `mcp_` for `simb::MCParticle*`

 - `sh_` for custom `ana::SortedHits`, a structure containing info on the sorted list of hits of a track
MichelAnalysis class properties are prefixed:
 - `as`  for art services
 - `geo` for geometric information
 - `in`  for input variable, taken from fhicl file

 - `ev`  for reco info about the full event (stored in `evTree`)
 - `mu`  for reco info about the track and michel electron (stored in `muTree`)
 - `tru` for truth info about the track (stored in `muTree`)
 - `mi`  for truth info about the michel electron (stored in `muTree`)

Top means X>0, Bot means X<0 (by reference to PDVD)
*/

namespace ana { class MichelAnalysis; }

class ana::MichelAnalysis: 
    public art::EDAnalyzer, 
    private ana::MichelModule 
{
public:
    explicit MichelAnalysis(fhicl::ParameterSet const& p);
    MichelAnalysis(MichelAnalysis const&) = delete;
    MichelAnalysis(MichelAnalysis&&) = delete;
    MichelAnalysis& operator=(MichelAnalysis const&) = delete;
    MichelAnalysis& operator=(MichelAnalysis&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;
private:
    // Geometry Information
    ana::Bounds<float> wireWindow;
    ana::Bounds3D<float> geoTop, geoBot;
    float geoCathodeGap; // cm

    // Input Parameters
    bool        inLog;
    bool        inKeepAll;
    float       inTrackLengthCut; // in cm
    float       inFiducialLength; // in cm
    float       inBarycenterRadius; // in cm
    float       inMichelRadius; // in cm
    unsigned    inRegN;

    // Output Variables
    TTree       *evTree;
    TTree       *muTree;

    // Event information
    unsigned                evRun;
    unsigned                evSubRun;
    unsigned                evEvent;
    unsigned                evIndex=0;
    unsigned                evMuonNumber;
    std::vector<unsigned>   evMuonIndices;
    bool                    evIsData;
    ana::Hits               evHits;

    // Track information: recob::Track
    unsigned                muIndex=0;
    float                   muLength;
    ana::Point              muStartPoint;
    ana::Point              muEndPoint;

    // Track information: recob::Hit
    ana::Hits               muHits;
    std::vector<float>      muHitCathodeX;
    std::vector<float>      muHitAnodeX;
    std::vector<float>      muHitY;
    ana::Hit                muStartHit; 
    float                   muStartHitCathodeX;
    float                   muStartHitY;
    bool                    muStartInCathodeX;
    bool                    muStartInY;
    bool                    muStartInZ;
    bool                    muStartInT;
    bool                    muStartInYZT;    
    ana::Hit                muEndHit;
    float                   muEndHitCathodeX; 
    float                   muEndHitY;
    bool                    muEndInCathodeX;
    bool                    muEndInY;
    bool                    muEndInZ;
    bool                    muEndInT;
    bool                    muEndInYZT;
    // bool                    muRegError;
    // ana::LinearRegression   muTopReg;
    // ana::LinearRegression   muBotReg;
    float                   muEndAngle;
    bool                    muCathodeCrossing;
    float                   muCathodeMisalignment;
    bool                    muAnodeCrossing;
    std::vector<float>      muHitdQds;
    ana::Hits               muSphereHits;
    float                   muSphereEnergy;
    float                   muSphereEnergyTP;
    std::vector<float>      muSphereHitMuonAngle;
    // float                   muSphereMaxShowerEnergy;

    // Track information: reconstructed Michel electron
    ana::Vec2   muBary;
    ana::Hits   muBaryHits;
    float       muBaryAngle;
    float       muBaryMuonAngle;
    bool        muSphereHasLongTrack;

    // Truth information: Muon
    int                     truPdg;
    std::string             truEndProcess;
    ana::Point              truStartPoint;
    ana::Point              truEndPoint;
    float                   truEndEnergy;
    // ana::Hit                truStartHit;
    // ana::Hit                truEndHit;
    // ana::LinearRegression   truReg;
    // float                   truEndAngle;
    enum EnumHasMichel: int { 
        kHasNoMichel        = 0, 
        kHasMichelOutside   = 1, 
        kHasMichelInside    = 2, 
        kHasMichelFiducial  = 3
    };
    EnumHasMichel           truHasMichel;

    // Truth information: Michel electron
    float               miTrueEnergy;
    float               miTrackLength;
    // float               miShowerLength;
    ana::Hits           miHits;
    std::vector<float>  miHitEnergyFrac;
    // std::vector<float>  miHitMuonAngle;
    float               miHitEnergy;

    // Truth information: Hits nearby muon's end
    // unsigned            miBaryNHit;
    // ana::Vec2           miBary;
    // float               miBaryAngle;
    // float               miBaryMuonAngle;

    void resetEvent(void);
    void resetMuon(void);
};

ana::MichelAnalysis::MichelAnalysis(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} 
    , MichelModule{p}
    , inLog(p.get<bool>("Log", true))
    , inKeepAll(p.get<bool>("KeepAll", true))
    , inTrackLengthCut(p.get<float>("TrackLengthCut", 30.F)) // in cm
    , inFiducialLength(p.get<float>("FiducialLength", 20.F)) // in cm
    , inBarycenterRadius(p.get<float>("BarycenterRadius", 10.F)) // in cm
    , inMichelRadius(p.get<float>("MichelRadius", 20.F)) //in cm
    , inRegN(p.get<unsigned>("RegN", 6))
{
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();
    wireWindow = ana::Bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    switch (geoDet) {
        case kPDVD:
            geoBot = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 0}).Min(),
                asGeo->TPC(geo::TPCID{0, 7}).Max()
            };
            geoTop = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 8}).Min(),
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoBot = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 5}).Max()
            };
            geoTop = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 2}).Min(),
                asGeo->TPC(geo::TPCID{0, 6}).Max()
            };
            break;
        case kPDSP:
            geoBot = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 9}).Max()
            };
            geoTop = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 2}).Min(),
                asGeo->TPC(geo::TPCID{0, 10}).Max()
            };
            break;
        default: break;
    }
    geoCathodeGap = geoTop.x.min - geoBot.x.max;

    std::cout << "MiAna: " "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Detector Geometry: " << asGeo->DetectorName()
        << "  (" << ana::det_name[geoDet] << ")" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  Top Bounds: " << geoTop << std::endl
        << "  Bot Bounds: " << geoBot << std::endl;
    std::cout << "MiAna: " "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Keep All Events: " << (inKeepAll ? "true" : "false") << std::endl
        << "  Track Length Cut: " << inTrackLengthCut << " cm" << std::endl
        << "  Fiducial Length: " << inFiducialLength << " cm" << std::endl
        << "  Barycenter Radius: " << inBarycenterRadius << " cm" << std::endl
        << "  Michel Space Radius: " << inMichelRadius << " cm" << std::endl
        << "  Smoothing Length: " << inRegN << " points" << std::endl;

    evTree = asFile->make<TTree>("event","");

    evTree->Branch("EventRun",      &evRun);
    evTree->Branch("EventSubRun",   &evSubRun);
    evTree->Branch("EventEvent",    &evEvent);
    evTree->Branch("Index",         &evIndex);
    evTree->Branch("MuonNumber",    &evMuonNumber);
    evTree->Branch("MuonIndices",   &evMuonIndices);
    evTree->Branch("IsData",        &evIsData);
    SetBranches(evTree, "",         &evHits);

    muTree = asFile->make<TTree>("muon","");

    // Event
    muTree->Branch("EventRun",      &evRun);
    muTree->Branch("EventSubRun",   &evSubRun);
    muTree->Branch("EventEvent",    &evEvent);
    muTree->Branch("IsData",        &evIsData);
    muTree->Branch("EventIndex",    &evIndex);
    muTree->Branch("IndexInEvent",  &evMuonNumber);
    muTree->Branch("Index",         &muIndex);

    // Track
    muTree->Branch("Length",        &muLength);
    SetBranches(muTree, "Start",    &muStartPoint);
    SetBranches(muTree, "End",      &muEndPoint);

    // Hit
    // muTree->Branch("RegError",                  &muRegError);
    muTree->Branch("CathodeCrossing",           &muCathodeCrossing);
    muTree->Branch("CathodeMisalignment",       &muCathodeMisalignment);
    muTree->Branch("AnodeCrossing",             &muAnodeCrossing);
    SetBranches(muTree, "Start",                &muStartHit);
    muTree->Branch("StartHitCathodeX",          &muStartHitCathodeX);
    muTree->Branch("StartHitY",                 &muStartHitY);
    muTree->Branch("StartInCathodeX",           &muStartInCathodeX);
    muTree->Branch("StartInY",                  &muStartInY);
    muTree->Branch("StartInZ",                  &muStartInZ);
    muTree->Branch("StartInT",                  &muStartInT);
    muTree->Branch("StartInYZT",                &muStartInYZT);
    SetBranches(muTree, "End",                  &muEndHit);
    muTree->Branch("EndHitCathodeX",            &muEndHitCathodeX);
    muTree->Branch("EndHitY",                   &muEndHitY);
    muTree->Branch("EndInCathodeX",             &muEndInCathodeX);
    muTree->Branch("EndInY",                    &muEndInY);
    muTree->Branch("EndInZ",                    &muEndInZ);
    muTree->Branch("EndInT",                    &muEndInT);
    muTree->Branch("EndInYZT",                  &muEndInYZT);
    muTree->Branch("EndAngle",                  &muEndAngle);
    SetBranches(muTree, "",                     &muHits);
    muTree->Branch("HitCathodeX",               &muHitCathodeX);
    muTree->Branch("HitAnodeX",                 &muHitAnodeX);
    muTree->Branch("HitY",                      &muHitY);
    // SetBranches(muTree, "Top",                  &muTopReg);
    // SetBranches(muTree, "Bot",                  &muBotReg);
    muTree->Branch("HitdQds",                   &muHitdQds);
    SetBranches(muTree, "Sphere",               &muSphereHits);
    muTree->Branch("SphereHitMuonAngle",        &muSphereHitMuonAngle);
    muTree->Branch("SphereEnergy",              &muSphereEnergy); // ADC
    muTree->Branch("SphereEnergyTP",            &muSphereEnergyTP); // ADC
    muTree->Branch("SphereHasLongTrack",        &muSphereHasLongTrack);
    // muTree->Branch("SphereMaxShowerEnergy",     &muSphereMaxShowerEnergy);

    SetBranches(muTree, "Bary",                 &muBaryHits);
    SetBranches(muTree, "Bary",                 &muBary);
    muTree->Branch("BaryAngle",                 &muBaryAngle);
    muTree->Branch("BaryMuonAngle",             &muBaryMuonAngle);

    // Truth
    muTree->Branch("TruePdg",               &truPdg);
    muTree->Branch("TrueEndProcess",        &truEndProcess);
    SetBranches(muTree, "TrueStart",        &truStartPoint);
    SetBranches(muTree, "TrueEnd",          &truEndPoint);
    muTree->Branch("TrueEndEnergy",         &truEndEnergy);
    // SetBranches(muTree, "TrueStart",        &truStartHit);
    // SetBranches(muTree, "TrueEnd",          &truEndHit);
    // SetBranches(muTree, "True",             &truReg);
    // muTree->Branch("TrueEndAngle",          &truEndAngle);
    muTree->Branch("TrueHasMichel",   (int*)&truHasMichel);
    muTree->Branch("MichelTrueEnergy",      &miTrueEnergy); // MeV
    muTree->Branch("MichelTrackLength",     &miTrackLength); // cm
    // muTree->Branch("MichelShowerLength",    &miShowerLength); // cm
    SetBranches(muTree, "Michel",           &miHits);
    muTree->Branch("MichelHitEnergyFrac",   &miHitEnergyFrac);
    // muTree->Branch("MichelHitMuonAngle",    &miHitMuonAngle);
    muTree->Branch("MichelHitEnergy",       &miHitEnergy); // ADC

    // muTree->Branch("MichelBaryNHit",        &miBaryNHit);
    // SetBranches(muTree, "MichelBary",       &miBary); 
    // muTree->Branch("MichelBaryAngle",       &miBaryAngle); // rad
    // muTree->Branch("MichelBaryMuonAngle",   &miBaryMuonAngle); // rad
}

void ana::MichelAnalysis::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) {
        std::cout << "MiAna: " "\033[1;91m" "No valid recob::Hit handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrHit vph_ev;
    art::fill_ptr_vector(vph_ev, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) {
        std::cout << "MiAna: " "\033[1;91m" "No valid recob::Track handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrTrk vpt_ev;
    art::fill_ptr_vector(vpt_ev, vh_trk);

    // auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
    // if (!vh_shw.isValid()) {
    //     std::cout << "MiAna: " "\033[1;91m" "No valid recob::Shower handle" "\033[0m" << std::endl;
    //     return;
    // }
    // VecPtrShw vps_ev;
    // art::fill_ptr_vector(vps_ev, vh_shw);

    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
    // art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
    // art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

    resetEvent();

    // dump event information
    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();
    evIsData = e.isRealData();
    for (PtrHit p_hit : vph_ev)
        if (p_hit->View() == geo::kW)
            evHits.push_back(GetHit(p_hit));

    // loop over tracks to find stopping muons
    for (PtrTrk const& pt_ev : vpt_ev) {
        if (inLog) std::cout << "e" << evIndex << "t" << pt_ev->ID() << "\r" << std::flush;
        resetMuon();

        // get hits and metadata associated to the track
        VecPtrHit vph_mu_all = fmp_trk2hit.at(pt_ev.key());
        std::vector<recob::TrackHitMeta const*> const& vhm_mu = fmp_trk2hit.data(pt_ev.key());
        // std::map<size_t, unsigned> map_hitkey2metaidx;
        std::map<size_t, size_t> map_hitkey2trkidx;
        ASSERT(!vph_mu_all.empty())
        ASSERT(vph_mu_all.size() == vhm_mu.size())

        // remove hits not associated to a track point (hit misassociated to the track?)
        VecPtrHit vph_mu;
        for (unsigned i=0; i<vph_mu_all.size(); i++) {
            if (vph_mu_all[i]->View() != geo::kW) continue;
            if (!pt_ev->HasValidPoint(vhm_mu[i]->Index())) continue;
            // map_hitkey2metaidx[vph_mu_all[i].key()] = i;
            map_hitkey2trkidx[vph_mu_all[i].key()] = vhm_mu[i]->Index();
            vph_mu.push_back(vph_mu_all[i]);
        }
        vph_mu_all.clear();
        ASSERT(!vph_mu.empty())

        bool track_is_up =  IsUpright(*pt_ev);
        geo::Point_t Start = track_is_up ? pt_ev->Start() : pt_ev->End();
        geo::Point_t End = track_is_up ? pt_ev->End() : pt_ev->Start();
        muStartPoint = ana::Point(Start);
        muEndPoint = ana::Point(End);

        // sort hits and retrieve some info: mainly if it crosses the cathode
        // don't specify X direction (Start.X()>End.X()? is not reliable)
        // ana::SortedHits sh_mu = GetSortedHits(vph_mu);
        // muRegError = !sh_mu;
        // ASSERT(!muRegError)

        if (inLog) std::cout << "\t" "\033[1;93m" "e" << evIndex << "m" << evMuonNumber << " (" << muIndex << ")" "\033[0m" << std::endl;

        // get truth information about the track's particle and potential michel electron
        simb::MCParticle const* mcp = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        simb::MCParticle const* mcp_mi = nullptr;
        VecPtrHit vph_mcp_mu, vph_mi;
        std::vector<float> energyFracs_mi;
        if (mcp) {
            mcp_mi = GetMichelMCP(mcp);
            vph_mcp_mu = ana::mcp2hits(mcp, vph_ev, clockData, false);
            vph_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, true, &energyFracs_mi);
        }

        // dump basic track information
        muLength = pt_ev->Length();
        // muChi2 = pt_ev->Chi2();
        // muChi2PerNdof = pt_ev->Chi2PerNdof();

        LOG(muLength >= inTrackLengthCut);
        if (!inKeepAll && muLength < inTrackLengthCut) continue;


        std::sort(vph_mu.begin(), vph_mu.end(), [&map_hitkey2trkidx](PtrHit const& ph1, PtrHit const& ph2) {
            return map_hitkey2trkidx.at(ph1.key()) < map_hitkey2trkidx.at(ph2.key());
        });

        bool is_in_bot=false, is_in_top=false;
        for (PtrHit const& ph : vph_mu) {
            if (is_in_bot && is_in_top) break;
            if (!is_in_bot && GetSide(ph) == kBot) is_in_bot = true;
            if (!is_in_top && GetSide(ph) == kTop) is_in_top = true;
        }
        bool is_cc = is_in_bot && is_in_top;

        bool is_up = false;
        if (geoDet == kPDVD) {
            Side_t front_side = GetSide(vph_mu.front());
            if (is_cc)
                is_up = front_side == kTop;
            else if (front_side == kTop)
                is_up = vph_mu.front()->PeakTime() < vph_mu.back()->PeakTime();
            else
                is_up = vph_mu.front()->PeakTime() > vph_mu.back()->PeakTime();
        } else if (geoDet == kPDHD || geoDet == kPDSP) {
            float front_y = pt_ev->LocationAtPoint(map_hitkey2trkidx.at(vph_mu.front().key())).Y();
            float back_y = pt_ev->LocationAtPoint(map_hitkey2trkidx.at(vph_mu.back().key())).Y();
            is_up = front_y > back_y;
        }
        if (!is_up)
            std::reverse(vph_mu.begin(), vph_mu.end());

        Side_t first_side = GetSide(vph_mu.front());
        VecPtrHit::iterator cathode_it = std::find_if(
            vph_mu.begin()+1, vph_mu.end(), [&](PtrHit const& ph) {
                return GetSide(ph) != first_side;
            }
        );
        if (is_cc) { ASSERT(cathode_it != vph_mu.end()) }
        PtrHit const cc_first = is_cc ? *(cathode_it-1) : PtrHit{};
        PtrHit const cc_second = is_cc ? *cathode_it : PtrHit{};

        // // get track orientation (assuming downward), update orientation of sorted hits
        // bool is_up = true;
        // if (geoDet == kPDVD) {
        //     Side_t front_side = GetSide(sh_mu.start());
        //     is_up = sh_mu.is_cc()
        //         ? front_side == kTop
        //         : ( front_side == kTop
        //             ? sh_mu.start()->PeakTime() < sh_mu.end()->PeakTime()
        //             : sh_mu.start()->PeakTime() > sh_mu.end()->PeakTime()
        //         );
        // } else if (geoDet == kPDHD) {
        //     // size_t front_hit_track_idx = vhm_mu.at(map_hitkey2metaidx[sh_mu.start().key()])->Index();
        //     size_t front_hit_track_idx = map_hitkey2trkidx.at(sh_mu.start().key());
        //     float front_hit_y = pt_ev->HasValidPoint(front_hit_track_idx)
        //         ? pt_ev->LocationAtPoint(front_hit_track_idx).Y()
        //         : util::kBogusF;
        //     // size_t back_hit_track_idx = vhm_mu.at(map_hitkey2metaidx[sh_mu.end().key()])->Index();
        //     size_t back_hit_track_idx = map_hitkey2trkidx.at(sh_mu.end().key());
        //     float back_hit_y = pt_ev->HasValidPoint(back_hit_track_idx)
        //         ? pt_ev->LocationAtPoint(back_hit_track_idx).Y()
        //         : util::kBogusF;
        //     is_up = front_hit_y > back_hit_y;
        // }
        // if (!is_up) sh_mu.reverse();

        // dump first hit of the track:
        // muStartHit = GetHit(sh_mu.start());
        muStartHit = GetHit(vph_mu.front());
        // size_t start_track_idx = vhm_mu.at(map_hitkey2metaidx[sh_mu.start().key()])->Index();
        // size_t start_track_idx = map_hitkey2trkidx.at(sh_mu.start().key());
        size_t start_track_idx = map_hitkey2trkidx.at(vph_mu.front().key());
        muStartHitY = pt_ev->HasValidPoint(start_track_idx)
            ? pt_ev->LocationAtPoint(start_track_idx).Y()
            : util::kBogusF;
        muStartInY = geoBot.y.isInside(muStartHitY, inFiducialLength) || geoTop.y.isInside(muStartHitY, inFiducialLength);
        muStartInZ = geoBot.z.isInside(muStartHit.space, inFiducialLength) || geoTop.z.isInside(muStartHit.space, inFiducialLength);
        muStartInT = wireWindow.isInside(muStartHit.tick, inFiducialLength / fTick2cm);

        // dump last hit of the track:
        muEndHit = GetHit(vph_mu.back());
        // size_t end_track_idx = vhm_mu.at(map_hitkey2metaidx[sh_mu.end().key()])->Index();
        size_t end_track_idx = map_hitkey2trkidx.at(vph_mu.back().key());
        muEndHitY = pt_ev->HasValidPoint(end_track_idx)
            ? pt_ev->LocationAtPoint(end_track_idx).Y()
            : util::kBogusF;
        muEndInY = geoBot.y.isInside(muEndHitY, inFiducialLength) || geoTop.y.isInside(muEndHitY, inFiducialLength);
        muEndInZ = geoBot.z.isInside(muEndHit.space, inFiducialLength) || geoTop.z.isInside(muEndHit.space, inFiducialLength);
        muEndInT = wireWindow.isInside(muEndHit.tick, inFiducialLength / fTick2cm);

        // dump linear regression of hits
        // muTopReg = sh_mu.regs[kTop];
        // muBotReg = sh_mu.regs[kBot];
        // LOG(muTopReg.r2 >= 0.5 && muBotReg.r2 >= 0.5);
        // if (!inKeepAll && !(muTopReg.r2 >= 0.5 && muBotReg.r2 >= 0.5)) continue;
        if (geoDet == kPDHD && !inKeepAll && vph_mu.back()->WireID().TPC == 1) continue;


        // dump dQ/ds
        muHitdQds.clear();
        std::vector<float> tmp;
        tmp = GetdQds(vph_mu.begin(), cathode_it, inRegN);
        muHitdQds.insert(muHitdQds.end(), tmp.begin(), tmp.end());
        tmp = GetdQds(cathode_it, vph_mu.end(), inRegN);
        muHitdQds.insert(muHitdQds.end(), tmp.begin(), tmp.end());
        tmp.clear();

        // dump cathode crossing boolean: defined as track with 4+ hits on both sides of the cathode
        muCathodeCrossing = is_cc;
        muCathodeMisalignment = is_cc ? abs(cc_first->PeakTime()-cc_second->PeakTime())*fTick2cm : util::kBogusF;

        LOG(muCathodeCrossing);
        // if (!inKeepAll && !muCathodeCrossing) continue;

        switch (geoDet) {
        case kPDVD: /* ASSUMS DOWNWARD MUON */
            muAnodeCrossing = GetSide(muStartHit.section) == kTop
                && geoTop.y.isInside(muStartHitY, inFiducialLength)
                && geoTop.z.isInside(muStartHit.space, inFiducialLength)
                && wireWindow.isInside(muStartHit.tick, inFiducialLength/fTick2cm);
            break;
        case kPDHD:
        case kPDSP:
            muAnodeCrossing =
                geoTop.y.isInside(muStartHitY, inFiducialLength)
                && geoTop.z.isInside(muStartHit.space, inFiducialLength)
                && wireWindow.isInside(muStartHit.tick, inFiducialLength/fTick2cm);
            break;
        }
        LOG(muAnodeCrossing);
        if (!inKeepAll && !(muAnodeCrossing || muCathodeCrossing)) continue;


        // dump hit X positions for cathode/anode crossing tracks
        if (muCathodeCrossing) {
            Side_t start_side = GetSide(muStartHit.tpc);
            Side_t end_side = GetSide(muEndHit.tpc);
            PtrHit const& cc_bot = start_side == kBot ? cc_first : cc_second;
            PtrHit const& cc_top = start_side == kTop ? cc_first : cc_second;

            muStartHitCathodeX = GetCathodeX(muStartHit, cc_bot, cc_top, geoCathodeGap);
            muStartInCathodeX = start_side == kBot
                ? geoBot.x.isInside(muStartHitCathodeX, inFiducialLength)
                : geoTop.x.isInside(muStartHitCathodeX, inFiducialLength);

            muEndHitCathodeX = GetCathodeX(muEndHit, cc_bot, cc_top, geoCathodeGap);
            muEndInCathodeX = end_side == kBot
                ? geoBot.x.isInside(muEndHitCathodeX, inFiducialLength)
                : geoTop.x.isInside(muEndHitCathodeX, inFiducialLength);

            for (PtrHit const& ph_mu : vph_mu) {
                muHitCathodeX.push_back(GetCathodeX(ph_mu, cc_bot, cc_top, geoCathodeGap));
            }

            // if (mcp) {
            //     ana::SortedHits sh_mcp = GetSortedHits(vph_mcp_mu, mcp->EndX() > mcp->Vx() ? 1 : -1);
            //     if (sh_mcp && sh_mcp.is_cc()) {
            //         Side_t tru_start_side = GetSide(sh_mcp.start());
            //         PtrHit const& tru_cc_bot = tru_start_side == kBot ? sh_mcp.cc_first() : sh_mcp.cc_second();
            //         PtrHit const& tru_cc_top = tru_start_side == kTop ? sh_mcp.cc_first() : sh_mcp.cc_second();
            //         std::cout << "TRUTH" << std::endl;
            //         std::cout << "\t" "(X, Y, Z): " << mcp->Vx() << ", " << mcp->Vy() << ", " << mcp->Vz() << " -> " << mcp->EndX() << ", " << mcp->EndY() << ", " << mcp->EndZ() << std::endl;
            //         std::cout << "\t" "start side: " << tru_start_side << std::endl;
            //         std::cout << "\t" "start_T: " << sh_mcp.start()->PeakTime() << "\t" "end_T: " << sh_mcp.end()->PeakTime() << std::endl;
            //         std::cout << "\t" "cc_bot_T: " << tru_cc_bot->PeakTime() << "\t" "cc_top_T: " << tru_cc_top->PeakTime() << std::endl;

            //         std::cout << "RECO" << std::endl;
            //         std::cout << "\t" "(X, Y, Z): " << Start.X() << ", " << Start.Y() << ", " << Start.Z() << " -> " << End.X() << ", " << End.Y() << ", " << End.Z() << std::endl;
            //         std::cout << "\t" "start_side: " << start_side << std::endl;
            //         std::cout << "\t" "start_T: " << vph_mu.front()->PeakTime() << "\t" "end_T: " << vph_mu.back()->PeakTime() << std::endl;
            //         std::cout << "\t" "start_X: " << muStartHitX << "\t" "end_X: " << muEndHitX << std::endl;
            //         std::cout << "\t" "cc_bot_T: " << cc_bot->PeakTime() << "\t" "cc_top_T: " << cc_top->PeakTime() << std::endl;
            //         std::cout << std::endl;
            //     }
            // }
        }
        if (muAnodeCrossing) {
            for (PtrHit const& ph_mu : vph_mu) {
                muHitAnodeX.push_back(GetAnodeX(ph_mu, vph_mu.front(), geoBot.x.min, geoTop.x.max));
            }
        }

        muStartInYZT = muStartInY && muStartInZ && muStartInT;
        muEndInYZT = muEndInY && muEndInZ && muEndInT;
        LOG(muStartInYZT);
        LOG(muEndInYZT);
        if (!inKeepAll && !muEndInYZT) continue;

        // dump hits
        for (PtrHit const& ph_mu : vph_mu) {
            if (ph_mu->View() != geo::kW) continue;
            ana::Hit hit = GetHit(ph_mu);
            muHits.push_back(hit);
            // size_t hit_track_idx = vhm_mu.at(map_hitkey2metaidx[ph_mu.key()])->Index();
            size_t hit_track_idx = map_hitkey2trkidx.at(ph_mu.key());
            muHitY.push_back(pt_ev->HasValidPoint(hit_track_idx)
                ? pt_ev->LocationAtPoint(hit_track_idx).Y()
                : util::kBogusF
            );
        }

        // get hits from the section of the last hit
        VecPtrHit vph_ev_endsec;
        for (PtrHit const& ph_ev : vph_ev) {
            if (ph_ev->View() != geo::kW) continue;
            if (GetSec(ph_ev->WireID().TPC) != GetSec(vph_mu.back())) continue;
            vph_ev_endsec.push_back(ph_ev);
        }

        // muEndAngle = sh_mu.regs
        //     .at(GetSide(muEndHit.section))
        //     .theta(muEndHit.space > muStartHit.space ? 1 : -1);
        
        
        recob::Track::Vector_t end_dir = pt_ev->EndDirection();
        muEndAngle = atan2(end_dir.x(), end_dir.z());
        // integrate charges around muon endpoint
        muSphereEnergy = 0;
        muSphereEnergyTP = 0;
        // muSphereMaxShowerEnergy = 0;
        for (PtrHit const& ph_ev : vph_ev_endsec) {
            float dist = GetDistance(ph_ev, vph_mu.back());
            if (dist > inMichelRadius) continue;

            PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
            if (pt_hit && pt_hit->Length() > inTrackLengthCut) {
                if (pt_hit.key() != pt_ev.key()) muSphereHasLongTrack = true;
                continue;
            }

            // PtrShw ps_hit = fop_hit2shw.at(ph_ev.key());
            // if (ps_hit)
            //     for (double energy : ps_hit->Energy())
            //         if (energy > muSphereMaxShowerEnergy)
            //             muSphereMaxShowerEnergy = energy;

            if (dist < inBarycenterRadius) {
                muBaryHits.push_back(GetHit(ph_ev));
            }

            ana::Hit hit = GetHit(ph_ev);
            muSphereEnergy += ph_ev->ROISummedADC();
            muSphereHits.push_back(hit);

            // float da = (hit.vec(fTick2cm) - muEndHit.vec(fTick2cm)).angle() - muEndAngle;
            // da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
            // muSphereHitMuonAngle.push_back(da);

            if (std::find_if(
                vph_mi.begin(), vph_mi.end(),
                [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
            ) != vph_mi.end())
                muSphereEnergyTP = ph_ev->ROISummedADC();
        }

        // Cone
        LOG(muBaryHits.size());
        if (muBaryHits.size()) {
            muBary = muBaryHits.barycenter(fTick2cm);
            ana::Vec2 end_bary = muBary - muEndHit.vec(fTick2cm);
            muBaryAngle = end_bary.angle();
            float da = muBaryAngle - muEndAngle;
            da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
            muBaryMuonAngle = da;

            // if (inCone) {
            //     // float angle = end_bary.angle();
            //     PandoraConeEnergy = 0;
            //     PandoraConeEnergyTP = 0;
            //     for (PtrHit const& ph_ev : vph_ev_endsec) {
            //         ana::Hit hit = GetHit(ph_ev);
            //         PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
            //         if (pt_hit && pt_hit->Length() > inTrackLengthCut) continue;

            //         float dist = GetDistance(ph_ev, sh_mu.end);
            //         if (dist > 30) continue;

            //         ana::Vec2 end_hit = hit.vec(fTick2cm) - muEndHit.vec(fTick2cm);
            //         float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

            //         if (dist > 5
            //             && cosa < cos(30.F * TMath::DegToRad())
            //         ) continue;

            //         bool tp = std::find_if(
            //             vph_mi.begin(), vph_mi.end(),
            //             [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
            //         ) != vph_mi.end();

            //         PandoraKeyholeEnergy += ph_ev->ROISummedADC();
            //         PandoraKeyholeHits.push_back(hit);

            //         if (tp) PandoraKeyholeEnergyTP += ph_ev->ROISummedADC();

            //         if (cosa < cos(30.F * TMath::DegToRad())) continue;

            //         PandoraConeEnergy += ph_ev->ROISummedADC();
            //         PandoraConeHits.push_back(hit);

            //         if (tp) PandoraConeEnergyTP += ph_ev->ROISummedADC();
            //     }
            // }
        }


        // LOG(!muEndSecHitdQds.empty());
        // if (!inKeepAll && muEndSecHitdQds.empty()) continue;

        // MIPdQds = 0;
        // if (!inBragg || std::distance(sh_mu.endsec_it(), sh_mu.vph.end()) < 2 * inBraggN) {
        //     BraggError = true;
        // } else {
        //     BraggError = false;

        //     VecPtrHit vph_mu_bragg(sh_mu.endsec_it(), sh_mu.vph.end() - 2*inBraggN);
        //     VecPtrHit vph_mu_tail(sh_mu.vph.end() - 2*inBraggN, sh_mu.vph.end());

        //     int n = inBraggN;
        //     while (n--) {
        //         float min_dist = std::numeric_limits<float>::max();
        //         PtrHit closest_hit;
        //         for (PtrHit const& ph_ev : vph_ev_endsec) {
        //             // not already in muon
        //             if (std::find_if(
        //                 vph_mu_bragg.begin(), vph_mu_bragg.end(),
        //                 [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
        //             ) != vph_mu_bragg.end()) continue;
        //             if (std::find_if(
        //                 vph_mu_tail.begin(), vph_mu_tail.end(),
        //                 [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
        //             ) != vph_mu_tail.end()) continue;

        //             float dist = GetDistance(ph_ev, vph_mu_tail.back());
        //             if (dist < min_dist) {
        //                 min_dist = dist;
        //                 closest_hit = ph_ev;
        //             }
        //         }
        //         if (closest_hit) vph_mu_tail.push_back(closest_hit);
        //         else break;
        //     }

            // unsigned i_dQds_max;
            // std::vector<float> bragg_dQds = GetdQds(vph_mu_tail, inRegN, &i_dQds_max);

            // BraggdQds = bragg_dQds[i_dQds_max];
            // for (unsigned i=0; i<=i_dQds_max; i++) {
            //     vph_mu_bragg.push_back(vph_mu_tail[i]);
            // }

            // VecPtrHit::iterator it_dQ_max = std::max_element(
            //     vph_mu_tail.begin(), vph_mu_tail.end(),
            //     [&](PtrHit const& a, PtrHit const& b) -> bool {
            //         return a->ROISummedADC() < b->ROISummedADC();
            //     }
            // );

            // MIPdQds = std::accumulate(
            //     sh_mu.endsec_it(), sh_mu.vph.end() - inBraggN, 0.F,
            //     [&](float sum, PtrHit const& ph) -> float {
            //         return sum + ph->ROISummedADC();
            //     }
            // ) / (sh_mu.vph.end() - inBraggN - sh_mu.endsec_it());
            // BraggdQds = (*it_dQ_max)->ROISummedADC();
            // for (VecPtrHit::iterator it=vph_mu_tail.begin(); it!=it_dQ_max+1; ++it)
            //     vph_mu_bragg.push_back(*it);



            // BraggEndHit = GetHit(vph_mu_bragg.back());
            // for (PtrHit const& ph_mu : vph_mu_bragg)
            //     BraggMuonHits.push_back(GetHit(ph_mu));

            // BraggBaryHasLongTrack = false;
            // for (PtrHit const& ph_ev : vph_ev_endsec) {
            //     double dist = GetDistance(ph_ev, vph_mu_bragg.back());
            //     if (dist > inMichelRadius) continue;

            //     PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
            //     // not from other long track
            //     if (pt_hit
            //         && pt_hit.key() != pt_ev.key()
            //         && pt_hit->Length() > inMichelRadius
            //     ) {
            //         if (dist < inBarycenterRadius) BraggBaryHasLongTrack = true;
            //         continue;
            //     }

            //     if (std::find_if(
            //         vph_mu_bragg.begin(), vph_mu_bragg.end(),
            //         [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
            //     ) != vph_mu_bragg.end()) continue;

            //     if (dist < inBarycenterRadius) {
            //         BraggBaryHits.push_back(GetHit(ph_ev));
            //     }

            //     // PtrShw ps_hit = fop_hit2shw.at(ph_ev.key());
            //     // BraggSphereHasShower = bool(ps_hit);

            //     ana::Hit hit = GetHit(ph_ev);
            //     BraggSphereEnergy += ph_ev->ROISummedADC();
            //     BraggSphereHits.push_back(hit);

            //     float da = (hit.vec(fTick2cm) - TrkEndHit.vec(fTick2cm)).angle() - TrkReg.theta(dirZ);
            //     da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
            //     BraggSphereHitMuonAngle.push_back(da);

            //     if (std::find_if(
            //         vph_mi.begin(), vph_mi.end(),
            //         [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
            //     ) != vph_mi.end())
            //         BraggSphereEnergyTP = ph_ev->ROISummedADC();
            // }

            // LOG(BraggBaryHits.size());
            // if (BraggBaryHits.size()) {
            //     BraggBary = BraggBaryHits.barycenter(fTick2cm);
            //     ana::Vec2 end_bary = BraggBary - BraggEndHit.vec(fTick2cm);
            //     BraggBaryAngle = end_bary.angle();
            //     float da = BraggBaryAngle - TrkReg.theta(dirZ);
            //     da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
            //     BraggBaryMuonAngle = da;
            // }

            // if (inCone) {

            // }
            // // Cone
            // for (PtrHit const& ph_near : vph_near) {
            //     if (GetDistance(ph_near, bragg.end) > 10) continue;
            //     NearbyBaryHits.push_back(GetHit(ph_near));
            // }

            // LOG(NearbyBaryHits.size());
            // if (NearbyBaryHits.size()) {
            //     NearbyBary = NearbyBaryHits.barycenter(fTick2cm);
            //     ana::Vec2 end_bary = NearbyBary - BraggEndHit.vec(fTick2cm);
            //     NearbyBaryAngle = end_bary.angle();
            //     float da = NearbyBaryAngle - TrkReg.theta(dirZ);
            //     da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
            //     NearbyBaryMuonAngle = da;

            //     // float angle = end_bary.angle();
            //     BraggConeEnergy = 0;
            //     BraggConeEnergyTP = 0;
            //     for (PtrHit const& ph_near : vph_near) {
            //         float dist = GetDistance(ph_near, bragg.end);
            //         if (dist > 30) continue;

            //         ana::Hit hit = GetHit(ph_near);
            //         ana::Vec2 end_hit = hit.vec(fTick2cm) - BraggEndHit.vec(fTick2cm);
            //         float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

            //         if (dist > 5
            //             && cosa < cos(30.F * TMath::DegToRad())
            //         ) continue;

            //         bool tp = std::find_if(
            //             vph_mi.begin(), vph_mi.end(),
            //             [&ph_near](PtrHit const& h) -> bool { return h.key() == ph_near.key(); }
            //         ) != vph_mi.end();

            //         BraggKeyholeEnergy += ph_near->ROISummedADC();
            //         BraggKeyholeHits.push_back(hit);

            //         if (tp) BraggKeyholeEnergyTP += ph_near->ROISummedADC();

            //         if (cosa < cos(30.F * TMath::DegToRad())) continue;

            //         BraggConeEnergy += ph_near->ROISummedADC();
            //         BraggConeHits.push_back(hit);

            //         if (tp) BraggConeEnergyTP += ph_near->ROISummedADC();
            //     }
            // }
        // }

        // Truth Information
        LOG(mcp);
        if (mcp) {
            truPdg = mcp->PdgCode();
            truEndProcess = mcp->EndProcess();
            truStartPoint = ana::Point(mcp->Position().Vect());
            truEndPoint = ana::Point(mcp->EndPosition().Vect());
            truEndEnergy = (mcp->EndE() - mcp->Mass()) * 1e3; // MeV

            LOG(mcp_mi);
            if (mcp_mi) {
                truHasMichel = (
                    geoTop.isInside(mcp_mi->Position().Vect(), 20.F)
                    || geoBot.isInside(mcp_mi->Position().Vect(), 20.F)
                ) ? kHasMichelFiducial : (
                    geoTop.isInside(mcp_mi->Position().Vect())
                    || geoBot.isInside(mcp_mi->EndPosition().Vect())
                    ? kHasMichelInside
                    : kHasMichelOutside
                );
                miTrueEnergy = (mcp_mi->E() - mcp_mi->Mass()) * 1e3;
                PtrTrk pt_mi = ana::mcp2trk(mcp_mi, vpt_ev, clockData, fmp_trk2hit);
                miTrackLength = pt_mi ? pt_mi->Length() : util::kBogusF;
                // PtrShw ps_mi = ana::mcp2shw(mcp_mi, vps_ev, clockData, fmp_shw2hit);
                // MichelShowerLength = ps_mi ? ps_mi->Length() : util::kBogusF;

                for (size_t i=0; i<vph_mi.size(); i++) {
                    PtrHit const& ph_mi = vph_mi[i];
                    float energyFrac = energyFracs_mi[i];
                    if (ph_mi->View() != geo::kW) continue;
                    miHits.push_back(GetHit(ph_mi));
                    miHitEnergyFrac.push_back(energyFrac);
                }
                miHitEnergy = std::accumulate(miHits.adc.begin(), miHits.adc.end(), 0.F);
            }

            // ana::SortedHits sh_mcp = GetSortedHits(vph_mcp_mu, mcp->EndX() > mcp->Vx() ? 1 : -1);
            // LOG(sh_mcp);
            // if (sh_mcp) {

            //     truStartHit = GetHit(sh_mcp.start());
            //     truEndHit = GetHit(sh_mcp.end());
            //     truReg = sh_mcp.end_reg(geoDet);

            //     LOG(mcp_mi);
            //     if (mcp_mi) {
            //         truHasMichel = (
            //             geoTop.isInside(mcp_mi->Position().Vect(), 20.F)
            //             || geoBot.isInside(mcp_mi->Position().Vect(), 20.F)
            //         ) ? kHasMichelFiducial : (
            //             geoTop.isInside(mcp_mi->Position().Vect())
            //             || geoBot.isInside(mcp_mi->EndPosition().Vect())
            //             ? kHasMichelInside
            //             : kHasMichelOutside
            //         );
            //         miTrueEnergy = (mcp_mi->E() - mcp_mi->Mass()) * 1e3;


            //         PtrTrk pt_mi = ana::mcp2trk(mcp_mi, vpt_ev, clockData, fmp_trk2hit);
            //         miTrackLength = pt_mi ? pt_mi->Length() : util::kBogusF;
            //         // PtrShw ps_mi = ana::mcp2shw(mcp_mi, vps_ev, clockData, fmp_shw2hit);
            //         // MichelShowerLength = ps_mi ? ps_mi->Length() : util::kBogusF;

            //         truEndAngle = sh_mcp.end_reg(geoDet).theta(mcp->EndZ() > mcp->Vz() ? 1 : -1);
            //         Hits bary_hits;
            //         for (size_t i=0; i<vph_mi.size(); i++) {
            //             PtrHit const& ph_mi = vph_mi[i];
            //             float energyFrac = energyFracs_mi[i];

            //             if (ph_mi->View() != geo::kW) continue;
            //             Hit hit = GetHit(ph_mi);
            //             miHits.push_back(hit);
            //             miHitEnergyFrac.push_back(energyFrac);

            //             if (hit.section != truEndHit.section) {
            //                 miHitMuonAngle.push_back(100);
            //             } else {
            //                 float da = (hit.vec(fTick2cm) - truEndHit.vec(fTick2cm)).angle() - truEndAngle;
            //                 da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
            //                 miHitMuonAngle.push_back(da);
            //             }

            //             if (GetDistance(ph_mi, sh_mcp.end()) > inBarycenterRadius) continue;
            //             bary_hits.push_back(GetHit(ph_mi));
            //         }
            //         miHitEnergy = miHits.energy();

            //         LOG(miBaryNHit);
            //         if (bary_hits.size()) {
            //             miBary = bary_hits.barycenter(fTick2cm);
            //             ana::Vec2 end_bary = miBary - truEndHit.vec(fTick2cm);
            //             miBaryAngle = end_bary.angle();
            //             float da = miBaryAngle - truEndAngle;
            //             da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
            //             miBaryMuonAngle = da;
            //         }
            //     }
            // }
        }

        muTree->Fill();
        evMuonIndices.push_back(muIndex);
        muIndex++;
        evMuonNumber++;
    } // end of loop over tracks
    evTree->Fill();
    evIndex++;
}

void ana::MichelAnalysis::beginJob() {}
void ana::MichelAnalysis::endJob() {}

void ana::MichelAnalysis::resetEvent() {
    evMuonNumber = 0;
    evMuonIndices.clear();
    evHits.clear();
}
void ana::MichelAnalysis::resetMuon() {
    muHits.clear();
    muHitCathodeX.clear();
    muHitAnodeX.clear();
    muHitY.clear();
    muStartHit = ana::Hit{};
    muStartHitCathodeX = util::kBogusF;
    muStartHitY = util::kBogusF;
    muStartInCathodeX = false;
    muStartInY = false;
    muStartInZ = false;
    muStartInT = false;
    muStartInYZT = false;
    muEndHit = ana::Hit{};
    muEndHitCathodeX = util::kBogusF;
    muEndHitY = util::kBogusF;
    muEndInCathodeX = false;
    muEndInY = false;
    muEndInZ = false;
    muEndInT = false;
    muEndInYZT = false;
    // muTopReg = ana::LinearRegression{};
    // muBotReg = ana::LinearRegression{};
    muEndAngle = util::kBogusF;
    muCathodeCrossing = false;
    muCathodeMisalignment = util::kBogusF;
    muAnodeCrossing = false;
    muHitdQds.clear();
    muSphereHits.clear();
    muSphereEnergy = util::kBogusF;
    muSphereEnergyTP = util::kBogusF;
    muSphereHitMuonAngle.clear();
    // muSphereMaxShowerEnergy = util::kBogusF;

    muBary = ana::Vec2{0,0};
    muBaryHits.clear();
    muBaryAngle = util::kBogusF;
    muBaryMuonAngle = util::kBogusF;
    muSphereHasLongTrack = false;

    truPdg = 0;
    truEndProcess = "";
    truStartPoint = ana::Point{};
    truEndPoint = ana::Point{};
    truEndEnergy = util::kBogusF;
    // truStartHit = ana::Hit{};
    // truEndHit = ana::Hit{};
    // truReg = ana::LinearRegression{};
    // truEndAngle = util::kBogusF;
    truHasMichel = kHasNoMichel;

    miTrueEnergy = util::kBogusF;
    miTrackLength = util::kBogusF;
    // miShowerLength = util::kBogusF;
    miHits.clear();
    miHitEnergyFrac.clear();
    // miHitMuonAngle.clear();
    miHitEnergy = util::kBogusF;

    // miBaryNHit = 0;
    // miBary = ana::Vec2{0,0};
    // miBaryAngle = util::kBogusF;
    // miBaryMuonAngle = util::kBogusF;
}

DEFINE_ART_MODULE(ana::MichelAnalysis)