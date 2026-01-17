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

Top means X>0, Bot means X<0 (by reference to PDVD)
*/

namespace ana { class MichelAnalysis; }

class ana::MichelAnalysis: 
    public art::EDAnalyzer, 
    private ana::MichelAnalyzer 
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
    ana::Bounds<float> wireWindow;
    ana::Bounds3D<float> geoTop, geoBot;
    float geoCathodeGap; // cm

    // Input Parameters
    bool        fLog;
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

    // Track information
    unsigned                muIndex=0;
    float                   muLength;
    ana::Point              muStartPoint;
    ana::Point              muEndPoint;

    // Track information: Hits
    ana::Hits               muHits;
    ana::Hit                muStartHit; 
    float                   muStartHitX;
    float                   muStartHitY;
    bool                    muStartInXYZT;    
    ana::Hit                muEndHit;
    float                   muEndHitX; 
    float                   muEndHitY;
    bool                    muEndInXYZT;
    bool                    muRegError;
    ana::LinearRegression   muTopReg;
    ana::LinearRegression   muBotReg;
    float                   muEndAngle;
    bool                    muCathodeCrossing;
    bool                    muCathodeMisaligned;
    bool                    muAnodeCrossing;
    std::vector<float>      muHitdQds;
    ana::Hits               muSphereHits;
    float                   muSphereEnergy;
    float                   muSphereEnergyTP;
    std::vector<float>      muSphereHitMuonAngle;
    // float                   muSphereMaxShowerEnergy;

    // Track information: Hits nearby end
    ana::Vec2   muBary;
    ana::Hits   muBaryHits;
    float       muBaryAngle;
    float       muBaryMuonAngle;
    bool        muSphereHasLongTrack;

    // Truth information
    int                     truPdg;
    std::string             truEndProcess;
    ana::Point              truStartPoint;
    ana::Point              truEndPoint;
    float                   truEndEnergy;
    ana::Hit                truStartHit;
    ana::Hit                truEndHit;
    ana::LinearRegression   truReg;
    float                   truEndAngle;
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
    std::vector<float>  miHitMuonAngle;
    float               miHitEnergy;

    // Truth information: Hits nearby muon's end
    unsigned            miBaryNHit;
    ana::Vec2           miBary;
    float               miBaryAngle;
    float               miBaryMuonAngle;

    void resetEvent(void);
    void resetMuon(void);
};

ana::MichelAnalysis::MichelAnalysis(fhicl::ParameterSet const& p) : 
    EDAnalyzer{p}, 
    MichelAnalyzer{p},
    fLog(p.get<bool>("Log", true)),
    inKeepAll(p.get<bool>("KeepAll", true)),
    inTrackLengthCut(p.get<float>("TrackLengthCut", 30.F)), // in cm
    inFiducialLength(p.get<float>("FiducialLength", 20.F)), // in cm
    inBarycenterRadius(p.get<float>("BarycenterRadius", 10.F)), // in cm
    inMichelRadius(p.get<float>("MichelRadius", 20.F)), //in cm
    inRegN(p.get<unsigned>("RegN", 6))
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
        default: break;
    }
    geoCathodeGap = geoTop.x.min - geoBot.x.max;

    std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Detector Geometry: " << asGeo->DetectorName()
        << "  (" << (!geoDet ? "PDVD" : "PDHD") << ")" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  Top Bounds: " << geoTop << std::endl
        << "  Bot Bounds: " << geoBot << std::endl;
    std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << inTrackLengthCut << " cm" << std::endl
        << "  Fiducial Length: " << inFiducialLength << " cm" << std::endl
        << "  Barycenter Radius: " << inBarycenterRadius << " cm" << std::endl
        << "  Michel Space Radius: " << inMichelRadius << " cm" << std::endl;

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
    muTree->Branch("RegError",                  &muRegError);
    muTree->Branch("CathodeCrossing",           &muCathodeCrossing);
    muTree->Branch("AnodeCrossing",             &muAnodeCrossing);
    SetBranches(muTree, "Start",                &muStartHit);
    muTree->Branch("StartHitX",                 &muStartHitX);
    muTree->Branch("StartHitY",                 &muStartHitY);
    muTree->Branch("StartInXYZT",               &muStartInXYZT);
    SetBranches(muTree, "End",                  &muEndHit);
    muTree->Branch("EndHitX",                   &muEndHitX);
    muTree->Branch("EndHitY",                   &muEndHitY);
    muTree->Branch("EndInXYZT",                 &muEndInXYZT);
    muTree->Branch("EndAngle",                  &muEndAngle);
    SetBranches(muTree, "",                     &muHits);
    SetBranches(muTree, "Top",                  &muTopReg);
    SetBranches(muTree, "Bot",                  &muBotReg);
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
    SetBranches(muTree, "TrueStart",        &truStartHit);
    SetBranches(muTree, "TrueEnd",          &truEndHit);
    SetBranches(muTree, "True",             &truReg);
    muTree->Branch("TrueEndAngle",          &truEndAngle);
    muTree->Branch("TrueHasMichel",   (int*)&truHasMichel);
    muTree->Branch("MichelTrueEnergy",      &miTrueEnergy); // MeV
    muTree->Branch("MichelTrackLength",     &miTrackLength); // cm
    // muTree->Branch("MichelShowerLength",    &miShowerLength); // cm
    SetBranches(muTree, "Michel",           &miHits);
    muTree->Branch("MichelHitEnergyFrac",   &miHitEnergyFrac);
    muTree->Branch("MichelHitMuonAngle",    &miHitMuonAngle);
    muTree->Branch("MichelHitEnergy",       &miHitEnergy); // ADC

    muTree->Branch("MichelBaryNHit",        &miBaryNHit);
    SetBranches(muTree, "MichelBary",       &miBary); 
    muTree->Branch("MichelBaryAngle",       &miBaryAngle); // rad
    muTree->Branch("MichelBaryMuonAngle",   &miBaryMuonAngle); // rad
}

void ana::MichelAnalysis::analyze(art::Event const& e) {
    auto const clockData = asDetClocks->DataFor(e);
    auto const detProp = asDetProp->DataFor(e,clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();

    auto const & vh_hit = e.getHandle<std::vector<recob::Hit>>(tag_hit);
    if (!vh_hit.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Hit handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrHit vph_ev;
    art::fill_ptr_vector(vph_ev, vh_hit);

    auto const & vh_trk = e.getHandle<std::vector<recob::Track>>(tag_trk);
    if (!vh_trk.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Track handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrTrk vpt_ev;
    art::fill_ptr_vector(vpt_ev, vh_trk);

    // auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
    // if (!vh_shw.isValid()) {
    //     std::cout << "\033[1;91m" "No valid recob::Shower handle" "\033[0m" << std::endl;
    //     return;
    // }
    // VecPtrShw vps_ev;
    // art::fill_ptr_vector(vps_ev, vh_shw);

    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
    // art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
    // art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

    resetEvent();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();
    evIsData = e.isRealData();

    for (PtrHit p_hit : vph_ev)
        if (p_hit->View() == geo::kW)
            evHits.push_back(GetHit(p_hit));

    // loop over tracks to find stopping muons
    for (PtrTrk const& pt_ev : vpt_ev) {
        if (fLog) std::cout << "e" << evIndex << "t" << pt_ev->ID() << "\r" << std::flush;
        resetMuon();

        VecPtrHit vph_mu = fmp_trk2hit.at(pt_ev.key());
        std::vector<recob::TrackHitMeta const*> const& vhm_mu = fmp_trk2hit.data(pt_ev.key());
        std::map<size_t, unsigned> map_hitkey2metaidx;
        ASSERT(vph_mu.size())

        ASSERT(vph_mu.size() == vhm_mu.size())
        std::vector<unsigned> bad_hit_indices;
        for (unsigned i=0; i<vph_mu.size(); i++) {
            if (vph_mu[i]->View() != geo::kW) continue;
            if (!pt_ev->HasValidPoint(vhm_mu[i]->Index())) {
                bad_hit_indices.push_back(i);
            } else {
                map_hitkey2metaidx[vph_mu[i].key()] = i;
            }
        }
        for (int i=bad_hit_indices.size()-1; i>=0; i--)
            vph_mu.erase(vph_mu.begin() + bad_hit_indices[i]);


        simb::MCParticle const* mcp = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        simb::MCParticle const* mcp_mi = nullptr;
        VecPtrHit vph_mcp_mu, vph_mi;
        std::vector<float> energyFracs_mi;
        if (mcp) {
            mcp_mi = GetMichelMCP(mcp);
            vph_mcp_mu = ana::mcp2hits(mcp, vph_ev, clockData, false);
            vph_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, true, &energyFracs_mi);
        }

        if (fLog) std::cout << "\t" "\033[1;93m" "e" << evIndex << "m" << evMuonNumber << " (" << muIndex << ")" "\033[0m" << std::endl;

        // ============================
        // Dump basic track information
        muLength = pt_ev->Length();
        // muChi2 = pt_ev->Chi2();
        // muChi2PerNdof = pt_ev->Chi2PerNdof();

        LOG(muLength >= inTrackLengthCut);
        if (!inKeepAll && muLength < inTrackLengthCut) continue;

        bool is_up =  IsUpright(*pt_ev);
        geo::Point_t Start = is_up ? pt_ev->Start() : pt_ev->End();
        geo::Point_t End = is_up ? pt_ev->End() : pt_ev->Start();
        muStartPoint = ana::Point(Start);
        muEndPoint = ana::Point(End);

        // ===================================================
        // Sorting hits according to give X direction
        // in PDVD that assums downward muons
        //
        // in case of hits on both sides (cathode crossing): 
        // sh_mu.cc.first -> last hit in the top volume (side1)
        // sh_mu.cc.second -> first hit in the bot volume (side0)

        ana::SortedHits sh_mu = geoDet == kPDVD
            ? GetSortedHits_dirX(vph_mu, -1) // PDVD: decreasing X <-> downward
            : GetSortedHits_dirX(vph_mu, End.X() > Start.X() ? 1 : -1); // PDHD
        muRegError = !sh_mu;


        // Sort Hits by their Index in the TrackHitMeta data
        // Never tested
        /*
        ana::SortedHits sh_mu;
        sh_mu.vph = vph_mu; // collection hits with valid points
        std::sort(sh_mu.vph.begin(), sh_mu.vph.end(),
            [&map_hitkey2metaidx](PtrHit const& a, PtrHit const& b) {
                return map_hitkey2metaidx.at(a.key()) < map_hitkey2metaidx.at(b.key());
            }
        );
        int prev_sec = -1, prev_side = -1;
        PtrHit& prev_ph = sh_mu.vph.front();
        for (PtrHit const& ph : sh_mu.vph) {
            int sec = ana::tpc2sec.at(geoDet).at(ph->WireID().TPC);
            int side = ana::tpc2side.at(geoDet).at(ph->WireID().TPC);
            if (prev_sec == -1 || sec != prev_sec) {
                sh_mu.secs.push_back(sec);
                sh_mu.sc.push_back(prev_ph);
                sh_mu.sc.push_back(ph);
            }
            prev_sec = sec;
            if (side != -1) {
                double z = GetSpace(ph->WireID());
                double t = ph->PeakTime() * fTick2cm;
                sh_mu.regs[side].add(z, t);
            }
            if (prev_side != -1 && side != prev_side) {
                sh_mu.cc = std::make_pair(prev_ph, ph);
            }
            prev_side = side;
            prev_ph = ph;
            // sh_mu.bot_index?
            // sh_mu.endsec_index?
        }
        sh_mu.regs[0].compute();
        sh_mu.regs[1].compute();

        sh_mu.start = sh_mu.vph.front();
        sh_mu.end = sh_mu.vph.back();
        */


        LOG(!muRegError);
        if (!inKeepAll && muRegError) continue;
        if (!muRegError) {
            // First hit of the track:
            muStartHit = GetHit(sh_mu.start());
            size_t start_track_idx = vhm_mu[map_hitkey2metaidx.at(sh_mu.start().key())]->Index();
            muStartHitY = pt_ev->HasValidPoint(start_track_idx)
                ? pt_ev->LocationAtPoint(start_track_idx).Y()
                : util::kBogusF;
            bool start_in_Y = geoBot.y.isInside(muStartHitY, inFiducialLength) || geoTop.y.isInside(muStartHitY, inFiducialLength);
            bool start_in_Z = geoBot.z.isInside(muStartHit.space, inFiducialLength) || geoTop.z.isInside(muStartHit.space, inFiducialLength);
            bool start_in_T = wireWindow.isInside(muStartHit.tick, inFiducialLength / fTick2cm);

            // Last hit of the track:
            muEndHit = GetHit(sh_mu.end());
            size_t end_track_idx = vhm_mu[map_hitkey2metaidx.at(sh_mu.end().key())]->Index();
            muEndHitY = pt_ev->HasValidPoint(end_track_idx)
                ? pt_ev->LocationAtPoint(end_track_idx).Y()
                : util::kBogusF;
            bool end_in_Y = geoBot.y.isInside(muEndHitY, inFiducialLength) || geoTop.y.isInside(muEndHitY, inFiducialLength);
            bool end_in_Z = geoBot.z.isInside(muEndHit.space, inFiducialLength) || geoTop.z.isInside(muEndHit.space, inFiducialLength);
            bool end_in_T = wireWindow.isInside(muEndHit.tick, inFiducialLength / fTick2cm);

            // Linear regression of hits, same side of the cathode as the end hit
            muTopReg = sh_mu.regs[1];
            muBotReg = sh_mu.regs[0];
            LOG(muTopReg.r2 >= 0.5 && muBotReg.r2 >= 0.5);
            if (!inKeepAll && !(muTopReg.r2 >= 0.5 && muBotReg.r2 >= 0.5)) continue;

            // muEndSecHitdQds = GetdQds(sh_mu.endsec_it(), sh_mu.vph.end(), inRegN);
            muHitdQds.clear();
            std::vector<float> tmp;
            tmp = GetdQds(sh_mu.vph.begin(), sh_mu.after_cathode_it(), inRegN);
            muHitdQds.insert(muHitdQds.end(), tmp.begin(), tmp.end());
            tmp = GetdQds(sh_mu.after_cathode_it(), sh_mu.vph.end(), inRegN);
            muHitdQds.insert(muHitdQds.end(), tmp.begin(), tmp.end());
            tmp.clear();

            // Cathode crossing <-> track has hits on both sides of the cathode
            muCathodeCrossing = sh_mu.is_cc();
            muCathodeMisaligned = sh_mu.is_cc()
                && abs(sh_mu.cc_first()->PeakTime()-sh_mu.cc_second()->PeakTime())*fTick2cm < 3 * geoCathodeGap;

            LOG(muCathodeCrossing);
            if (!inKeepAll && !muCathodeCrossing) continue;

            // switch (geoDet) {
            // case kPDVD: /* ASSUMS DOWNWARD MUON */
            //     muAnodeCrossing = muStartHit.section < 4
            //         && geoTop.z.isInside(muStartHit.space, inFiducialLength)
            //         && wireWindow.isInside(muStartHit.tick, inFiducialLength/fTick2cm);
            //     break;
            // case kPDHD:
            //     muAnodeCrossing =
            //         geoTop.z.isInside(muStartHit.space, inFiducialLength)
            //         && wireWindow.isInside(muStartHit.tick, inFiducialLength/fTick2cm);
            //     break;
            // }
            // LOG(muCathodeCrossing || muAnodeCrossing);
            // if (!inKeepAll && (!muAnodeCrossing && !muCathodeCrossing)) continue;


            // Calculate hit X positions for cathode crossing tracks
            bool start_in_X=false;
            bool end_in_X=false;
            if (muCathodeCrossing) {
                int start_side = ana::tpc2side.at(geoDet).at(muStartHit.tpc);
                int end_side = ana::tpc2side.at(geoDet).at(muEndHit.tpc);

                // side0: X<0 | side1: X>0
                muStartHitX = start_side == 0
                    ? -(geoCathodeGap/2) - (sh_mu.cc_second()->PeakTime() - muStartHit.tick) * fTick2cm
                    : +(geoCathodeGap/2) + (sh_mu.cc_first()->PeakTime() - muStartHit.tick) * fTick2cm;
                start_in_X = start_side == 0
                    ? geoBot.x.isInside(muStartHitX, inFiducialLength)
                    : geoTop.x.isInside(muStartHitX, inFiducialLength);

                muEndHitX = end_side == 0
                    ? -(geoCathodeGap/2) - (sh_mu.cc_second()->PeakTime() - muEndHit.tick) * fTick2cm
                    : +(geoCathodeGap/2) + (sh_mu.cc_first()->PeakTime() - muEndHit.tick) * fTick2cm;
                end_in_X = end_side == 0
                    ? geoBot.x.isInside(muEndHitX, inFiducialLength)
                    : geoTop.x.isInside(muEndHitX, inFiducialLength);
            }

            muStartInXYZT = start_in_X && start_in_Y && start_in_Z && start_in_T;
            muEndInXYZT = end_in_X && end_in_Y && end_in_Z && end_in_T;
            LOG(muStartInXYZT);
            LOG(muEndInXYZT);
            if (!inKeepAll && !muEndInXYZT) continue;
            for (PtrHit const& ph_mu : sh_mu.vph) {
                if (ph_mu->View() != geo::kW) continue;
                ana::Hit hit = GetHit(ph_mu);
                muHits.push_back(hit);
            }

            VecPtrHit vph_ev_endsec;
            for (PtrHit const& ph_ev : vph_ev) {
                if (ph_ev->View() != geo::kW) continue;
                int sec = ana::tpc2sec.at(geoDet).at(ph_ev->WireID().TPC);
                if (sec != sh_mu.end_sec()) continue;
                vph_ev_endsec.push_back(ph_ev);
            }

            muEndAngle = sh_mu.regs
                [ana::sec2side.at(geoDet).at(muEndHit.section)]
                .theta(muEndHit.space > muStartHit.space ? 1 : -1);
            // integrate charges around muon endpoint
            muSphereEnergy = 0;
            muSphereEnergyTP = 0;
            muSphereHasLongTrack = false;
            // muSphereMaxShowerEnergy = 0;
            for (PtrHit const& ph_ev : vph_ev_endsec) {
                float dist = GetDistance(ph_ev, sh_mu.end());
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

                float da = (hit.vec(fTick2cm) - muEndHit.vec(fTick2cm)).angle() - muEndAngle;
                da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                muSphereHitMuonAngle.push_back(da);

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
        }

        // Truth Information
        LOG(mcp);
        if (mcp) {
            truPdg = mcp->PdgCode();
            truEndProcess = mcp->EndProcess();
            truStartPoint = ana::Point(mcp->Position().Vect());
            truEndPoint = ana::Point(mcp->EndPosition().Vect());
            truEndEnergy = (mcp->EndE() - mcp->Mass()) * 1e3; // MeV

            ana::SortedHits sh_mcp = GetSortedHits_dirX(vph_mcp_mu, mcp->EndX() > mcp->Vx() ? 1 : -1);
            LOG(sh_mcp);
            if (sh_mcp) {
                truStartHit = GetHit(sh_mcp.start());
                truEndHit = GetHit(sh_mcp.end());
                truReg = sh_mcp.end_reg(geoDet);

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
                    miTrackLength = pt_mi ? pt_mi->Length() : -1.F;
                    // PtrShw ps_mi = ana::mcp2shw(mcp_mi, vps_ev, clockData, fmp_shw2hit);
                    // MichelShowerLength = ps_mi ? ps_mi->Length() : -1.F;

                    truEndAngle = sh_mcp.end_reg(geoDet).theta(mcp->EndZ() > mcp->Vz() ? 1 : -1);
                    Hits bary_hits;
                    for (size_t i=0; i<vph_mi.size(); i++) {
                        PtrHit const& ph_mi = vph_mi[i];
                        float energyFrac = energyFracs_mi[i];

                        if (ph_mi->View() != geo::kW) continue;
                        Hit hit = GetHit(ph_mi);
                        miHits.push_back(hit);
                        miHitEnergyFrac.push_back(energyFrac);

                        if (hit.section != truEndHit.section) {
                            miHitMuonAngle.push_back(100);
                        } else {
                            float da = (hit.vec(fTick2cm) - truEndHit.vec(fTick2cm)).angle() - truEndAngle;
                            da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                            miHitMuonAngle.push_back(da);
                        }

                        if (GetDistance(ph_mi, sh_mcp.end()) > inBarycenterRadius) continue;
                        bary_hits.push_back(GetHit(ph_mi));
                    }
                    miHitEnergy = miHits.energy();

                    LOG(miBaryNHit);
                    if (bary_hits.size()) {
                        miBary = bary_hits.barycenter(fTick2cm);
                        ana::Vec2 end_bary = miBary - truEndHit.vec(fTick2cm);
                        miBaryAngle = end_bary.angle();
                        float da = miBaryAngle - truEndAngle;
                        da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                        miBaryMuonAngle = da;
                    }
                }
            }
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
    muStartHit = ana::Hit{};
    muStartHitX = util::kBogusF;
    muStartHitY = util::kBogusF;
    muEndHit = ana::Hit{};
    muEndHitX = util::kBogusF;
    muEndHitY = util::kBogusF;
    muTopReg = ana::LinearRegression{};
    muBotReg = ana::LinearRegression{};
    muEndAngle = util::kBogusF;
    muCathodeCrossing = false;
    muCathodeMisaligned = false;
    muAnodeCrossing = false;
    muHitdQds.clear();
    muSphereHits.clear();
    muSphereEnergy = -1.F;
    muSphereEnergyTP = -1.F;
    muSphereHitMuonAngle.clear();
    // muSphereMaxShowerEnergy = -1.F;

    muBary = ana::Vec2{0,0};
    muBaryHits.clear();
    muBaryAngle = util::kBogusF;
    muBaryMuonAngle = util::kBogusF;
    muSphereHasLongTrack = false;

    truPdg = 0;
    truEndProcess = "";
    truStartPoint = ana::Point{};
    truEndPoint = ana::Point{};
    truEndEnergy = -1.F;
    truStartHit = ana::Hit{};
    truEndHit = ana::Hit{};
    truReg = ana::LinearRegression{};
    truEndAngle = util::kBogusF;
    truHasMichel = kHasNoMichel;

    miTrueEnergy = -1.F;
    miTrackLength = -1.F;
    // miShowerLength = -1.F;
    miHits.clear();
    miHitEnergyFrac.clear();
    miHitMuonAngle.clear();
    miHitEnergy = -1.F;

    miBaryNHit = 0;
    miBary = ana::Vec2{0,0};
    miBaryAngle = util::kBogusF;
    miBaryMuonAngle = util::kBogusF;
}

DEFINE_ART_MODULE(ana::MichelAnalysis)