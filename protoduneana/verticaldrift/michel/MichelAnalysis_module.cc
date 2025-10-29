#include "utils.h"

namespace ana {
    class MichelAnalysis;
}

// using HitPtr = art::Ptr<recob::Hit>;
// using HitPtrVec = std::vector<art::Ptr<recob::Hit>>;
// using HitPtrPair = std::pair<art::Ptr<recob::Hit>, art::Ptr<recob::Hit>>;

class ana::MichelAnalysis : public art::EDAnalyzer, private ana::MichelAnalyzer {
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
    ana::Bounds3D<float> geoHighX, geoLowX;
    float fCathodeGap; // cm

    // Input Parameters
    bool fLog;
    float fTrackLengthCut; // in cm
    float fFiducialLength; // in cm
    float fBarycenterRadius; // in cm
    float fMichelRadius; // in cm
    float fNearbyRadius; // in cm
    float fBodyDistance; // in cm
    unsigned fRegN;
    float fBraggThreshold; // in MIP

    // Output Variables
    TTree *tEvent, *tMuon;

    // Event information
    unsigned evRun, evSubRun, evEvent;
    bool EventIsReal;
    unsigned iEvent=0;
    unsigned EventNMuon;
    std::vector<unsigned> EventiMuon;
    ana::Hits EventHits;

    // Track information
    unsigned iMuon=0;

    float TrkLength;
    bool TrkIsUpright; // Supposition: Muon is downward
    ana::Point TrkStartPoint, TrkEndPoint;
    bool TrkEndInVolumeYZ;
    bool TrkCathodeCrossing;
    bool TrkAnodeCrossing;
    // bool AgnoTagTrkEndLowN,
    //      AgnoTagTrkEndBadCC;

    // Hit from track information
    bool TrkHitError;
    int TrkHitCathodeCrossing;
    enum EnumCathodeCrossing { kNoCC, kHitOnBothSides, kAlignedHitOnBothSides };
    bool TrkHitAnodeCrossing;
    ana::Hit TrkStartHit, TrkEndHit;
    float TrkEndHitX;
    int TrkRegDirZ;
    ana::LinearRegression TrkReg;
    bool TrkHitEndInVolumeX;
    bool TrkHitEndInWindow;
    ana::Hits TrkHits;
    std::vector<float> PandoraTrkHitdQdx;
    ana::Hits PandoraSphereHits;
    std::vector<float> PandoraSphereHitMuonAngle;
    bool PandoraSphereHasShower;
    bool PandoraBaryHasShower;
    float PandoraSphereEnergy;
    float PandoraSphereEnergyTP;
    // Cone
    ana::Hits PandoraBaryHits;
    ana::Vec2 PandoraBary;
    float PandoraBaryAngle;
    float PandoraBaryMuonAngle;
    ana::Hits PandoraConeHits;
    float PandoraConeEnergy;
    float PandoraConeEnergyTP;
    ana::Hits PandoraKeyholeHits;
    float PandoraKeyholeEnergy;
    float PandoraKeyholeEnergyTP;


    // Bragg information
    int BraggError;
    enum EnumBraggError { kNoError, kEndNotFound, kSmallBody };
    float MIPdQdx;
    float BraggdQdx;
    ana::Hit BraggEndHit;
    ana::Hits BraggMuonHits;

    // Sphere
    ana::Hits BraggSphereHits;
    std::vector<float> BraggSphereHitMuonAngle;
    float BraggSphereEnergy;
    float BraggSphereEnergyTP;

    // Cone
    ana::Hits NearbyBaryHits;
    ana::Vec2 NearbyBary;
    float NearbyBaryAngle;
    float NearbyBaryMuonAngle;
    ana::Hits BraggConeHits;
    float BraggConeEnergy;
    float BraggConeEnergyTP;
    ana::Hits BraggKeyholeHits;
    float BraggKeyholeEnergy;
    float BraggKeyholeEnergyTP;


    // Truth information
    int TruePdg;
    std::string TrueEndProcess;
    ana::Point MuonTrueStartPoint;
    ana::Point MuonTrueEndPoint;
    float MuonTrueEndEnergy;
    bool TrueDownward;
    ana::Hit MuonTrueStartHit;
    ana::Hit MuonTrueEndHit;
    int MuonTrueRegDirZ;
    ana::LinearRegression MuonTrueReg;

    int TrueHasMichel;
    enum EnumHasMichel { kNoMichel, kHasMichelOutside, kHasMichelInside, kHasMichelFiducial };
    float MichelTrueEnergy;
    float MichelTrackLength, MichelShowerLength;
    ana::Hits MichelHits;
    std::vector<float> MichelHitMuonAngle;
    float MichelHitEnergy;

    unsigned MichelBaryNHit;
    ana::Vec2 MichelBary;
    float MichelBaryAngle;
    float MichelBaryMuonAngle;
    float MichelConeEnergy;
    float MichelKeyholeEnergy;

    void resetEvent();
    void resetMuon();

    // HitPtrPair GetTrackEndsHits(
    //     HitPtrVec const& vp_hit,
    //     HitPtrPair *pp_cathode_crossing = nullptr,
    //     HitPtrVec *vp_tpc_crossing = nullptr,
    //     HitPtrVec *vp_sorted_hit = nullptr,
    //     geo::View_t view = geo::kW
    // );
};

ana::MichelAnalysis::MichelAnalysis(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}, MichelAnalyzer{p},
    fLog(p.get<bool>("Log", true)),
    fTrackLengthCut(p.get<float>("TrackLengthCut", 20.F)), // in cm
    fFiducialLength(p.get<float>("FiducialLength", 10.F)), // in cm
    fBarycenterRadius(p.get<float>("BarycenterRadius", 10.F)), // in cm
    fMichelRadius(p.get<float>("MichelRadius", 20.F)), //in cm
    fNearbyRadius(p.get<float>("NearbyRadius", 40.F)), //in cm
    fBodyDistance(p.get<float>("BodyDistance", 20.F)), //in cm
    fRegN(p.get<unsigned>("RegN", 6)),
    fBraggThreshold(p.get<float>("BraggThreshold", 1.7)) // in MIP dE/dx
{
    auto const clockData = asDetClocks->DataForJob();
    auto const detProp = asDetProp->DataForJob(clockData);
    fTick2cm = detinfo::sampling_rate(clockData) * 1e-3 * detProp.DriftVelocity();
    wireWindow = ana::Bounds<float>{0.F, (float) detProp.ReadOutWindowSize()};
    switch (geoDet) {
        case kPDVD:
            geoLowX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 0}).Min(), 
                asGeo->TPC(geo::TPCID{0, 7}).Max()
            };
            geoHighX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 8}).Min(),
                asGeo->TPC(geo::TPCID{0, 15}).Max()
            };
            break;
        case kPDHD:
            geoLowX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 1}).Min(),
                asGeo->TPC(geo::TPCID{0, 5}).Max()
            };
            geoHighX = ana::Bounds3D<float>{
                asGeo->TPC(geo::TPCID{0, 2}).Min(),
                asGeo->TPC(geo::TPCID{0, 6}).Max()
            };
            break;
        default: break;
    }
    fCathodeGap = geoHighX.x.min - geoLowX.x.max;

    std::cout << "\033[1;93m" "Detector Properties:" "\033[0m" << std::endl
        << "  Detector Geometry: " << asGeo->DetectorName()
        << "  (" << (!geoDet ? "PDVD" : "PDHD") << ")" << std::endl
        << "  Tick Window: " << wireWindow << std::endl
        << "  HighX Bounds: " << geoHighX << std::endl
        << "  LowX Bounds: " << geoLowX << std::endl;
    std::cout << "\033[1;93m" "Analysis Parameters:" "\033[0m" << std::endl
        << "  Track Length Cut: " << fTrackLengthCut << " cm" << std::endl
        << "  Fiducial Length: " << fFiducialLength << " cm" << std::endl
        << "  Barycenter Radius: " << fBarycenterRadius << " cm" << std::endl
        << "  Michel Space Radius: " << fMichelRadius << " cm" << std::endl;

    tEvent = asFile->make<TTree>("event","");

    tEvent->Branch("eventRun", &evRun);
    tEvent->Branch("eventSubRun", &evSubRun);
    tEvent->Branch("eventEvent", &evEvent);
    tEvent->Branch("isReal", &EventIsReal);
    tEvent->Branch("iEvent", &iEvent);
    tEvent->Branch("NMuon", &EventNMuon);
    tEvent->Branch("iMuon", &EventiMuon);
    EventHits.SetBranches(tEvent);

    tMuon = asFile->make<TTree>("muon","");

    // Event
    tMuon->Branch("eventRun", &evRun);
    tMuon->Branch("eventSubRun", &evSubRun);
    tMuon->Branch("eventEvent", &evEvent);
    tEvent->Branch("isReal", &EventIsReal);
    tMuon->Branch("iEvent", &iEvent);
    tMuon->Branch("iMuon", &iMuon);
    tMuon->Branch("iMuonInEvent", &EventNMuon);

    // Track
    tMuon->Branch("TrkLength", &TrkLength);
    tMuon->Branch("TrkIsUpright", &TrkIsUpright);
    TrkStartPoint.SetBranches(tMuon, "Start");
    TrkEndPoint.SetBranches(tMuon, "End");
    tMuon->Branch("TrkEndInVolumeYZ", &TrkEndInVolumeYZ);
    tMuon->Branch("TrkCathodeCrossing", &TrkCathodeCrossing);
    tMuon->Branch("TrkAnodeCrossing", &TrkAnodeCrossing);
    // tMuon->Branch("AgnoTagTrkEndLowN", &AgnoTagTrkEndLowN);
    // tMuon->Branch("AgnoTagTrkEndBadCC", &AgnoTagTrkEndBadCC);

    // Hit
    tMuon->Branch("TrkHitError", &TrkHitError);
    tMuon->Branch("TrkHitCathodeCrossing", &TrkHitCathodeCrossing);
    tMuon->Branch("TrkHitAnodeCrossing", &TrkHitAnodeCrossing);
    TrkStartHit.SetBranches(tMuon, "Start");
    TrkEndHit.SetBranches(tMuon, "End");
    tMuon->Branch("EndHitX", &TrkEndHitX);
    tMuon->Branch("RegDirZ", &TrkRegDirZ);
    TrkReg.SetBranches(tMuon, "");
    tMuon->Branch("TrkHitEndInVolumeX", &TrkHitEndInVolumeX);
    tMuon->Branch("TrkHitEndInWindow", &TrkHitEndInWindow);
    TrkHits.SetBranches(tMuon, "");
    tMuon->Branch("HitdQdx", &PandoraTrkHitdQdx);
    PandoraSphereHits.SetBranches(tMuon, "PandoraSphere");
    tMuon->Branch("PandoraSphereHitMuonAngle", &PandoraSphereHitMuonAngle);
    tMuon->Branch("PandoraSphereHasShower", &PandoraSphereHasShower);
    tMuon->Branch("PandoraBaryHasShower", &PandoraBaryHasShower);
    tMuon->Branch("PandoraSphereEnergy", &PandoraSphereEnergy); // ADC
    tMuon->Branch("PandoraSphereEnergyTP", &PandoraSphereEnergyTP); // ADC

    // Cone
    PandoraBaryHits.SetBranches(tMuon, "PandoraBary");
    PandoraBary.SetBranches(tMuon, "PandoraBary");
    tMuon->Branch("PandoraBaryAngle", &PandoraBaryAngle);
    tMuon->Branch("PandoraBaryMuonAngle", &PandoraBaryMuonAngle);
    PandoraConeHits.SetBranches(tMuon, "PandoraCone");
    tMuon->Branch("PandoraConeEnergy", &PandoraConeEnergy); // ADC
    tMuon->Branch("PandoraConeEnergyTP", &PandoraConeEnergyTP); // ADC
    PandoraKeyholeHits.SetBranches(tMuon, "PandoraCone");
    tMuon->Branch("PandoraKeyholeEnergy", &PandoraKeyholeEnergy); // ADC
    tMuon->Branch("PandoraKeyholeEnergyTP", &PandoraKeyholeEnergyTP); // ADC


    // Bragg
    tMuon->Branch("BraggError", &BraggError);
    tMuon->Branch("MIPdQdx", &MIPdQdx);
    tMuon->Branch("BraggdQdx", &BraggdQdx);
    BraggEndHit.SetBranches(tMuon, "BraggEnd");
    BraggMuonHits.SetBranches(tMuon, "BraggMuon");

    // Sphere
    BraggSphereHits.SetBranches(tMuon, "BraggSphere");
    tMuon->Branch("BraggSphereHitMuonAngle", &BraggSphereHitMuonAngle);
    tMuon->Branch("BraggSphereEnergy", &BraggSphereEnergy); // ADC
    tMuon->Branch("BraggSphereEnergyTP", &BraggSphereEnergyTP); // ADC

    // Cone
    NearbyBaryHits.SetBranches(tMuon, "NearbyBary");
    NearbyBary.SetBranches(tMuon, "NearbyBary");
    tMuon->Branch("NearbyBaryAngle", &NearbyBaryAngle);
    tMuon->Branch("NearbyBaryMuonAngle", &NearbyBaryMuonAngle);
    BraggConeHits.SetBranches(tMuon, "BraggCone");
    tMuon->Branch("BraggConeEnergy", &BraggConeEnergy); // ADC
    tMuon->Branch("BraggConeEnergyTP", &BraggConeEnergyTP); // ADC
    BraggKeyholeHits.SetBranches(tMuon, "BraggCone");
    tMuon->Branch("BraggKeyholeEnergy", &BraggKeyholeEnergy); // ADC
    tMuon->Branch("BraggKeyholeEnergyTP", &BraggKeyholeEnergyTP); // ADC

    // Truth
    tMuon->Branch("TruePdg", &TruePdg);
    tMuon->Branch("TrueEndProcess", &TrueEndProcess);
    MuonTrueStartPoint.SetBranches(tMuon, "TrueStart");
    MuonTrueEndPoint.SetBranches(tMuon, "TrueEnd");
    tMuon->Branch("TrueEndEnergy", &MuonTrueEndEnergy); 
    tMuon->Branch("TrueDownward", &TrueDownward);
    MuonTrueStartHit.SetBranches(tMuon, "TrueStart");
    MuonTrueEndHit.SetBranches(tMuon, "TrueEnd");
    tMuon->Branch("TrueRegDirZ", &MuonTrueRegDirZ);
    MuonTrueReg.SetBranches(tMuon, "True");

    tMuon->Branch("TrueHasMichel", &TrueHasMichel);
    tMuon->Branch("MichelTrueEnergy", &MichelTrueEnergy); // MeV
    tMuon->Branch("MichelTrackLength", &MichelTrackLength); // cm
    tMuon->Branch("MichelShowerLength", &MichelShowerLength); // cm
    MichelHits.SetBranches(tMuon, "Michel");
    tMuon->Branch("MichelHitMuonAngle", &MichelHitMuonAngle);
    tMuon->Branch("MichelHitEnergy", &MichelHitEnergy); // ADC

    tMuon->Branch("MichelBaryNHit", &MichelBaryNHit);
    MichelBary.SetBranches(tMuon, "MichelBary");
    tMuon->Branch("MichelBaryAngle", &MichelBaryAngle); // rad
    tMuon->Branch("MichelBaryMuonAngle", &MichelBaryMuonAngle); // rad
    tMuon->Branch("MichelConeEnergy", &MichelConeEnergy); // ADC
    tMuon->Branch("MichelKeyholeEnergy", &MichelKeyholeEnergy); // ADC


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

    auto const & vh_shw = e.getHandle<std::vector<recob::Shower>>(tag_shw);
    if (!vh_shw.isValid()) {
        std::cout << "\033[1;91m" "No valid recob::Shower handle" "\033[0m" << std::endl;
        return;
    }
    VecPtrShw vps_ev;
    art::fill_ptr_vector(vps_ev, vh_shw);

    art::FindManyP<recob::Hit> fmp_trk2hit(vh_trk, e, tag_trk);
    art::FindOneP<recob::Track> fop_hit2trk(vh_hit, e, tag_trk);
    art::FindManyP<recob::Hit> fmp_shw2hit(vh_shw, e, tag_shw);
    art::FindOneP<recob::Shower> fop_hit2shw(vh_hit, e, tag_shw);

    resetEvent();

    evRun = e.run();
    evSubRun = e.subRun();
    evEvent = e.event();
    EventIsReal = e.isRealData();

    for (PtrHit p_hit : vph_ev)
        if (p_hit->View() == geo::kW)
            EventHits.push_back(GetHit(p_hit));

    // loop over tracks to find muons
    for (PtrTrk const& pt_ev : vpt_ev) {
        if (fLog) std::cout << "e" << iEvent << "t" << pt_ev->ID() << "\r" << std::flush;
        resetMuon();

        VecPtrHit vph_mu = fmp_trk2hit.at(pt_ev.key());
        ASSERT(vph_mu.size())

        simb::MCParticle const* mcp = ana::trk2mcp(pt_ev, clockData, fmp_trk2hit);
        simb::MCParticle const* mcp_mi = nullptr;
        VecPtrHit vph_mcp_mu, vph_mi;
        if (mcp) {
            mcp_mi = GetMichelMCP(mcp);
            vph_mcp_mu = ana::mcp2hits(mcp, vph_ev, clockData, false);
            vph_mi = ana::mcp2hits(mcp_mi, vph_ev, clockData, true);
        }

        if (fLog) std::cout << "\t" "\033[1;93m" "e" << iEvent << "m" << EventNMuon << " (" << iMuon << ")" "\033[0m" << std::endl;
        EventiMuon.push_back(iMuon);

        TrkLength = pt_ev->Length();

        TrkIsUpright =  IsUpright(*pt_ev);
        geo::Point_t Start = TrkIsUpright ? pt_ev->Start() : pt_ev->End();
        geo::Point_t End = TrkIsUpright ? pt_ev->End() : pt_ev->Start();
        TrkStartPoint = ana::Point(Start);
        TrkEndPoint = ana::Point(End);

        TrkEndInVolumeYZ = geoHighX.isInsideYZ(End, fFiducialLength);

        TrkCathodeCrossing = (
            geoLowX.isInside(Start)
            && geoHighX.isInside(End)
        ) || (
            geoHighX.isInside(Start)
            && geoLowX.isInside(End)
        );

        // Anode Crossing: SUPPOSITION: Muon is downward
        if (geoDet == kPDVD)
            TrkAnodeCrossing = geoHighX.isInsideYZ(Start, fFiducialLength);
        else if (geoDet == kPDHD)
            TrkAnodeCrossing = geoHighX.isInsideYZ(Start, fFiducialLength) || geoLowX.isInsideYZ(Start, fFiducialLength);

        
        /* COMPARE TO AGNOCHECKS
        std::vector<ana::LinearRegression> reg_side(2);
        std::vector<VecPtrHit> vph_side(2);
        for (PtrHit const& ph_hit : vph_mu) {
            if (ph_hit->View() != geo::kW) continue;
            int side = ana::tpc2side[geoDet][ph_hit->WireID().TPC];
            if (side == -1) continue;
            double z = GetSpace(ph_hit->WireID());
            double t = ph_hit->PeakTime() * fTick2cm;
            reg_side[side].add(z, t);
            vph_side[side].push_back(ph_hit);
        }
        AgnoTagTrkEndLowN = reg_side[0].n < ana::LinearRegression::nmin
            && reg_side[1].n < ana::LinearRegression::nmin;

        if (reg_side[0].n >= ana::LinearRegression::nmin
            && reg_side[1].n >= ana::LinearRegression::nmin
        ) {
            reg_side[0].compute();
            reg_side[1].compute();
            std::vector<std::pair<VecPtrHit::iterator, VecPtrHit::iterator>> ends_side(2);
            ends_side[0] = std::minmax_element(
                vph_side[0].begin(), vph_side[0].end(),
                [&](PtrHit const& a, PtrHit const& b) -> bool {
                    double sa = reg_side[0].projection(GetSpace(a->WireID()), a->PeakTime() * fTick2cm);
                    double sb = reg_side[0].projection(GetSpace(b->WireID()), b->PeakTime() * fTick2cm);
                    return sa < sb;
                } 
            );
            ends_side[1] = std::minmax_element(
                vph_side[1].begin(), vph_side[1].end(),
                [&](PtrHit const& a, PtrHit const& b) -> bool {
                    double sa = reg_side[1].projection(GetSpace(a->WireID()), a->PeakTime() * fTick2cm);
                    double sb = reg_side[1].projection(GetSpace(b->WireID()), b->PeakTime() * fTick2cm);
                    return sa < sb;
                } 
            );
            std::vector<std::pair<PtrHit, PtrHit>> pairs = {
                { *ends_side[0].first, *ends_side[1].first },
                { *ends_side[0].first, *ends_side[1].second },
                { *ends_side[0].second, *ends_side[1].first },
                { *ends_side[0].second, *ends_side[1].second }
            };
            AgnoTagTrkEndBadCC = std::find_if(
                pairs.begin(), pairs.end(),
                [&](std::pair<PtrHit, PtrHit> const& p) -> bool {
                    return GetDistance(p.first, p.second) < 2 * fCathodeGap;
                }
            ) == pairs.end();
        } else AgnoTagTrkEndBadCC = false;
        // End of Agnochecks
        */


        // SUPPOSITION: Muon is downward
        TrkRegDirZ = End.Z() > Start.Z() ? 1 : -1;
        ana::SortedHits sh_mu = GetSortedHits(vph_mu, TrkRegDirZ);
        TrkHitError = !sh_mu;

        LOG(!TrkHitError);
        if (!TrkHitError) {
            TrkStartHit = GetHit(sh_mu.start);
            TrkEndHit = GetHit(sh_mu.end);
            TrkReg = sh_mu.end_reg(geoDet);
            TrkHitEndInWindow = wireWindow.isInside(TrkEndHit.tick, fFiducialLength / fTick2cm);

            if (sh_mu.is_cc()) {
                if (abs(sh_mu.cc.first->PeakTime()-sh_mu.cc.second->PeakTime())*fTick2cm < 3 * fCathodeGap)
                    TrkHitCathodeCrossing = kAlignedHitOnBothSides;
                else
                    TrkHitCathodeCrossing = kHitOnBothSides;
            } else
                TrkHitCathodeCrossing = kNoCC;

            /* ASSUMS DOWNWARD MUONS FOR PDVD */
            if (geoDet == kPDVD)
                TrkHitAnodeCrossing = TrkStartHit.section < 4 
                    && geoHighX.z.isInside(TrkStartHit.space, fFiducialLength)
                    && wireWindow.isInside(TrkStartHit.tick, fFiducialLength/fTick2cm);
            else if (geoDet == kPDHD)
                TrkHitAnodeCrossing = 
                    geoHighX.z.isInside(TrkStartHit.space, fFiducialLength)
                    && wireWindow.isInside(TrkStartHit.tick, fFiducialLength/fTick2cm);
            
            LOG(TrkHitCathodeCrossing == kAlignedHitOnBothSides);
            LOG(TrkHitAnodeCrossing);
            if (TrkHitCathodeCrossing == kAlignedHitOnBothSides) {
                if (geoDet == kPDVD) {
                    TrkEndHitX = geoLowX.x.max - abs(TrkEndHit.tick - sh_mu.cc.second->PeakTime()) * fTick2cm;
                    TrkHitEndInVolumeX = geoLowX.x.isInside(TrkEndHitX, fFiducialLength);
                } else if (geoDet == kPDHD) {
                    int cc_sec = ana::tpc2sec[geoDet][sh_mu.cc.second->WireID().TPC];
                    if (cc_sec == 0) {
                        TrkEndHitX = geoLowX.x.max - abs(TrkEndHit.tick - sh_mu.cc.second->PeakTime()) * fTick2cm;
                        TrkHitEndInVolumeX = geoLowX.x.isInside(TrkEndHitX, fFiducialLength);
                    } else if (cc_sec == 1) {
                        TrkEndHitX = geoHighX.x.min + abs(TrkEndHit.tick - sh_mu.cc.second->PeakTime()) * fTick2cm;
                        TrkHitEndInVolumeX = geoHighX.x.isInside(TrkEndHitX, fFiducialLength);
                    }
                }
            } else if (TrkHitAnodeCrossing) {
                if (geoDet == kPDVD) {
                    TrkEndHitX = geoHighX.x.max - abs(TrkEndHit.tick - TrkStartHit.tick) * fTick2cm;
                    TrkHitEndInVolumeX = geoHighX.x.isInside(TrkEndHitX, fFiducialLength);
                } else if (geoDet == kPDHD) {
                    if (TrkStartHit.section == 0) {
                        TrkEndHitX = geoLowX.x.min + abs(TrkEndHit.tick - TrkStartHit.tick) * fTick2cm;
                        TrkHitEndInVolumeX = geoLowX.x.isInside(TrkEndHitX, fFiducialLength);
                    } else if (TrkStartHit.section == 1) {
                        TrkEndHitX = geoHighX.x.max - abs(TrkEndHit.tick - TrkStartHit.tick) * fTick2cm;
                        TrkHitEndInVolumeX = geoHighX.x.isInside(TrkEndHitX, fFiducialLength);
                    }
                }
            }

            VecPtrHit vph_mu_sec;
            for (PtrHit const& ph_mu : sh_mu.vph) {
                if (ph_mu->View() != geo::kW) continue;
                ana::Hit hit = GetHit(ph_mu);
                TrkHits.push_back(hit);

                if (hit.section != sh_mu.end_sec()) continue;
                vph_mu_sec.push_back(ph_mu);
            }
                    

            /* COMPARE TO AGNOCHECKS
            HitPtrPair ends = GetTrackEndsHits(vph_mu);
            if (!LOG(ends.first && ends.second)) continue;
            bool increasing_z = IsUpright(*pt_ev)
                ? pt_ev->End().Z() > pt_ev->Start().Z()
                : pt_ev->Start().Z() > pt_ev->End().Z();
            int dir_z = increasing_z ? 1 : -1;
            float fz = GetSpace(ends.first->WireID());
            float sz = GetSpace(ends.second->WireID());
            HitPtr end = (sz-fz) * dir_z > 0 ? ends.second : ends.first;
            TrkEndHit = GetHit(end);
            */

            VecPtrHit vph_end_sec;
            for (PtrHit const& ph_ev : vph_ev) {
                if (ph_ev->View() != geo::kW) continue;
                if (ana::tpc2sec[geoDet][ph_ev->WireID().TPC] != TrkEndHit.section) continue;
                vph_end_sec.push_back(ph_ev);
            }

            // integrate charges around muon endpoint
            PandoraSphereEnergy = 0;
            for (PtrHit const& ph_ev : vph_end_sec) {
                if (GetDistance(ph_ev, sh_mu.end) > fMichelRadius) continue;

                PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                if (pt_hit && pt_hit->Length() > fTrackLengthCut) continue;

                PtrShw ps_hit = fop_hit2shw.at(ph_ev.key());
                if (ps_hit) PandoraSphereHasShower = true;

                ana::Hit hit = GetHit(ph_ev);
                PandoraSphereEnergy += ph_ev->Integral();
                PandoraSphereHits.push_back(hit);

                float da = (hit.vec(fTick2cm) - TrkEndHit.vec(fTick2cm)).angle() - TrkReg.theta(TrkRegDirZ);
                da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                PandoraSphereHitMuonAngle.push_back(da);

                if (std::find_if(
                    vph_mi.begin(), vph_mi.end(),
                    [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
                ) != vph_mi.end()) 
                    PandoraSphereEnergyTP = ph_ev->Integral();
            }

            // Cone
            for (PtrHit const& ph_ev : vph_end_sec) {
                PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                if (pt_hit && pt_hit->Length() > fTrackLengthCut) continue;
                if (GetDistance(ph_ev, sh_mu.end) > 10) continue;

                PtrShw ps_hit = fop_hit2shw.at(ph_ev.key());
                if (ps_hit) PandoraBaryHasShower = true;

                PandoraBaryHits.push_back(GetHit(ph_ev));
            }

            LOG(PandoraBaryHits.size());
            if (PandoraBaryHits.size()) {
                PandoraBary = PandoraBaryHits.barycenter(fTick2cm);
                ana::Vec2 end_bary = PandoraBary - TrkEndHit.vec(fTick2cm);
                PandoraBaryAngle = end_bary.angle();
                float da = PandoraBaryAngle - TrkReg.theta(TrkRegDirZ);
                da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                PandoraBaryMuonAngle = da;

                // float angle = end_bary.angle();
                PandoraConeEnergy = 0;
                PandoraConeEnergyTP = 0;
                for (PtrHit const& ph_ev : vph_end_sec) {
                    ana::Hit hit = GetHit(ph_ev);
                    PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                    if (pt_hit && pt_hit->Length() > fTrackLengthCut) continue;

                    float dist = GetDistance(ph_ev, sh_mu.end);
                    if (dist > 30) continue;

                    ana::Vec2 end_hit = hit.vec(fTick2cm) - TrkEndHit.vec(fTick2cm);
                    float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

                    if (dist > 5
                        && cosa < cos(30.F * TMath::DegToRad())
                    ) continue;

                    bool tp = std::find_if(
                        vph_mi.begin(), vph_mi.end(),
                        [&ph_ev](PtrHit const& h) -> bool { return h.key() == ph_ev.key(); }
                    ) != vph_mi.end();

                    PandoraKeyholeEnergy += ph_ev->Integral();
                    PandoraKeyholeHits.push_back(hit);

                    if (tp) PandoraKeyholeEnergyTP += ph_ev->Integral();

                    if (cosa < cos(30.F * TMath::DegToRad())) continue;

                    PandoraConeEnergy += ph_ev->Integral();
                    PandoraConeHits.push_back(hit);

                    if (tp) PandoraConeEnergyTP += ph_ev->Integral();
                }
            }

            // End dQdx
            PandoraTrkHitdQdx = GetdQdx(vph_mu_sec, fRegN);

            ana::Bragg bragg = GetBragg(
                sh_mu.vph,
                sh_mu.end,
                // vph_mu,
                // end,
                pt_ev,
                vph_ev,
                fop_hit2trk,
                // TrkReg,
                { fBodyDistance, fRegN, fTrackLengthCut, fNearbyRadius }
            );
            BraggError = bragg.error;

            LOG(BraggError == kNoError);
            if (BraggError == kNoError) {
                MIPdQdx = bragg.mip_dQdx;
                BraggdQdx = bragg.max_dQdx;
                BraggEndHit = GetHit(bragg.end);

                for (PtrHit const& ph_mu : bragg.vph_mu)
                    BraggMuonHits.push_back(GetHit(ph_mu));

                VecPtrHit vph_near;
                for (PtrHit const& ph_ev : vph_end_sec) {
                    if (GetDistance(ph_ev, bragg.end) > fNearbyRadius) continue;

                    // not from mu
                    if (std::find_if(
                        bragg.vph_mu.begin(), bragg.vph_mu.end(),
                        [&ph_ev](PtrHit const& ph) { return ph.key() == ph_ev.key(); }
                    ) != bragg.vph_mu.end()) continue;
                    
                    // not from other long track
                    PtrTrk pt_hit = fop_hit2trk.at(ph_ev.key());
                    if (pt_hit
                        && pt_hit.key() != pt_ev.key()
                        && pt_hit->Length() > fTrackLengthCut
                    ) continue;

                    vph_near.push_back(ph_ev);
                }

                // Sphere
                BraggSphereEnergy = 0;
                BraggSphereEnergyTP = 0;
                for (PtrHit const& ph_near : vph_near) {
                    if (GetDistance(ph_near, bragg.end) > fMichelRadius) continue;
                    // PtrShw ps_hit = fop_hit2shw.at(iph->key());
                    BraggSphereEnergy += ph_near->Integral();
                    ana::Hit hit = GetHit(ph_near);
                    BraggSphereHits.push_back(hit);

                    float da = (hit.vec(fTick2cm) - BraggEndHit.vec(fTick2cm)).angle() - TrkReg.theta(TrkRegDirZ);
                    da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                    BraggSphereHitMuonAngle.push_back(da);

                    if (std::find_if(
                        vph_mi.begin(), vph_mi.end(),
                        [&ph_near](PtrHit const& h) -> bool { return h.key() == ph_near.key(); }
                    ) != vph_mi.end()) 
                        BraggSphereEnergyTP = ph_near->Integral();
                }

                // Cone
                for (PtrHit const& ph_near : vph_near) {
                    if (GetDistance(ph_near, bragg.end) > 10) continue;
                    NearbyBaryHits.push_back(GetHit(ph_near));
                }

                LOG(NearbyBaryHits.size());
                if (NearbyBaryHits.size()) {
                    NearbyBary = NearbyBaryHits.barycenter(fTick2cm);
                    ana::Vec2 end_bary = NearbyBary - BraggEndHit.vec(fTick2cm);
                    NearbyBaryAngle = end_bary.angle();
                    float da = NearbyBaryAngle - TrkReg.theta(TrkRegDirZ);
                    da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                    NearbyBaryMuonAngle = da;

                    // float angle = end_bary.angle();
                    BraggConeEnergy = 0;
                    BraggConeEnergyTP = 0;
                    for (PtrHit const& ph_near : vph_near) {
                        float dist = GetDistance(ph_near, bragg.end);
                        if (dist > 30) continue;

                        ana::Hit hit = GetHit(ph_near);
                        ana::Vec2 end_hit = hit.vec(fTick2cm) - BraggEndHit.vec(fTick2cm);
                        float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

                        if (dist > 5
                            && cosa < cos(30.F * TMath::DegToRad())
                        ) continue;

                        bool tp = std::find_if(
                            vph_mi.begin(), vph_mi.end(),
                            [&ph_near](PtrHit const& h) -> bool { return h.key() == ph_near.key(); }
                        ) != vph_mi.end();

                        BraggKeyholeEnergy += ph_near->Integral();
                        BraggKeyholeHits.push_back(hit);

                        if (tp) BraggKeyholeEnergyTP += ph_near->Integral();

                        if (cosa < cos(30.F * TMath::DegToRad())) continue;

                        BraggConeEnergy += ph_near->Integral();
                        BraggConeHits.push_back(hit);

                        if (tp) BraggConeEnergyTP += ph_near->Integral();
                    }
                }
            }
        }

        // Truth Information
        LOG(mcp);
        if (mcp) {
            TruePdg = mcp->PdgCode();
            TrueEndProcess = mcp->EndProcess();
            MuonTrueStartPoint = ana::Point(mcp->Position().Vect());
            MuonTrueEndPoint = ana::Point(mcp->EndPosition().Vect());
            MuonTrueEndEnergy = (mcp->EndE() - mcp->Mass()) * 1e3; // MeV

            if (geoDet == kPDVD)
                TrueDownward = mcp->Position(0).X() > mcp->EndPosition().X();
            else if (geoDet == kPDHD)
                TrueDownward = mcp->Position(0).Y() > mcp->EndPosition().Y();

            MuonTrueRegDirZ = mcp->EndZ() > mcp->Vz() ? 1 : -1;
            ana::SortedHits sh_mcp = GetSortedHits(vph_mcp_mu, MuonTrueRegDirZ);

            LOG(sh_mcp);
            if (sh_mcp) {
                MuonTrueStartHit = GetHit(sh_mcp.start);
                MuonTrueEndHit = GetHit(sh_mcp.end);
                MuonTrueReg = sh_mcp.end_reg(geoDet);

                LOG(mcp_mi);
                if (mcp_mi) {
                    TrueHasMichel = (
                        geoHighX.isInside(mcp_mi->Position().Vect(), 20.F)
                        || geoLowX.isInside(mcp_mi->Position().Vect(), 20.F)
                    ) ? kHasMichelFiducial : (
                        geoHighX.isInside(mcp_mi->Position().Vect())
                        || geoLowX.isInside(mcp_mi->EndPosition().Vect())
                        ? kHasMichelInside
                        : kHasMichelOutside
                    );
                    MichelTrueEnergy = (mcp_mi->E() - mcp_mi->Mass()) * 1e3;


                    PtrTrk pt_mi = ana::mcp2trk(mcp_mi, vpt_ev, clockData, fmp_trk2hit);
                    MichelTrackLength = pt_mi ? pt_mi->Length() : -1.F;
                    PtrShw ps_mi = ana::mcp2shw(mcp_mi, vps_ev, clockData, fmp_shw2hit);
                    MichelShowerLength = ps_mi ? ps_mi->Length() : -1.F;



                    float mu_end_angle = sh_mcp.end_reg(geoDet).theta(MuonTrueRegDirZ);
                    for (PtrHit const& ph_mi : vph_mi) {
                        if (ph_mi->View() != geo::kW) continue;
                        Hit hit = GetHit(ph_mi);
                        MichelHits.push_back(hit);
                        if (hit.section != MuonTrueEndHit.section) {
                            MichelHitMuonAngle.push_back(100);
                            continue;
                        }
                        float mu_hit_angle = (hit.vec(fTick2cm) - MuonTrueEndHit.vec(fTick2cm)).angle() - mu_end_angle;
                        mu_hit_angle = abs(mu_hit_angle) > TMath::Pi()
                            ? mu_hit_angle - (mu_hit_angle>0 ? 1 : -1) * 2 * TMath::Pi()
                            : mu_hit_angle;
                        MichelHitMuonAngle.push_back(mu_hit_angle);
                    }
                    MichelHitEnergy = MichelHits.energy();

                    // Cone
                    Hits bary_hits;
                    for (PtrHit const& ph_mi : vph_mi) {
                        if (ph_mi->View() != geo::kW) continue;
                        if (GetDistance(ph_mi, sh_mcp.end) > 10) continue;
                        bary_hits.push_back(GetHit(ph_mi));
                    }
                    MichelBaryNHit = bary_hits.size();

                    LOG(MichelBaryNHit);
                    if (bary_hits.size()) {
                        MichelBary = bary_hits.barycenter(fTick2cm);
                        ana::Vec2 end_bary = MichelBary - MuonTrueEndHit.vec(fTick2cm);
                        MichelBaryAngle = end_bary.angle();
                        float da = MichelBaryAngle - mu_end_angle;
                        da = abs(da) > M_PI ? da - (da>0 ? 1 : -1) * 2 * M_PI : da;
                        MichelBaryMuonAngle = da;

                        // float angle = end_bary.angle();
                        MichelConeEnergy = 0;
                        for (PtrHit const& ph_mi : vph_mi) {
                            if (ph_mi->View() != geo::kW) continue;
                            float dist = GetDistance(ph_mi, sh_mcp.end);
                            if (dist > 30) continue;

                            ana::Vec2 end_hit = GetHit(ph_mi).vec(fTick2cm) - MuonTrueEndHit.vec(fTick2cm);
                            float cosa = end_bary.dot(end_hit) / (end_bary.norm() * end_hit.norm());

                            if (dist > 5 
                                && cosa < cos(30.F * TMath::DegToRad())
                            ) continue;
                            
                            MichelKeyholeEnergy += ph_mi->Integral();

                            if (cosa < cos(30.F * TMath::DegToRad())) continue;

                            MichelConeEnergy += ph_mi->Integral();
                        }
                    }
                }
            }
        }


        tMuon->Fill();
        iMuon++;
        EventNMuon++;
    } // end of loop over tracks
    tEvent->Fill();
    iEvent++;
}

void ana::MichelAnalysis::beginJob() {}
void ana::MichelAnalysis::endJob() {}

void ana::MichelAnalysis::resetEvent() {
    EventNMuon = 0;
    EventiMuon.clear();
    EventHits.clear();
}
void ana::MichelAnalysis::resetMuon() {
    // Hit
    TrkHitCathodeCrossing = -1;
    TrkHitAnodeCrossing = false;
    TrkStartHit = ana::Hit{};
    TrkEndHit = ana::Hit{};
    TrkEndHitX = -500.F;
    TrkHitEndInVolumeX = false;
    TrkHitEndInWindow = false;
    TrkRegDirZ = 0;
    TrkReg = ana::LinearRegression{};
    TrkHits.clear();
    PandoraTrkHitdQdx.clear();
    PandoraSphereHits.clear();
    PandoraSphereHitMuonAngle.clear();
    PandoraSphereHasShower = false;
    PandoraBaryHasShower = false;
    PandoraSphereEnergy = -1.F;
    PandoraSphereEnergyTP = -1.F;

    // Cone
    PandoraBaryHits.clear();
    PandoraBary = ana::Vec2{0,0};
    PandoraBaryAngle = -10.F;
    PandoraBaryMuonAngle = -10.F;
    PandoraConeHits.clear();
    PandoraConeEnergy = -1.F;
    PandoraConeEnergyTP = -1.F;
    PandoraKeyholeHits.clear();
    PandoraKeyholeEnergy = -1.F;
    PandoraKeyholeEnergyTP = -1.F;


    // Bragg
    ana::Bragg bragg;
    BraggError = bragg.error;
    MIPdQdx = bragg.mip_dQdx;
    BraggdQdx = bragg.max_dQdx;
    BraggEndHit = ana::Hit{};
    BraggMuonHits.clear();
    BraggSphereHits.clear();
    BraggSphereHitMuonAngle.clear();
    BraggSphereEnergy = -1.F;
    BraggSphereEnergyTP = -1.F;

    NearbyBaryHits.clear();
    NearbyBary = ana::Vec2{0,0};
    NearbyBaryAngle = -10.F;
    NearbyBaryMuonAngle = -10.F;
    BraggConeHits.clear();
    BraggConeEnergy = -1.F;
    BraggConeEnergyTP = -1.F;
    BraggKeyholeHits.clear();
    BraggKeyholeEnergy = -1.F;
    BraggKeyholeEnergyTP = -1.F;

    // Muon Truth
    TruePdg = 0;
    TrueEndProcess = "";
    MuonTrueStartPoint = ana::Point{};
    MuonTrueEndPoint = ana::Point{};
    MuonTrueEndEnergy = -1.F;
    TrueDownward = false;
    MuonTrueStartHit = ana::Hit{};
    MuonTrueEndHit = ana::Hit{};
    MuonTrueRegDirZ = 0;
    MuonTrueReg = ana::LinearRegression{};

    // Michel Truth
    TrueHasMichel = kNoMichel;
    MichelTrueEnergy = -1.F;
    MichelTrackLength = -1.F;
    MichelShowerLength = -1.F;
    MichelHits.clear();
    MichelHitMuonAngle.clear();
    MichelHitEnergy = -1.F;

    MichelBaryNHit = 0;
    MichelBary = ana::Vec2{0,0};
    MichelBaryAngle = -10.F;
    MichelBaryMuonAngle = -10.F;
    MichelConeEnergy = -1.F;
    MichelKeyholeEnergy = -1.F;
}

// HitPtrPair ana::MichelAnalysis::GetTrackEndsHits(
//     HitPtrVec const& vp_hit,
//     HitPtrPair *pp_cathode_crossing,
//     HitPtrVec *vp_tpc_crossing,
//     HitPtrVec *vp_sorted_hit,
//     geo::View_t view
// ) {
//     // minimum number of hits to perform a linear regression
//     unsigned const nmin = ana::LinearRegression::nmin;

//     // split volume at de cathode
//     auto cathodeSide =
//         geoDet == kPDVD
//         ? [](geo::TPCID::TPCID_t tpc) -> int {
//                 return tpc >= 8 ? 1 : 0;
//             }
//             // geoDet == kPDHD
//         : [](geo::TPCID::TPCID_t tpc) -> int {
//                 return (tpc == 1 || tpc == 5)
//                     ? 0
//                     : ((tpc == 2 || tpc == 6) ? 1 : -1);
//             };


//     // linear regression on each side to have a curvilinear coordinate of each hit inside a track
//     // z = m*t + p
//     // struct LinearRegression {
//     //     unsigned n=0;
//     //     double mz=0, mt=0, mz2=0, mt2=0, mzt=0;
//     //     void add(double z, double t) {
//     //         mz+=z; mt+=t; mz2+=z*z; mt2+=t*t; mzt+=z*t; n++;
//     //     }
//     //     void normalize() {
//     //         mz/=n; mt/=n; mz2/=n; mt2/=n; mzt/=n;
//     //     }
//     //     double cov() const { return mzt - mz*mt; }
//     //     double varz() const { return mz2 - mz*mz; }
//     //     double vart() const { return mt2 - mt*mt; }
//     //     double m() const { return n<nmin ? 0 : cov()/vart(); }
//     //     double p() const { return mz - m()*mt; }
//     //     double r2() const { return n<nmin ? 0 : cov()*cov() / (varz()*vart()); }
//     //     double projection(double z, double t) const {
//     //         return (t + m()*(z-p())) / (1 + m()*m());
//     //     }
//     // };

//     std::vector<ana::LinearRegression> side_reg(2);
//     for (HitPtr const& p_hit : vp_hit) {
//         if (p_hit->View() != view) continue;
//         int side = cathodeSide(p_hit->WireID().TPC);
//         if (side == -1) continue; // skip hits on the other side of the anodes
//         double z = GetSpace(p_hit->WireID());
//         double t = p_hit->PeakTime() * fTick2cm;
//         side_reg[side].add(z, t);
//     }

//     // if not enough hits on both sides, return empty pair
//     if (side_reg[0].n < nmin && side_reg[1].n < nmin) return {};

//     // compute average from sum
//     for (ana::LinearRegression& reg : side_reg) reg.compute();

//     // find the track ends on each side of the cathode
//     std::vector<HitPtrPair> side_ends(2);
//     std::vector<Bounds<double>> side_mimmax(2);
//     for (HitPtr const& p_hit : vp_hit) {
//         if (p_hit->View() != view) continue;
//         int side = cathodeSide(p_hit->WireID().TPC);
//         if (side == -1) continue; // skip hits on the other side of the anodes
//         double s = side_reg[side].projection(
//             GetSpace(p_hit->WireID()),
//             p_hit->PeakTime() * fTick2cm
//         );
//         if (s > side_mimmax[side].max) {
//             side_mimmax[side].max = s;
//             side_ends[side].second = p_hit;
//         }
//         if (s < side_mimmax[side].min) {
//             side_mimmax[side].min = s;
//             side_ends[side].first = p_hit;
//         }
//     }

//     // if hits are all on one side, and no other info is requested
//     if (!vp_tpc_crossing && !vp_sorted_hit) {
//         if (side_reg[0].n < nmin)
//             return side_ends[1];
//         else if (side_reg[1].n < nmin)
//             return side_ends[0];
//     }
    
//     // given the ends of two pieces of track, find the closest ends
//     auto closestHits = [&](
//         HitPtrPair const& pph1,
//         HitPtrPair const& pph2,
//         double dmin,
//         HitPtrPair *otherHits = nullptr
//     ) -> HitPtrPair {

//         // all combinations of pairs
//         std::vector<HitPtrPair> pairs = {
//             { pph1.first, pph2.first },
//             { pph1.first, pph2.second },
//             { pph1.second, pph2.second },
//             { pph1.second, pph2.first }
//         };

//         // distance squared between all pairs
//         std::vector<double> d2s(4, 0);
//         for (unsigned i=0; i<4; i++) {
//             double zf = GetSpace(pairs[i].first->WireID());
//             double tf = pairs[i].first->PeakTime() * fTick2cm;
//             double zs = GetSpace(pairs[i].second->WireID());
//             double ts = pairs[i].second->PeakTime() * fTick2cm;
//             d2s[i] = pow(zf-zs,2) + pow(tf-ts,2);
//         }

//         // find all distances under dmin threshold
//         std::vector<unsigned> candidates_idx;
//         std::vector<double>::iterator it = d2s.begin();
//         while ((it = std::find_if(
//                 it,
//                 d2s.end(),
//                 [dmin](double d2) { return d2 < dmin*dmin; }
//             )) != d2s.end())
//             candidates_idx.push_back(std::distance(d2s.begin(), it++));
        
//         // no candidates found
//         if (candidates_idx.empty()) return {};

//         // get the closest pair
//         unsigned closest_idx = *std::min_element(
//             candidates_idx.begin(),
//             candidates_idx.end(),
//             [&d2s](unsigned i, unsigned j) { return d2s[i] < d2s[j]; });
        
//         // if outermost hits are requested, get the outermost pair
//         if (otherHits) {
//             unsigned other_idx = (closest_idx+2) % 4; // opposite pair
//             otherHits->first = pairs[other_idx].first;
//             otherHits->second = pairs[other_idx].second;
//         }
//         return pairs[closest_idx];
//     };

//     HitPtrPair trk_ends, cathode_crossing;
//     if (side_reg[0].n < nmin)
//         trk_ends = side_ends[1];
//     else if (side_reg[1].n < nmin)
//         trk_ends = side_ends[0];
//     else
//         cathode_crossing = closestHits(
//             side_ends[0],
//             side_ends[1],
//             2*fCathodeGap,
//             &trk_ends
//         );

//     // if cathode crossing info is requested
//     if (pp_cathode_crossing)
//         *pp_cathode_crossing = cathode_crossing;
    
//     // if no tpc crossing info is needed
//     if (geoDet == kPDHD || !vp_tpc_crossing) {
//         return trk_ends;
//     }

//     std::vector<HitPtrPair> per_sec_ends(ana::n_sec[geoDet]);
    
//     // if a sorted list of hits is requested
//     // if (vvp_sec_sorted_hits) {
//     //     // get a sorted list of hits for each section (ie. pair of TPCs)
//     //     vvp_sec_sorted_hits->clear();
//     //     vvp_sec_sorted_hits->resize(ana::n_sec[geoDet]);
//     //     for (HitPtr const& p_hit : vp_hit) {
//     //         if (p_hit->View() != view) continue;
//     //         int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
//     //         if (s == -1) continue;
//     //         vvp_sec_sorted_hits->at(s).push_back(p_hit);
//     //     }

//     //     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//     //         int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
//     //         std::sort(
//     //             vvp_sec_sorted_hits->at(s).begin(), 
//     //             vvp_sec_sorted_hits->at(s).end(),
//     //             [&, &reg=side_reg[side]](
//     //                 HitPtr const& h1, HitPtr const& h2
//     //             ) -> bool {
//     //                 double const s1 = reg.projection(
//     //                     GetSpace(h1->WireID()),
//     //                     h1->PeakTime() * fTick2cm
//     //                 );
//     //                 double const s2 = reg.projection(
//     //                     GetSpace(h2->WireID()),
//     //                     h2->PeakTime() * fTick2cm
//     //                 );
//     //                 return s1 < s2;
//     //             }
//     //         );
//     //     }

//     //     // get the track ends for each section
//     //     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//     //         if (vvp_sec_sorted_hits->at(s).size() < nmin) continue;
//     //         per_sec_ends[s].first = vvp_sec_sorted_hits->at(s).front();
//     //         per_sec_ends[s].second = vvp_sec_sorted_hits->at(s).back();
//     //     }

//     std::vector<HitPtrVec> vp_sec_hit(ana::n_sec[geoDet]);
//     for (HitPtr const& p_hit : vp_hit) {
//         if (p_hit->View() != view) continue;
//         int s = ana::tpc2sec[geoDet][p_hit->WireID().TPC];
//         if (s == -1) continue;
//         vp_sec_hit[s].push_back(p_hit);
//     }
//     // THIS CAUSES A SEGFAULT FOR SOME REASON???
//     // auto side_sort = [&](int side) {
//     //     return [&](HitPtr const& h1, HitPtr const& h2) -> bool {
//     //         double const s1 = side_reg[side].projection(
//     //             GetSpace(h1->WireID()),
//     //             h1->PeakTime() * fTick2cm
//     //         );
//     //         double const s2 = side_reg[side].projection(
//     //             GetSpace(h2->WireID()),
//     //             h2->PeakTime() * fTick2cm
//     //         );
//     //         return s1 < s2;
//     //     };
//     // };

//     if (vp_sorted_hit) {
//         vp_sorted_hit->clear();

//         for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//             HitPtrVec& vp_sec_sorted = vp_sec_hit[s];
//             if (vp_sec_sorted.size() < nmin) continue;

//             int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
//             std::sort(
//                 vp_sec_sorted.begin(), 
//                 vp_sec_sorted.end(),
//                 [&, &reg=side_reg[side]](
//                     HitPtr const& h1, HitPtr const& h2
//                 ) -> bool {
//                     double const s1 = reg.projection(
//                         GetSpace(h1->WireID()),
//                         h1->PeakTime() * fTick2cm
//                     );
//                     double const s2 = reg.projection(
//                         GetSpace(h2->WireID()),
//                         h2->PeakTime() * fTick2cm
//                     );
//                     return s1 < s2;
//                 }
//             );

//             // get the track ends for each section
//             per_sec_ends[s].first = vp_sec_sorted.front();
//             per_sec_ends[s].second = vp_sec_sorted.back();

//             vp_sorted_hit->insert(
//                 vp_sorted_hit->end(),
//                 vp_sec_sorted.begin(), vp_sec_sorted.end()
//             );
//         }
//     } else { // only get the minmax ends of each section
//         for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//             if (vp_sec_hit[s].size() < nmin) continue;
//             int side = s >= ana::n_sec[geoDet]/2 ? 1 : 0;
//             auto minmax = std::minmax_element(
//                 vp_sec_hit[s].begin(),
//                 vp_sec_hit[s].end(),
//                 [&, &reg=side_reg[side]](HitPtr const& h1, HitPtr const& h2) -> bool {
//                     double const s1 = reg.projection(
//                         GetSpace(h1->WireID()),
//                         h1->PeakTime() * fTick2cm
//                     );
//                     double const s2 = reg.projection(
//                         GetSpace(h2->WireID()),
//                         h2->PeakTime() * fTick2cm
//                     );
//                     return s1 < s2;
//                 }
//             );
//             per_sec_ends[s].first = *minmax.first;
//             per_sec_ends[s].second = *minmax.second;
//         }
//     }


//     // get the hits that are at the boundaries of two sections
//     HitPtrVec tpc_crossing;
//     bool prev = false;
//     for (unsigned s=0; s<ana::n_sec[geoDet]; s++) {
//         if (per_sec_ends[s].first.isNull()) {
//             prev = false;
//             continue;
//         }
//         if (prev) {
//             HitPtrPair const pp_tpc_crossing = closestHits(
//                 per_sec_ends[s-1], per_sec_ends[s], 2
//             );
//             if (pp_tpc_crossing.first.isNonnull()) {
//                 tpc_crossing.push_back(pp_tpc_crossing.first);
//                 tpc_crossing.push_back(pp_tpc_crossing.second);
//             }
//         }
//         if ((geoDet == kPDVD && s == 3) || (geoDet == kPDHD && s == 1))
//             prev = false; // cathode crossing
//         else
//             prev = true;
//     }

//     *vp_tpc_crossing = tpc_crossing;
//     return trk_ends;
// }

DEFINE_ART_MODULE(ana::MichelAnalysis)