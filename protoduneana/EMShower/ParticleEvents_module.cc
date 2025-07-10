////////////////////////////////////////////////////////////////////////
// Class:       ParticleEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyseEvents_module.cc
//
// Generated at Tue Mar 25 09:54:37 2025 by Victor Chalamet using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

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
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
// #include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcore/Geometry/WireReadout.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/ArtDataHelper/MVAReader.h"
// #include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// ROOT includes
#include <TTree.h>

// C++ includes
#include <vector>
#include <string>
#include <iterator>
#include <cmath>

namespace analyse {
    class ParticleEvents;
}

class analyse::ParticleEvents : public art::EDAnalyzer {
public:
    explicit ParticleEvents(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    ParticleEvents(ParticleEvents const&) = delete;
    ParticleEvents(ParticleEvents&&) = delete;
    ParticleEvents& operator=(ParticleEvents const&) = delete;
    ParticleEvents& operator=(ParticleEvents&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

private:
    TTree *fTree;

    static const int NMaxParticles=50000;
    static const int NMaxMCParticles=3000000;
    static const int NMaxBGParticles=7000;
    unsigned int fEventID;
    // TODO: If too big switch to float data type
    // === MC information (SimEnergyDeposit) ===
    unsigned int fNMCParticles;
    double fMCParticleEnergy[NMaxMCParticles];
    double fMCTotalEnergy;
    double fMCParticleStartPositionX[NMaxMCParticles];
    double fMCParticleStartPositionY[NMaxMCParticles];
    double fMCParticleStartPositionZ[NMaxMCParticles];
    double fMCParticleEndPositionX[NMaxMCParticles];
    double fMCParticleEndPositionY[NMaxMCParticles];
    double fMCParticleEndPositionZ[NMaxMCParticles];
    int fMCParticlePdgCode[NMaxMCParticles];
    
    // === Reco information ===
    unsigned int fNParticles;
    double fEnergy;

    // === Track reco information ===
    double fTrackLength[NMaxParticles];
    double fTotalTrackLength;
    double fTotalTrackEnergy;
    double fTrackStartX[NMaxParticles];
    double fTrackStartY[NMaxParticles];
    double fTrackStartZ[NMaxParticles];
    double fTrackEndX[NMaxParticles];
    double fTrackEndY[NMaxParticles];
    double fTrackEndZ[NMaxParticles];
    double fTrackZenithAngle[NMaxParticles];
    double fTrackAzimuthAngle[NMaxParticles];

    // === Shower reco information ===
    double fShowerEnergy[NMaxParticles];
    double fTotalShowerEnergy;
    double fShowerStartX[NMaxParticles];
    double fShowerStartY[NMaxParticles];
    double fShowerStartZ[NMaxParticles];
    double fShowerEndX[NMaxParticles];
    double fShowerEndY[NMaxParticles];
    double fShowerEndZ[NMaxParticles];
    double fShowerLength[NMaxParticles];
    double fShowerDirectionX[NMaxParticles];
    double fShowerDirectionY[NMaxParticles];
    double fShowerDirectionZ[NMaxParticles];
    double fShowerOpenAngle[NMaxParticles];

    // === Background info ===
    // Ar39
    unsigned int fNAr39Particles;
    double fAr39Energy[NMaxBGParticles];
    double fAr39StartX[NMaxBGParticles];
    double fAr39StartY[NMaxBGParticles];
    double fAr39StartZ[NMaxBGParticles];
    int fAr39PdgCode[NMaxBGParticles];
    // Ar42
    unsigned int fNAr42Particles;
    double fAr42Energy[NMaxBGParticles];
    double fAr42StartX[NMaxBGParticles];
    double fAr42StartY[NMaxBGParticles];
    double fAr42StartZ[NMaxBGParticles];
    int fAr42PdgCode[NMaxBGParticles];
    // Kr85
    unsigned int fNKr85Particles;
    double fKr85Energy[NMaxBGParticles];
    double fKr85StartX[NMaxBGParticles];
    double fKr85StartY[NMaxBGParticles];
    double fKr85StartZ[NMaxBGParticles];
    int fKr85PdgCode[NMaxBGParticles];
    // Cosmic
    unsigned int fNCosmicParticles;
    double fCosmicEnergy[NMaxBGParticles];
    double fCosmicStartX[NMaxBGParticles];
    double fCosmicStartY[NMaxBGParticles];
    double fCosmicStartZ[NMaxBGParticles];
    int fCosmicPdgCode[NMaxBGParticles];

    // === Fcl parameters ===
    std::string fTruthLabel;
    std::string fHitLabel;
    std::string fTrackLabel;
    std::string fShowerLabel;
    std::string fPFParticleLabel;
    std::string fEdepLabel;
    std::string fAr39Label;
    std::string fAr42Label;
    std::string fKr85Label;
    std::string fCosmicLabel;

    geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

    void reset();
};

analyse::ParticleEvents::ParticleEvents(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} {
    // Call appropriate consumes<>() for any products to be retrieved by this module.
    fTruthLabel = p.get<std::string>("TruthLabel");
    fHitLabel = p.get<std::string>("HitLabel");
    fPFParticleLabel = p.get<std::string>("PFParticleLabel");
    fTrackLabel = p.get<std::string>("TrackLabel");
    fShowerLabel = p.get<std::string>("ShowerLabel");
    fEdepLabel = p.get<std::string>("EdepLabel");
    fAr39Label = p.get<std::string>("Ar39Label");
    fAr42Label = p.get<std::string>("Ar42Label");
    fKr85Label = p.get<std::string>("Kr85Label");
    fCosmicLabel = p.get<std::string>("CosmicLabel");
}

void analyse::ParticleEvents::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;

    // === Create the output TTree ===
    fTree = tfs->make<TTree>("tree", "Output TTree");

    // === Set the branches ===
    fTree->Branch("eventID", &fEventID, "EventID/i");
    fTree->Branch("nParticles", &fNParticles,"nParticles/i");
    fTree->Branch("nMCParticles", &fNMCParticles, "nMCParticles/i");
    fTree->Branch("nAr39Particles", &fNAr39Particles, "nAr39Particles/i");
    fTree->Branch("nAr42Particles", &fNAr42Particles, "nAr42Particles/i");
    fTree->Branch("nKr85Particles", &fNKr85Particles, "nKr85Particles/i");
    fTree->Branch("nCosmicParticles", &fNCosmicParticles, "nCosmicParticles/i");
    fTree->Branch("Energy", &fEnergy, "Energy/D");

    // === MC information (SimEnergyDeposit) ===
    fTree->Branch("MCParticleEnergy", fMCParticleEnergy, "MCParticleEnergy[nMCParticles]/D");
    fTree->Branch("MCTotalEnergy", &fMCTotalEnergy, "MCTotalEnergy/D");
    fTree->Branch("MCParticleStartPositionX", fMCParticleStartPositionX, "MCParticleStartPositionX[nMCParticles]/D");
    fTree->Branch("MCParticleStartPositionY", fMCParticleStartPositionY, "MCParticleStartPositionY[nMCParticles]/D");
    fTree->Branch("MCParticleStartPositionZ", fMCParticleStartPositionZ, "MCParticleStartPositionZ[nMCParticles]/D");
    fTree->Branch("MCParticleEndPositionX", fMCParticleEndPositionX, "MCParticleEndPositionX[nMCParticles]/D");
    fTree->Branch("MCParticleEndPositionY", fMCParticleEndPositionY, "MCParticleEndPositionY[nMCParticles]/D");
    fTree->Branch("MCParticleEndPositionZ", fMCParticleEndPositionZ, "MCParticleEndPositionZ[nMCParticles]/D");
    fTree->Branch("MCParticlePdgCode", fMCParticlePdgCode, "MCParticlePdgCode[nMCParticles]/i");
    
    // === Track reco information ===
    fTree->Branch("TrackLength", fTrackLength, "TrackLength[nParticles]/D");
    fTree->Branch("TotalTrackLength", &fTotalTrackLength, "TotalTrackLength/D");
    fTree->Branch("TotalTrackEnergy", &fTotalTrackEnergy, "TotalTrackEnergy/D");
    fTree->Branch("TrackStartX", fTrackStartX, "TrackStartX[nParticles]/D");
    fTree->Branch("TrackStartY", fTrackStartY, "TrackStartY[nParticles]/D");
    fTree->Branch("TrackStartZ", fTrackStartZ, "TrackStartZ[nParticles]/D");
    fTree->Branch("TrackEndX", fTrackEndX, "TrackEndX[nParticles]/D");
    fTree->Branch("TrackEndY", fTrackEndY, "TrackEndY[nParticles]/D");
    fTree->Branch("TrackEndZ", fTrackEndZ, "TrackEndZ[nParticles]/D");
    fTree->Branch("TrackZenithAngle", fTrackZenithAngle, "TrackZenithAngle[nParticles]/D");
    fTree->Branch("TrackAzimuthAngle", fTrackAzimuthAngle, "TrackAzimuthAngle[nParticles]/D");
    
    // === Shower reco information ===
    fTree->Branch("ShowerEnergy", fShowerEnergy, "ShowerEnergy[nParticles]/D");
    fTree->Branch("TotalShowerEnergy", &fTotalShowerEnergy, "TotalShowerEnergy/D");
    fTree->Branch("ShowerStartX", fShowerStartX, "ShowerStartX[nParticles]/D");
    fTree->Branch("ShowerStartY", fShowerStartY, "ShowerStartY[nParticles]/D");
    fTree->Branch("ShowerStartZ", fShowerStartZ, "ShowerStartZ[nParticles]/D");
    fTree->Branch("ShowerEndX", fShowerEndX, "ShowerEndX[nParticles]/D");
    fTree->Branch("ShowerEndY", fShowerEndY, "ShowerEndY[nParticles]/D");
    fTree->Branch("ShowerEndZ", fShowerEndZ, "ShowerEndZ[nParticles]/D");
    fTree->Branch("ShowerLength", fShowerLength, "ShowerLength[nParticles]/D");
    fTree->Branch("ShowerDirectionX", fShowerDirectionX, "ShowerDirectionX[nParticles]/D");
    fTree->Branch("ShowerDirectionY", fShowerDirectionY, "ShowerDirectionY[nParticles]/D");
    fTree->Branch("ShowerDirectionZ", fShowerDirectionZ, "ShowerDirectionZ[nParticles]/D");
    fTree->Branch("ShowerOpenAngle", fShowerOpenAngle, "ShowerOpenAngle[nParticles]/D");

    // === BG info ===
    // Ar39
    fTree->Branch("Ar39Energy", fAr39Energy, "Ar39Energy[nAr39Particles]/D");
    fTree->Branch("Ar39StartX", fAr39StartX, "Ar39StartX[nAr39Particles]/D");
    fTree->Branch("Ar39StartY", fAr39StartY, "Ar39StartY[nAr39Particles]/D");
    fTree->Branch("Ar39StartZ", fAr39StartZ, "Ar39StartZ[nAr39Particles]/D");
    fTree->Branch("Ar39PdgCode", fAr39PdgCode, "Ar39PdgCode[nAr39Particles]/i");
    // Ar42
    fTree->Branch("Ar42Energy", fAr42Energy, "Ar42Energy[nAr42Particles]/D");
    fTree->Branch("Ar42StartX", fAr42StartX, "Ar42StartX[nAr42Particles]/D");
    fTree->Branch("Ar42StartY", fAr42StartY, "Ar42StartY[nAr42Particles]/D");
    fTree->Branch("Ar42StartZ", fAr42StartZ, "Ar42StartZ[nAr42Particles]/D");
    fTree->Branch("Ar42PdgCode", fAr42PdgCode, "Ar42PdgCode[nAr42Particles]/i");
    // Kr85
    fTree->Branch("Kr85Energy", fKr85Energy, "Kr85Energy[nKr85Particles]/D");
    fTree->Branch("Kr85StartX", fKr85StartX, "Kr85StartX[nKr85Particles]/D");
    fTree->Branch("Kr85StartY", fKr85StartY, "Kr85StartY[nKr85Particles]/D");
    fTree->Branch("Kr85StartZ", fKr85StartZ, "Kr85StartZ[nKr85Particles]/D");
    fTree->Branch("Kr85PdgCode", fKr85PdgCode, "Kr85PdgCode[nKr85Particles]/i");
    // Cosmic
    fTree->Branch("CosmicEnergy", fCosmicEnergy, "CosmicEnergy[nCosmicParticles]/D");
    fTree->Branch("CosmicStartX", fCosmicStartX, "CosmicStartX[nCosmicParticles]/D");
    fTree->Branch("CosmicStartY", fCosmicStartY, "CosmicStartY[nCosmicParticles]/D");
    fTree->Branch("CosmicStartZ", fCosmicStartZ, "CosmicStartZ[nCosmicParticles]/D");
    fTree->Branch("CosmicPdgCode", fCosmicPdgCode, "CosmicPdgCode[nCosmicParticles]/i");
}

void analyse::ParticleEvents::analyze(art::Event const& e) {
    fEventID=e.id().event();
    reset();
    
    // === Track and Shower reco information ===
    // https://github.com/DUNE/duneana/blob/6ca410b3e06a19586e209e521299b232148fa9f1/duneana/PandoraAnalysis/PandoraAnalysis_module.cc
    const std::vector<art::Ptr<recob::PFParticle>> PFParticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(e, fPFParticleLabel);

    fNParticles = PFParticleVect.size();
    if(!fNParticles) {
        std::cout << "No PFParticles found!" << std::endl;
        return;
    }
  
    int iReco(0);
    for (const art::Ptr<recob::PFParticle> &pfp: PFParticleVect) {
        // === Get the track information ===
        if(dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, e, fPFParticleLabel, fTrackLabel)) {
            art::Ptr<recob::Track> track=dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, e, fPFParticleLabel, fTrackLabel);
            fTrackStartX[iReco] = track->Start().X();
            fTrackStartY[iReco] = track->Start().Y();
            fTrackStartZ[iReco] = track->Start().Z();
            fTrackEndX[iReco] = track->End().X();
            fTrackEndY[iReco] = track->End().Y();
            fTrackEndZ[iReco] = track->End().Z();
            fTrackLength[iReco] = track->Length();
            
            // Total track length
            if (fTrackLength[iReco]!=0) {
            fTrackZenithAngle[iReco] = track->ZenithAngle();
            fTrackAzimuthAngle[iReco] = track->AzimuthAngle();
            fTotalTrackLength += fTrackLength[iReco];
            }
        }
        // === Get the shower information ===
        if(dune_ana::DUNEAnaPFParticleUtils::IsShower(pfp, e, fPFParticleLabel, fShowerLabel)) {
            art::Ptr<recob::Shower> shower=dune_ana::DUNEAnaPFParticleUtils::GetShower(pfp, e, fPFParticleLabel, fShowerLabel);
            fShowerLength[iReco] = shower->Length();
            fShowerStartX[iReco] = shower->ShowerStart().X();
            fShowerStartY[iReco] = shower->ShowerStart().Y();
            fShowerStartZ[iReco] = shower->ShowerStart().Z();
            fShowerDirectionX[iReco] = shower->Direction().X();
            fShowerDirectionY[iReco] = shower->Direction().Y();
            fShowerDirectionZ[iReco] = shower->Direction().Z();
            fShowerEndX[iReco] = fShowerStartX[iReco]+fShowerLength[iReco]*fShowerDirectionX[iReco];
            fShowerEndY[iReco] = fShowerStartY[iReco]+fShowerLength[iReco]*fShowerDirectionY[iReco];
            fShowerEndZ[iReco] = fShowerStartZ[iReco]+fShowerLength[iReco]*fShowerDirectionZ[iReco];
            fShowerOpenAngle[iReco] = shower->OpenAngle();

            // Total shower energy
            const std::vector<double> &energy{shower->Energy()};
            fShowerEnergy[iReco] = energy[2]/1000.;
            if (fShowerEnergy[iReco]>0) fTotalShowerEnergy += fShowerEnergy[iReco];
        }
        iReco++;
    }

    // === Total energy computation ===
    fTotalTrackEnergy = fTotalTrackLength*2./1000.;
    fEnergy = fTotalShowerEnergy+fTotalTrackEnergy;
    
    // === MC information (SimEnergyDeposit) ===
    art::Handle<std::vector<sim::SimEnergyDeposit>> eDepHandle;
    if (!e.getByLabel(fEdepLabel, eDepHandle)) {
        std::cout << "Cannot load any energy deposits. Failing" << std::endl;
        mf::LogWarning("AnalyseEvents_module") << "Cannot load any energy deposits. Failing";
        return;
    }
    unsigned int iMC(0);
    for (auto const& eDep : *eDepHandle) {
        fMCParticleStartPositionX[iMC] = eDep.StartX();
        fMCParticleStartPositionY[iMC] = eDep.StartY();
        fMCParticleStartPositionZ[iMC] = eDep.StartZ();
        fMCParticleEndPositionX[iMC] = eDep.EndX();
        fMCParticleEndPositionY[iMC] = eDep.EndY();
        fMCParticleEndPositionZ[iMC] = eDep.EndZ();
        fMCParticleEnergy[iMC] = eDep.Energy()/1000.;
        fMCTotalEnergy += eDep.Energy()/1000.;
        fMCParticlePdgCode[iMC] = eDep.PdgCode();
        iMC++;
    }
    fNMCParticles = iMC;

    // === Background informartion ===
    // Ar39
    art::Handle<std::vector<simb::MCTruth>> Ar39Handle;
    if (!e.getByLabel(fAr39Label, Ar39Handle)) {
        std::cout << "Cannot load any ar39 background. Failing" << std::endl;
        mf::LogWarning("AnalyseEvents_module") << "Cannot load any ar39 background. Failing";
        return;
    }
    for (auto const& ar39 : *Ar39Handle) {
        unsigned const int nBGParticles = ar39.NParticles();
        for (unsigned int i=0; i<nBGParticles; i++) {
            simb::MCParticle const& particle = ar39.GetParticle(i);
            fAr39Energy[i] = particle.E();
            fAr39StartX[i] = particle.Position().X();
            fAr39StartY[i] = particle.Position().Y();
            fAr39StartZ[i] = particle.Position().Z();
            fAr39PdgCode[i] = particle.PdgCode();
            fNAr39Particles++;
        }
    }

    // Ar42
    art::Handle<std::vector<simb::MCTruth>> Ar42Handle;
    if (!e.getByLabel(fAr42Label, Ar42Handle)) {
        std::cout << "Cannot load any ar42 background. Failing" << std::endl;
        mf::LogWarning("AnalyseEvents_module") << "Cannot load any ar42 background. Failing";
        return;
    }
    for (auto const& ar42 : *Ar42Handle) {
        unsigned const int nBGParticles = ar42.NParticles();
        for (unsigned int i=0; i<nBGParticles; i++) {
            simb::MCParticle const& particle = ar42.GetParticle(i);
            fAr42Energy[i] = particle.E();
            fAr42StartX[i] = particle.Position().X();
            fAr42StartY[i] = particle.Position().Y();
            fAr42StartZ[i] = particle.Position().Z();
            fAr42PdgCode[i] = particle.PdgCode();
            fNAr42Particles++;
        }
    }

    // Kr85
    art::Handle<std::vector<simb::MCTruth>> Kr85Handle;
    if (!e.getByLabel(fKr85Label, Kr85Handle)) {
        std::cout << "Cannot load any kr85 background. Failing" << std::endl;
        mf::LogWarning("AnalyseEvents_module") << "Cannot load any kr85 background. Failing";
        return;
    }
    for (auto const& kr85 : *Kr85Handle) {
        unsigned const int nBGParticles = kr85.NParticles();
        for (unsigned int i=0; i<nBGParticles; i++) {
            simb::MCParticle const& particle = kr85.GetParticle(i);
            fKr85Energy[i] = particle.E();
            fKr85StartX[i] = particle.Position().X();
            fKr85StartY[i] = particle.Position().Y();
            fKr85StartZ[i] = particle.Position().Z();
            fKr85PdgCode[i] = particle.PdgCode();
            fNKr85Particles++;
        }
    }

    // Cosmics
    art::Handle<std::vector<simb::MCTruth>> CosmicHandle;
    if (!e.getByLabel(fCosmicLabel, CosmicHandle)) {
        std::cout << "Cannot load any cosmic background. Failing" << std::endl;
        mf::LogWarning("AnalyseEvents_module") << "Cannot load any cosmic background. Failing";
        return;
    }
    for (auto const& cosmic : *CosmicHandle) {
        unsigned const int nBGParticles = cosmic.NParticles();
        for (unsigned int i=0; i<nBGParticles; i++) {
            simb::MCParticle const& particle = cosmic.GetParticle(i);
            fCosmicEnergy[i] = particle.E();
            fCosmicStartX[i] = particle.Position().X();
            fCosmicStartY[i] = particle.Position().Y();
            fCosmicStartZ[i] = particle.Position().Z();
            fCosmicPdgCode[i] = particle.PdgCode();
            fNCosmicParticles++;
        }
    }

    // === Fill the tree ===
    fTree->Fill();
}

void analyse::ParticleEvents::reset() {
    // === Reset particle counters ===
    fNParticles = 0;
    fNMCParticles = 0;
    fNAr39Particles = 0;
    fNAr42Particles = 0;
    fNKr85Particles = 0;
    fNCosmicParticles = 0;
    
    // === Reset energies ===
    fEnergy = 0;
    fMCTotalEnergy = 0;
    
    // === Reset track information ===
    fTotalTrackEnergy = 0;
    fTotalTrackLength = 0;
    
    // === Reset shower information ===
    fTotalShowerEnergy = 0;
    
    // === Reset reco information ===
    for (int i=0; i<NMaxParticles; ++i) {
        // === Track ===
        fTrackLength[i] = 0;
        fTrackStartX[i] = 0;
        fTrackStartY[i] = 0;
        fTrackStartZ[i] = 0;
        fTrackEndX[i] = 0;
        fTrackEndY[i] = 0;
        fTrackEndZ[i] = 0;
        fTrackZenithAngle[i] = 0;
        fTrackAzimuthAngle[i] = 0;
        
        // === Shower ===
        fShowerEnergy[i] = 0;
        fShowerStartX[i] = 0;
        fShowerStartY[i] = 0;
        fShowerStartZ[i] = 0;
        fShowerEndX[i] = 0;
        fShowerEndY[i] = 0;
        fShowerEndZ[i] = 0;
        fShowerLength[i] = 0;
        fShowerDirectionX[i] = 0;
        fShowerDirectionY[i] = 0;
        fShowerDirectionZ[i] = 0;
        fShowerOpenAngle[i] = 0;
    }
    // === Reset MC information ===
    for (int i=0; i<NMaxMCParticles; ++i) {
        fMCParticleEnergy[i] = 0;
        fMCParticleStartPositionX[i] = 0;
        fMCParticleStartPositionY[i] = 0;
        fMCParticleStartPositionZ[i] = 0;
        fMCParticleEndPositionX[i] = 0;
        fMCParticleEndPositionY[i] = 0;
        fMCParticleEndPositionZ[i] = 0;
        fMCParticlePdgCode[i] = 0;
    }
    
    // === Reset BG info ===
    for (int i=0; i<NMaxBGParticles; i++) {
        // Ar39
        fAr39Energy[i] = 0;
        fAr39StartX[i] = 0;
        fAr39StartY[i] = 0;
        fAr39StartZ[i] = 0;
        fAr39PdgCode[i] = 0;
        // Ar42
        fAr42Energy[i] = 0;
        fAr42StartX[i] = 0;
        fAr42StartY[i] = 0;
        fAr42StartZ[i] = 0;
        fAr42PdgCode[i] = 0;
        // Kr85
        fKr85Energy[i] = 0;
        fKr85StartX[i] = 0;
        fKr85StartY[i] = 0;
        fKr85StartZ[i] = 0;
        fKr85PdgCode[i] = 0;
        // Cosmic
        fCosmicEnergy[i] = 0;
        fCosmicStartX[i] = 0;
        fCosmicStartY[i] = 0;
        fCosmicStartZ[i] = 0;
        fCosmicPdgCode[i] = 0;
    }
}

void analyse::ParticleEvents::endJob() {
}

DEFINE_ART_MODULE(analyse::ParticleEvents)