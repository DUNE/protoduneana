////////////////////////////////////////////////////////////////////////
//
// \file CAFMakerPDUNE_module.cc
//
// Chris Marshall's version
// Largely based on historical FDSensOpt/CAFMakerPDUNE_module.cc
// Overhauled by Pierre Granger to adapt it to the new CAF format
//
///////////////////////////////////////////////////////////////////////

#ifndef CAFMakerPDUNE_H
#define CAFMakerPDUNE_H

// Generic C++ includes
#include <iostream>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"

#include "duneanaobj/StandardRecord/Flat/FlatRecord.h"

//#include "Utils/AppInit.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dunereco/CVN/func/InteractionType.h"
#include "dunereco/CVN/func/Result.h"
#include "dunereco/RegCNN/func/RegCNNResult.h"
#include "dunereco/FDSensOpt/FDSensOptData/AngularRecoOutput.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larcore/Geometry/Geometry.h"
#include "nugen/EventGeneratorBase/GENIE/GENIE2ART.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardata/ArtDataHelper/MVAReader.h"


// root
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

// pdg
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

// genie
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/GHEP/GHepParticle.h"


namespace caf {

  class CAFMakerPDUNE : public art::EDAnalyzer {

    public:

      explicit CAFMakerPDUNE(fhicl::ParameterSet const& pset);
      virtual ~CAFMakerPDUNE();
      void beginJob() override;
      void endJob() override;
      void beginSubRun(const art::SubRun& sr) override;
      void endSubRun(const art::SubRun& sr) override;
      void analyze(art::Event const & evt) override;


    private:
      void PreLoadMCParticlesInfo(art::Event const& evt);
      void FillTruthInfo(caf::SRTruthBranch& sr,
                         std::vector<simb::MCTruth> const& mctruth,
                         art::Event const& evt);

      void FillMetaInfo(caf::SRDetectorMeta &meta, art::Event const& evt) const;
      void FillBeamInfo(caf::SRBeamBranch &beam, const art::Event &evt) const;
      void GetMVAResults(caf::SRPFP & output_pfp,
          const recob::PFParticle & pfp, 
          const art::Event &evt, anab::MVAReader<recob::Hit,4> * hitResults, 
          int planeid, 
          bool charge_weighted) const;
      void FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, const art::Event &evt) const;
      void FillRecoParticlesInfo(caf::SRRecoParticlesBranch &recoParticlesBranch, caf::SRFD &fdBranch, const art::Event &evt) const;
      void FillDirectionInfo(caf::SRDirectionBranch &dirBranch, const art::Event &evt) const;
      void FillTruthMatchingAndOverlap(const recob::PFParticle & pfp, const art::Event &evt, std::vector<TrueParticleID> &truth, std::vector<float> &truthOverlap) const;
      void FillPFPMetadata(caf::SRPFP &pfpBranch, art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const;
      // double GetSingleHitsEnergy(art::Event const& evt, int plane) const;

      std::string fMCTruthLabel;
      std::string fPandoraLabel;
      std::string fMVALabel;
      std::string fParticleIDLabel;
      std::string fTrackLabel;
      std::string fShowerLabel;
      std::string fSpacePointLabel;
      std::string fHitLabel;
      std::string fG4Label;
      protoana::ProtoDUNEPFParticleUtils pfpUtil;

      std::map<int, std::tuple<art::Ptr<simb::MCParticle>, bool, int>> fMCParticlesMap; //[tid] = (MCParticle, isPrimary, SRParticle ID)
      uint fNprimaries = 0;
      uint fNsecondaries = 0;

      TTree* fTree = nullptr;
      TTree* fMetaTree = nullptr;

      std::unique_ptr<TFile> fFlatFile;
      TTree* fFlatTree = 0x0; //Ownership will be managed directly by ROOT
      std::unique_ptr<flat::Flat<caf::StandardRecord>> fFlatRecord;

      genie::NtpMCEventRecord *fEventRecord = nullptr;

      double fMetaPOT;
      int fMetaRun, fMetaSubRun, fMetaVersion;

      const geo::Geometry* fGeom;
      std::vector<double> fActiveBounds;
      std::vector<double> fVertexFiducialVolumeCut;

      // calo::CalorimetryAlg fCalorimetryAlg;                    ///< the calorimetry algorithm
      // double fRecombFactor; ///< recombination factor for the isolated hits
      
      // const std::map<simb::Generator_t, caf::Generator> fgenMap = {
      //   {simb::Generator_t::kUnknown, caf::Generator::kUnknownGenerator},
      //   {simb::Generator_t::kGENIE,   caf::Generator::kGENIE},
      //   {simb::Generator_t::kCRY,     caf::Generator::kCRY},
      //   {simb::Generator_t::kGIBUU,   caf::Generator::kGIBUU},
      //   {simb::Generator_t::kNuWro,   caf::Generator::kNuWro},
      //   {simb::Generator_t::kMARLEY,  caf::Generator::kMARLEY},
      //   {simb::Generator_t::kNEUT,    caf::Generator::kNEUT},
      //   {simb::Generator_t::kCORSIKA, caf::Generator::kCORSIKA},
      //   {simb::Generator_t::kGEANT,   caf::Generator::kGEANT}
      // };


  }; // class CAFMakerPDUNE


  //------------------------------------------------------------------------------
  CAFMakerPDUNE::CAFMakerPDUNE(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset),
      fMCTruthLabel(pset.get<std::string>("MCTruthLabel")),
      fPandoraLabel(pset.get< std::string >("PandoraLabel")),
      fMVALabel(pset.get< std::string >("MVALabel")),
      fParticleIDLabel(pset.get< std::string >("ParticleIDLabel")),
      fTrackLabel(pset.get< std::string >("TrackLabel")),
      fShowerLabel(pset.get< std::string >("ShowerLabel")),
      fSpacePointLabel(pset.get< std::string >("SpacePointLabel")),
      fHitLabel(pset.get< std::string >("HitLabel")),
      fG4Label(pset.get< std::string >("G4Label")),
      fGeom(&*art::ServiceHandle<geo::Geometry>()) //,
      // fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
      // fRecombFactor(pset.get<double>("RecombFactor"))
  {

    if(pset.get<bool>("CreateFlatCAF", false)){
      // LZ4 is the fastest format to decompress. I get 3x faster loading with
      // this compared to the default, and the files are only slightly larger.
      fFlatFile = std::make_unique<TFile>("flatcaf.root", "RECREATE", "",
                            ROOT::CompressionSettings(ROOT::kLZ4, 1));
    }

  }

  //------------------------------------------------------------------------------
  caf::CAFMakerPDUNE::~CAFMakerPDUNE()
  {
  }

  //------------------------------------------------------------------------------
  void CAFMakerPDUNE::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("cafTree", "cafTree");

    // Create the branch. We will update the address before we write the tree
    caf::StandardRecord* rec = 0;
    fTree->Branch("rec", "caf::StandardRecord", &rec);

    fMetaTree = tfs->make<TTree>("meta", "meta");

    fMetaTree->Branch("pot", &fMetaPOT, "pot/D");
    fMetaTree->Branch("run", &fMetaRun, "run/I");
    fMetaTree->Branch("subrun", &fMetaSubRun, "subrun/I");
    fMetaTree->Branch("version", &fMetaVersion, "version/I");

    fMetaPOT = 0.;
    fMetaVersion = 1;

    if(fFlatFile){
      fFlatFile->cd();
      fFlatTree = new TTree("cafTree", "cafTree");

      fFlatRecord = std::make_unique<flat::Flat<caf::StandardRecord>>(fFlatTree, "rec", "", nullptr);
    }

  }

  //------------------------------------------------------------------------------

  void CAFMakerPDUNE::PreLoadMCParticlesInfo(art::Event const& evt)
  {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    //Preloading the MCParticles info to be able to access them later
    std::vector<art::Ptr< simb::MCParticle>> mcparticles = dune_ana::DUNEAnaEventUtils::GetMCParticles(evt, fG4Label);

    //Creating a map of the MC particles to easily access them by their TrackId
    fMCParticlesMap.clear();
    //Not the most efficient way, but I prefer good readibility over performance here
    fNprimaries = 0;
    fNsecondaries = 0;
    for (art::Ptr<simb::MCParticle> const& mcpart: mcparticles) {

      auto origin = pi_serv->TrackIdToMCTruth_P(mcpart->TrackId())->Origin();
      if (origin != 4) continue; //Skip non-beam/non-decay for now 

      // std::cout << mcpart->PdgCode() << " " << mcpart->TrackId() << " " << mcpart->Process() << std::endl;

      // if(isPrimary) {
        fMCParticlesMap[mcpart->TrackId()] = {mcpart, true, 0};
        // fNprimaries++;
      // } else {
      //   fMCParticlesMap[mcpart->TrackId()] = {mcpart, false, fNsecondaries};
      //   fNsecondaries++;
      // } 
    }
  }

  //------------------------------------------------------------------------------

  void CAFMakerPDUNE::FillTruthInfo(caf::SRTruthBranch& truthBranch,
                              std::vector<simb::MCTruth> const& mctruth,
                              art::Event const& evt) {
    art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
    const sim::ParticleList & plist = pi_serv->ParticleList();

  
    if (mctruth.size() != 1)
      throw std::runtime_error("MCTruth vector size is not 1 -- this is unexpected");
    
    auto this_mctruth = mctruth[0];
    caf::SRTrueInteraction inter;
    // std::cout << this_mctruth.Origin() << std::endl;
    //   //Looping on the particles in the MC truth to get the beam particles only
    std::deque<unsigned int> secondaries_to_add;
    for (int p = 0; p < this_mctruth.NParticles(); p++) {
        const simb::MCParticle &mcpart = this_mctruth.GetParticle(p);

        caf::SRTrueParticle part;
        int pdg = mcpart.PdgCode();
        part.pdg = pdg;
        part.G4ID = abs(mcpart.TrackId());
        part.time = mcpart.T();
        part.p = caf::SRLorentzVector(mcpart.Momentum());
        part.start_pos = caf::SRVector3D(mcpart.Position().Vect());
        part.end_pos = caf::SRVector3D(mcpart.EndPosition().Vect());
        part.parent = mcpart.Mother();
        // std::cout << "\t" << part.G4ID << " " << " " << pdg << " " << mcpart.Process() << std::endl;

        const simb::MCParticle * the_g4_part = nullptr;

        for (auto const g4_part : plist) {
          if ((abs(mcpart.TrackId()) == g4_part.second->TrackId())){
            // std::cout << "FOUND IN G4" << std::endl;
            the_g4_part = g4_part.second;
            break;
          }
        }
        if (!the_g4_part) continue;
        for(int i = 0; i < the_g4_part->NumberDaughters(); i++){
          int daughter = the_g4_part->Daughter(i);
          part.daughters.push_back(daughter);
          secondaries_to_add.push_back(daughter);
        }

        inter.prim.push_back(std::move(part));
    }

    while (secondaries_to_add.size()) {
        auto this_sec_ID = secondaries_to_add.front();
        if (fMCParticlesMap.count(this_sec_ID) > 0) {
            auto mcpart = std::get<0>(fMCParticlesMap[this_sec_ID]);
            caf::SRTrueParticle part;
            int pdg = mcpart->PdgCode();
            part.pdg = pdg;
            part.G4ID = abs(mcpart->TrackId());
            part.time = mcpart->T();
            part.p = caf::SRLorentzVector(mcpart->Momentum());
            part.start_pos = caf::SRVector3D(mcpart->Position().Vect());
            part.end_pos = caf::SRVector3D(mcpart->EndPosition().Vect());
            part.parent = mcpart->Mother();

            for(int i = 0; i < mcpart->NumberDaughters(); i++){
                int daughter = mcpart->Daughter(i);
                part.daughters.push_back(daughter);
                secondaries_to_add.push_back(daughter);
            }
            inter.sec.push_back(std::move(part));
        }
        secondaries_to_add.pop_front();
    }
    truthBranch.nu.push_back(std::move(inter));

  }

  //------------------------------------------------------------------------------

  void CAFMakerPDUNE::FillPFPMetadata(caf::SRPFP &pfpBranch, art::Ptr<recob::PFParticle> const& pfp, const art::Event &evt) const {
    art::Ptr<larpandoraobj::PFParticleMetadata> metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfp, evt, fPandoraLabel);
      if(metadata.isNull()){
        mf::LogWarning("CAFMakerPDUNE") << "No metadata found for PFP with ID " << pfp->Self();
        return;
      }

      std::map<std::string, float> properties = metadata->GetPropertiesMap();

      std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fPandoraLabel);
      std::vector<art::Ptr<recob::Hit>> hits_U = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, 0);
      std::vector<art::Ptr<recob::Hit>> hits_V = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, 1);
      std::vector<art::Ptr<recob::Hit>> hits_W = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, 2);
      std::vector<art::Ptr<recob::SpacePoint>> sps = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfp, evt, fSpacePointLabel);

      pfpBranch.nhits_U = hits_U.size();
      pfpBranch.nhits_V = hits_V.size();
      pfpBranch.nhits_W = hits_W.size();
      pfpBranch.nhits_3D = sps.size();

      std::map<std::string, float*> property_mapping = {
        // {"LArPfoHierarchyFeatureTool_DaughterParentHitRatio", &pfpBranch.daughter_parent_hit_ratio},
        // {"LArPfoHierarchyFeatureTool_NDaughterHits3D", &pfpBranch.ndaughters_hit_3d},
        // {"LArThreeDChargeFeatureTool_EndFraction", &pfpBranch.charge_end_fraction},
        // {"LArThreeDChargeFeatureTool_FractionalSpread", &pfpBranch.charge_fractional_spread},
        // {"LArThreeDLinearFitFeatureTool_DiffStraightLineMean", &pfpBranch.diff_straight_line_mean},
        // {"LArThreeDLinearFitFeatureTool_Length", &pfpBranch.line_length},
        // {"LArThreeDLinearFitFeatureTool_MaxFitGapLength", &pfpBranch.max_fit_gap_length},
        // {"LArThreeDLinearFitFeatureTool_SlidingLinearFitRMS", &pfpBranch.sliding_linear_fit_rms},
        // {"LArThreeDOpeningAngleFeatureTool_AngleDiff", &pfpBranch.angle_diff_3d},
        // {"LArThreeDPCAFeatureTool_SecondaryPCARatio", &pfpBranch.secondary_pca_ratio},
        // {"LArThreeDPCAFeatureTool_TertiaryPCARatio", &pfpBranch.tertiary_pca_ratio},
        // {"LArThreeDVertexDistanceFeatureTool_VertexDistance", &pfpBranch.vertex_distance},
        // {"TrackScore", &pfpBranch.track_score},
      };

      for(const auto& [property_name, property_value] : property_mapping) {
        if(properties.find(property_name) != properties.end()) {
          *(property_value) = properties[property_name];
        } else {
          mf::LogWarning("CAFMakerPDUNE") << "Unknown property '" << property_name << "' for PFP with ID " << pfp->Self();
        }
      }

      //TODO -- configure debug
      // for (const auto & [prop, val] : properties) {
      //   std::cout << prop << " " << val << std::endl;
      // }

  }


  //------------------------------------------------------------------------------

  void CAFMakerPDUNE::FillRecoInfo(caf::SRCommonRecoBranch &recoBranch, caf::SRFD &fdBranch, const art::Event &evt) const {
    SRInteractionBranch &ixn = recoBranch.ixn;

    //Only filling with Pandora Reco for the moment
    std::vector<SRInteraction> &pandora = ixn.pandora;
    
    lar_pandora::PFParticleVector particleVector;
    lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraLabel, particleVector);
    // lar_pandora::VertexVector vertexVector;
    // lar_pandora::PFParticlesToVertices particlesToVertices;
    // lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraLabel, vertexVector, particlesToVertices);

    SRInteraction reco;

    reco.vtx = SRVector3D(-999, -999, -999); //Setting an unambiguous default value if no vertex is found

    // //Retrieving the reco vertex
    // lar_pandora::PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
    // if (particlesToVertices.end() != vIter) {
    //   const lar_pandora::VertexVector &vertexVector = vIter->second;
    //   if (vertexVector.size() == 1) {
    //     const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
    //     double xyz[3] = {0.0, 0.0, 0.0} ;
    //     vertex->XYZ(xyz);
    //     reco.vtx = SRVector3D(xyz[0], xyz[1], xyz[2]);
    //   }
    // }

    //List of reconstructed particles
    SRRecoParticlesBranch &part = reco.part;
    FillRecoParticlesInfo(part, fdBranch, evt);

    //Assuming a single TrueInteraction for now. TODO: Change this if several interactions end up being simulated in the same event
    reco.truth = {0};
    reco.truthOverlap = {1.};

    pandora.emplace_back(reco);

    ixn.npandora = pandora.size();
    ixn.ndlp = ixn.dlp.size();
    std::cout << "Done reco" << std::endl;
  }


  //------------------------------------------------------------------------------


  void CAFMakerPDUNE::FillTruthMatchingAndOverlap(
    const recob::PFParticle & pfp,
    const art::Event &evt, std::vector<TrueParticleID> &truth,
    std::vector<float> &truthOverlap) const {
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    //First getting all the hits belonging to that PFP
    // std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fPandoraLabel);
    const auto & hits = pfpUtil.GetPFParticleHits_Ptrs(pfp, evt, fPandoraLabel);
    TruthMatchUtils::IDToEDepositMap idToEDepositMap;
    for (const art::Ptr<recob::Hit>& pHit : hits){
      TruthMatchUtils::FillG4IDToEnergyDepositMap(idToEDepositMap, clockData, pHit, true);
    }

    float totalEDeposit = 0;
    for (const auto& [id, eDeposit] : idToEDepositMap) {
      totalEDeposit += eDeposit;
    }

    if (totalEDeposit <= 0) {
      mf::LogWarning("CAFMakerPDUNE") << "No energy deposit found for PFP with ID " << pfp.Self() << ". Skipping truth matching.";
      return;
    }

    for (const auto& [id, eDeposit] : idToEDepositMap) {
      if (fMCParticlesMap.count(id) == 0) {
        mf::LogWarning("CAFMakerPDUNE") << "No MCParticle found with ID " << id << " for PFP with ID " << pfp.Self() << ". Skipping.";
        continue;
      }

      art::Ptr<simb::MCParticle> mcpart;
      bool isPrimary;
      int srID;
      std::tie(mcpart, isPrimary, srID) = fMCParticlesMap.at(id);

      TrueParticleID truePart;
      truePart.type = isPrimary ? TrueParticleID::kPrimary : TrueParticleID::kSecondary;
      truePart.ixn = 0; //Assuming a single interaction for now
      truePart.part = srID;

      truth.push_back(truePart);
      truthOverlap.push_back(eDeposit / totalEDeposit);
    }

  }

  //------------------------------------------------------------------------------

  // double CAFMakerPDUNE::GetSingleHitsEnergy(art::Event const& evt, int plane) const{
  //   std::vector<art::Ptr<recob::Hit>> hits = dune_ana::DUNEAnaEventUtils::GetHits(evt, fHitLabel);

  //   std::vector<art::Ptr<recob::Hit>> collection_plane_hits = dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(hits, plane);

  //   const art::FindManyP<recob::SpacePoint> sp_assoc(hits, evt, fSpacePointLabel);

  //   if (!sp_assoc.isValid()) {
  //     mf::LogWarning("CAFMakerPDUNE") << "No space points found with label '" << fSpacePointLabel << "'";
  //     return 0;
  //   }

  //   auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  //   auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

  //   double charge = 0;

  //   for (uint i = 0; i < collection_plane_hits.size(); i++){
  //     std::vector<art::Ptr<recob::SpacePoint>> matching_sps = sp_assoc.at(collection_plane_hits[i].key());
  //     if(!matching_sps.empty()){ // If there are some matched spacepoints, this means that the hit is associated to some PFP and we don't want it here
  //       continue;
  //     }
  //     // Get the energy deposited in the hit
  //     charge += dune_ana::DUNEAnaHitUtils::LifetimeCorrection(clockData, detProp, collection_plane_hits[i])*collection_plane_hits[i]->Integral();
  //   }

  //   return fCalorimetryAlg.ElectronsFromADCArea(charge,2)*1./fRecombFactor/util::kGeVToElectrons;

  // }

  // //------------------------------------------------------------------------------
 
 
  void CAFMakerPDUNE::beginSubRun(const art::SubRun& sr)
  {
    fMetaRun = sr.id().subRun();
    fMetaSubRun = sr.id().run();

  }

  //------------------------------------------------------------------------------

  void CAFMakerPDUNE::FillMetaInfo(caf::SRDetectorMeta &meta, const art::Event &evt) const
  {
    meta.enabled = true;
    meta.run = evt.id().run();
    meta.subrun = evt.id().subRun();
    meta.event = evt.id().event();
    meta.subevt = 0; //Hardcoded to 0, only makes sense in ND where multiple interactions can occur in the same event

    //Nothing is filled about the trigger for the moment
  }

  //------------------------------------------------------------------------------

  void CAFMakerPDUNE::FillBeamInfo(caf::SRBeamBranch &beam, const art::Event &evt) const
  {
    //This part will only be relevant when working on real data with real beam.
    beam.ismc = true; //Hardcoded to true at the moment.
  }

  //------------------------------------------------------------------------------

  void CAFMakerPDUNE::GetMVAResults(
    caf::SRPFP & output_pfp,
    const recob::PFParticle & pfp, 
    const art::Event &evt, anab::MVAReader<recob::Hit,4> * hitResults, 
    int planeid, 
    bool charge_weighted) const {

    if (planeid < -1 || planeid > 2) {
      std::stringstream ss;
      ss << "Unknown planeid (" << planeid << ") provided to GetMVAResults";
      throw std::runtime_error(
        ss.str()
      );
    }

    //First getting all the hits belonging to that PFP
    const auto & hits = (
      planeid == -1 ?
      pfpUtil.GetPFParticleHits_Ptrs(pfp, evt, fPandoraLabel) :
      pfpUtil.GetPFParticleHitsFromPlane_Ptrs(pfp, evt, fPandoraLabel, planeid)
    );

    float denom = 0.;
    output_pfp.cnn_stem_scores.charge_weighted = charge_weighted;
    output_pfp.cnn_stem_scores.plane_ID = planeid;
    for (const auto & hit : hits){
      auto output = hitResults->getOutput(hit);

      float scale = (charge_weighted ? hit->Integral() : 1.);

      output_pfp.cnn_stem_scores.track_score += scale*output[hitResults->getIndex("track")];
      output_pfp.cnn_stem_scores.shower_score += scale*output[hitResults->getIndex("em")];
      output_pfp.cnn_stem_scores.empty_score += scale*output[hitResults->getIndex("none")];
      output_pfp.cnn_stem_scores.michel_score += scale*output[hitResults->getIndex("michel")];
      denom += scale;
    }

    if (denom > 0.) {
      output_pfp.cnn_stem_scores /= denom;
    }
    else {
      output_pfp.cnn_stem_scores.track_score = output_pfp.NaN;
      output_pfp.cnn_stem_scores.shower_score = output_pfp.NaN;
      output_pfp.cnn_stem_scores.empty_score = output_pfp.NaN;
      output_pfp.cnn_stem_scores.michel_score = output_pfp.NaN;
    }
  }

  //------------------------------------------------------------------------------

  void CAFMakerPDUNE::FillRecoParticlesInfo(caf::SRRecoParticlesBranch &recoParticlesBranch, caf::SRFD &fdBranch, const art::Event &evt) const
  {

    //Getting Ar density in g/cm3
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
    double lar_density = detProp.Density();

    //Getting all the PFParticles from the event
    // lar_pandora::PFParticleVector particleVector;
    // lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraLabel, particleVector);

    //Get all PFParticles in the event
    auto pfpVec
        = evt.getValidHandle<std::vector<recob::PFParticle>>(fPandoraLabel);
    
    anab::MVAReader<recob::Hit,4> * hitResults = new anab::MVAReader<recob::Hit, 4>(evt, fMVALabel/*"emtrkmichelid:emtrkmichel"*/ );
    
    const auto sliceMap = pfpUtil.GetPFParticleSliceMap(evt, fPandoraLabel);
    std::cout << sliceMap.size() << " Slices" << std::endl;
    int n_test_beam_slices = 0;
    for (const auto & [slice, particle_vector] : sliceMap) {
      // std::cout << "NParticles: " << particle_vector.size() << std::endl;
      if (particle_vector.size() == 0)
        throw std::runtime_error("Found empty particle vector in a slice");
      // if (particle_vector.size() == 0) continue;

      //Creating the FD interaction record where we are going to save the tracks/showers in parallel to the PFPs objects
      caf::SRFDInt fdIxn;

      bool is_test_beam = false;
      auto slice_particle = particle_vector[0];
      is_test_beam = pfpUtil.IsBeamParticle(*slice_particle, evt, fPandoraLabel);
      if(!is_test_beam) continue;

      std::deque<const recob::PFParticle*> to_add{slice_particle};

      while (to_add.size() > 0) {

        const auto * particle = to_add.front();
        to_add.pop_front();
        
        // std::cout << "Found test beam slice. N particles in slice: " << slice.second.size() << std::endl;

        //Creating the particle record for this PFP
        caf::SRRecoParticle particle_record;
        FillTruthMatchingAndOverlap(*particle, evt, particle_record.truth, particle_record.truthOverlap);

        //Getting the track and shower objects associated to the PFP
        const recob::Track * track = 0x0;
        try {
          track = pfpUtil.GetPFParticleTrack(*(particle), evt, fPandoraLabel, fTrackLabel);
        }
        catch (const cet::exception &e) {
          mf::LogInfo("CAFMakerPDUNE") << "No associated track object. Moving on";
        }
        
        const recob::Shower * shower = 0x0;
        try {
          shower = pfpUtil.GetPFParticleShower(*(particle), evt, fPandoraLabel, fShowerLabel);
        }
        catch (const cet::exception &e) {
          mf::LogInfo("CAFMakerPDUNE") << "No associated shower object. Moving on";
        }

        //Seeing which option Pandora prefers
        // bool isTrack = lar_pandora::LArPandoraHelper::IsTrack(particle);

        //For every PFP we create a track and a shower object and save it, independently of the existence of a track object to keep the PFP/Track/Shower parallel indexing
        SRTrack srtrack;
        SRShower srshower;

        //This variable will be updated correctly during the recob::Track processing and will be used to fill the reco particle energy method if isTrack.
        // caf::PartEMethod trackErecoMethod = caf::PartEMethod::kUnknownMethod;

        if(track){
          // std::cout << "Getting track stuff" << std::endl;
          srtrack.start.SetX(track->Start().X());
          srtrack.start.SetY(track->Start().Y());
          srtrack.start.SetZ(track->Start().Z());

          srtrack.end.SetX(track->End().X());
          srtrack.end.SetY(track->End().Y());
          srtrack.end.SetZ(track->End().Z());

          srtrack.dir.SetX(track->StartDirection().X());
          srtrack.dir.SetY(track->StartDirection().Y());
          srtrack.dir.SetZ(track->StartDirection().Z());

          srtrack.enddir.SetX(track->EndDirection().X());
          srtrack.enddir.SetY(track->EndDirection().Y());
          srtrack.enddir.SetZ(track->EndDirection().Z());

          //srtrack.qual TODO: Not sure we have anything relevant to put on the FD side for this at the moment

          srtrack.len_gcm2 = track->Length() * lar_density; //Length in g/cm2
          srtrack.len_cm = track->Length();

          //TODO: I would prefer to use some unified module that the user can setup and that will decide how to compute the energy rather than making a specific choice here
          //Putting Evis as placeholder to not confuse the user too much
          // trackErecoMethod = caf::PartEMethod::kCalorimetry; //Using the visible energy of the PFP

          //Truth matching already filled at the PFP level, no need to do it again here
          srtrack.truth = particle_record.truth;
          srtrack.truthOverlap = particle_record.truthOverlap;

        }

        if(shower){
          // std::cout << "Getting shower stuff" << std::endl;
          //Filling the shower information
          srshower.start.SetX(shower->ShowerStart().X());
          srshower.start.SetY(shower->ShowerStart().Y());
          srshower.start.SetZ(shower->ShowerStart().Z());

          srshower.direction.SetX(shower->Direction().X());
          srshower.direction.SetY(shower->Direction().Y());
          srshower.direction.SetZ(shower->Direction().Z());
          
          //Truth matching already filled at the PFP level, no need to do it again here
          srshower.truth = particle_record.truth;        
          srshower.truthOverlap = particle_record.truthOverlap;
        }

        // if(isTrack){
        //   if(track){ //I hope this condition is always fullfilled is the particle is tagged at track, but who knows...
        //     // std::cout << "Getting track stuff2" << std::endl;
        //     particle_record.start = SRVector3D(track->Start().X(), track->Start().Y(), track->Start().Z());
        //     particle_record.end = SRVector3D(track->End().X(), track->End().Y(), track->End().Z());
        //     particle_record.E = srtrack.E;
        //     particle_record.E_method = trackErecoMethod;
        //   }
        //   particle_record.origRecoObjType = caf::RecoObjType::kTrack;
        // }
        // else{
        //   if(shower){ //I hope this condition is always fullfilled is the particle is tagged at shower, but who knows...
        //     // std::cout << "Getting shower stuff2" << std::endl;
        //     particle_record.start = SRVector3D(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
        //     //Only filling the start, no defined end for a shower

        //     particle_record.E = srshower.Evis; //Using the visible energy of the PFP
        //     particle_record.E_method = caf::PartEMethod::kCalorimetry;
            
        //   }
        //   particle_record.origRecoObjType = caf::RecoObjType::kShower;
        // }

        //Saving the track record
        fdIxn.tracks.push_back(std::move(srtrack));
        fdIxn.ntracks++;

        //Saving the shower record
        fdIxn.showers.push_back(std::move(srshower));
        fdIxn.nshowers++;

        //Also saving PFP metadata
        caf::SRPFP pfp_metarecord;
        GetMVAResults(
          pfp_metarecord,
          *particle,
          evt,
          hitResults, 
          2, //TODO -- make configurable
          true
        );
        pfp_metarecord.parent = particle->Parent();
        pfp_metarecord.daughters.insert(
          pfp_metarecord.daughters.begin(),
          particle->Daughters().begin(),particle->Daughters().end()
        );
        pfp_metarecord.truth = particle_record.truth;
        pfp_metarecord.truthOverlap = particle_record.truthOverlap;
        // FillPFPMetadata(pfp_metarecord, particle, evt);
        fdIxn.pfps.push_back(std::move(pfp_metarecord));
        fdIxn.npfps++;
          
        //Saving the particle record for this PFP
        recoParticlesBranch.pandora.push_back(std::move(particle_record));
        recoParticlesBranch.npandora++;

        //Iterate over all daughters and add to the queue
        for (size_t daughterID : particle->Daughters()) {
          to_add.push_back(&(pfpVec->at(daughterID)));
          std::cout << "Added daughter " << daughterID << std::endl;
          // const recob::PFParticle * daughterPFP = &(pfpVec->at(daughterID));
        }
      }

      if (!is_test_beam) continue;
      ++n_test_beam_slices;
      fdBranch.pandora.push_back(std::move(fdIxn));
      fdBranch.npandora++;
    }
    
    std::cout << "N test beam slices " << n_test_beam_slices << std::endl;

    // //Iterating on all the PFParticles to fill the reco particles
    // int ntest = 0;
    // for (unsigned int n = 0; n < particleVector.size(); ++n) {
    //   std::cout << "Processing part " << n << "/" << particleVector.size() << std::endl;
    //   const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
    //   const auto * part_ptr = particle.get();
    //   // if(particle->Self() == nuID){ //Skipping the neutrino that is not a "real" reco particle
    //   //   continue;
    //   // }

    //   // std::cout << "IsBeam: " << pfpUtil.IsBeamParticle(*(part_ptr), evt, fPandoraLabel) <<
    //   //              " TestBeamScore " << pfpUtil.GetBeamCosmicScore(*part_ptr, evt, fPandoraLabel) << std::endl;
    //   if (pfpUtil.IsBeamParticle(*(part_ptr), evt, fPandoraLabel))
    //     ++ntest;


      

    //   //Getting the track and shower objects associated to the PFP
    //   const recob::Track * track = 0x0;
    //   try {
    //     track = pfpUtil.GetPFParticleTrack(*(part_ptr), evt, fPandoraLabel, fTrackLabel);
    //   }
    //   catch (const cet::exception &e) {
    //     mf::LogInfo("CAFMakerPDUNE") << "No associated track object. Moving on";
    //   }
      
    //   const recob::Shower * shower = 0x0;
    //   try {
    //     shower = pfpUtil.GetPFParticleShower(*(part_ptr), evt, fPandoraLabel, fShowerLabel);
    //   }
    //   catch (const cet::exception &e) {
    //     mf::LogInfo("CAFMakerPDUNE") << "No associated shower object. Moving on";
    //   }

    //   //Seeing which option Pandora prefers
    //   bool isTrack = lar_pandora::LArPandoraHelper::IsTrack(particle);

    //   //For every PFP we create a track and a shower object and save it, independently of the existence of a track object to keep the PFP/Track/Shower parallel indexing
    //   SRTrack srtrack;
    //   SRShower srshower;

    //   //This variable will be updated correctly during the recob::Track processing and will be used to fill the reco particle energy method if isTrack.
    //   caf::PartEMethod trackErecoMethod = caf::PartEMethod::kUnknownMethod;

    //   if(track){
    //     // std::cout << "Getting track stuff" << std::endl;
    //     srtrack.start.SetX(track->Start().X());
    //     srtrack.start.SetY(track->Start().Y());
    //     srtrack.start.SetZ(track->Start().Z());

    //     srtrack.end.SetX(track->End().X());
    //     srtrack.end.SetY(track->End().Y());
    //     srtrack.end.SetZ(track->End().Z());

    //     srtrack.dir.SetX(track->StartDirection().X());
    //     srtrack.dir.SetY(track->StartDirection().Y());
    //     srtrack.dir.SetZ(track->StartDirection().Z());

    //     srtrack.enddir.SetX(track->EndDirection().X());
    //     srtrack.enddir.SetY(track->EndDirection().Y());
    //     srtrack.enddir.SetZ(track->EndDirection().Z());

    //     //srtrack.qual TODO: Not sure we have anything relevant to put on the FD side for this at the moment

    //     srtrack.len_gcm2 = track->Length() * lar_density; //Length in g/cm2
    //     srtrack.len_cm = track->Length();

    //     //TODO: I would prefer to use some unified module that the user can setup and that will decide how to compute the energy rather than making a specific choice here
    //     //Putting Evis as placeholder to not confuse the user too much
    //     trackErecoMethod = caf::PartEMethod::kCalorimetry; //Using the visible energy of the PFP

    //     //Truth matching already filled at the PFP level, no need to do it again here
    //     srtrack.truth = particle_record.truth;
    //     srtrack.truthOverlap = particle_record.truthOverlap;

    //   }

    //   if(shower){
    //     // std::cout << "Getting shower stuff" << std::endl;
    //     //Filling the shower information
    //     srshower.start.SetX(shower->ShowerStart().X());
    //     srshower.start.SetY(shower->ShowerStart().Y());
    //     srshower.start.SetZ(shower->ShowerStart().Z());

    //     srshower.direction.SetX(shower->Direction().X());
    //     srshower.direction.SetY(shower->Direction().Y());
    //     srshower.direction.SetZ(shower->Direction().Z());
        
    //     //Truth matching already filled at the PFP level, no need to do it again here
    //     srshower.truth = particle_record.truth;        
    //     srshower.truthOverlap = particle_record.truthOverlap;
    //   }

    //   if(isTrack){
    //     if(track){ //I hope this condition is always fullfilled is the particle is tagged at track, but who knows...
    //       // std::cout << "Getting track stuff2" << std::endl;
    //       particle_record.start = SRVector3D(track->Start().X(), track->Start().Y(), track->Start().Z());
    //       particle_record.end = SRVector3D(track->End().X(), track->End().Y(), track->End().Z());
    //       particle_record.E = srtrack.E;
    //       particle_record.E_method = trackErecoMethod;
    //     }
    //     particle_record.origRecoObjType = caf::RecoObjType::kTrack;
    //   }
    //   else{
    //     if(shower){ //I hope this condition is always fullfilled is the particle is tagged at shower, but who knows...
    //       // std::cout << "Getting shower stuff2" << std::endl;
    //       particle_record.start = SRVector3D(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    //       //Only filling the start, no defined end for a shower

    //       particle_record.E = srshower.Evis; //Using the visible energy of the PFP
    //       particle_record.E_method = caf::PartEMethod::kCalorimetry;
          
    //     }
    //     particle_record.origRecoObjType = caf::RecoObjType::kShower;
    //   }

    //   //Saving the track record
    //   fdIxn.tracks.push_back(std::move(srtrack));
    //   fdIxn.ntracks++;

    //   //Saving the shower record
    //   fdIxn.showers.push_back(std::move(srshower));
    //   fdIxn.nshowers++;

    //   //Also saving PFP metadata
    //   caf::SRPFP pfp_metarecord;
    //   FillPFPMetadata(pfp_metarecord, particle, evt);
    //   fdIxn.pfps.push_back(std::move(pfp_metarecord));
    //   fdIxn.npfps++;
        
    //   //Saving the particle record for this PFP
    //   recoParticlesBranch.pandora.push_back(std::move(particle_record));
    //   recoParticlesBranch.npandora++;

    // }
    // std::cout << "N test beam candidates: " << ntest << std::endl;

    // //Saving the FD interaction record
    // fdBranch.pandora.push_back(std::move(fdIxn));
    // fdBranch.npandora++;
    // std::cout << "Done" << std::endl;
    // //Adding some extra record with all the leftover hits not associated to any particle
    // caf::SRRecoParticle single_hits;
    // single_hits.primary = false;
    // single_hits.pdg = 0; //Not a real particle
    // single_hits.tgtA = 40; //Interaction on Ar40.
    // // single_hits.E = GetSingleHitsEnergy(evt, 2); //Using the collection plane for now
    // single_hits.origRecoObjType = caf::RecoObjType::kHitCollection;

    // recoParticlesBranch.pandora.push_back(std::move(single_hits));
    // recoParticlesBranch.npandora++;

    // particle_record.origRecoObjType
  }


  //------------------------------------------------------------------------------


  void CAFMakerPDUNE::analyze(art::Event const & evt)
  {
    caf::StandardRecord sr;
    caf::StandardRecord* psr = &sr;

    PreLoadMCParticlesInfo(evt);
    

    if(fTree){
      fTree->SetBranchAddress("rec", &psr);
    }

    std::string geoName = fGeom->DetectorName();

    mf::LogInfo("CAFMakerPDUNE") << "Geo name is: " << geoName;

    SRDetectorMeta *detector;
    SRFD *fdBranch;

    //TODO -- SWITCH TO PROTODUNE
    // if(geoName.find("dunevd10kt") != std::string::npos){ // ProtoDUNE Strings?
      detector = &(sr.meta.pd_hd);
      fdBranch = &(sr.fd.pd_hd);
    // }
    
    FillMetaInfo(*detector, evt);

    // FillBeamInfo(sr.beam, evt);
    art::Handle<std::vector<simb::MCTruth>> mct = evt.getHandle< std::vector<simb::MCTruth> >(fMCTruthLabel);
    if ( !mct ) {
      mf::LogWarning("CAFMakerPDUNE") << "No MCTruth. SRTruthBranch will be empty!";
    }
    else {
      FillTruthInfo(sr.mc, *mct, evt);
    }

    FillRecoInfo(sr.common, *fdBranch, evt);

    if(fTree){
      fTree->Fill();
    }

    if(fFlatTree){
      fFlatRecord->Clear();
      fFlatRecord->Fill(sr);
      fFlatTree->Fill();
    }
    std::cout << "Done event" << std::endl;
  }

  //------------------------------------------------------------------------------

  //------------------------------------------------------------------------------
  void CAFMakerPDUNE::endSubRun(const art::SubRun& sr){
  }

  void CAFMakerPDUNE::endJob()
  {
    fMetaTree->Fill();

    if(fFlatFile){
      fFlatFile->cd();
      fFlatTree->Write();
      fMetaTree->CloneTree()->Write();
      fFlatFile->Close();
    }

    delete fEventRecord; //Making this a unique_pointer requires too many circonvolutions because of TTree->Branch requiring a pointer to a pointer

  }

  DEFINE_ART_MODULE(CAFMakerPDUNE)

} // namespace caf

#endif // CAFMakerPDUNE_H
