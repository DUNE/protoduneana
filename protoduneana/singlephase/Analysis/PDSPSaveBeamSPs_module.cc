////////////////////////////////////////////////////////////////////////
// Class:       PDSPSaveBeamSPs
// Plugin Type: analyzer (Unknown Unknown)
// File:        PDSPSaveBeamSPs_module.cc
//
// Generated at Fri Aug 18 18:23:20 2023 by Jacob Calcutt using cetskelgen
// from  version .
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

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "hep_hpc/hdf5/File.hpp"
#include "TTree.h"
#include "art_root_io/TFileService.h"

namespace pduneana {
  class PDSPSaveBeamSPs;

  struct SPNode {
    SPNode(double x, double y, double z, double c0, double c1, double c2,
           double rms0, double rms1, double rms2)
     : fX(x), fY(y), fZ(z), fC0(c0), fC1(c1), fC2(c2),
       fRMS0(rms0), fRMS1(rms1), fRMS2(rms2) {};
    double fX, fY, fZ, fC0, fC1, fC2, fRMS0, fRMS1, fRMS2;
  };

  struct SPEdge {
    //SPEdge(size_t i, size_t j, std::vector<SPNode> & node_vec)
    SPEdge(size_t i, size_t j, SPNode & node_i, SPNode & node_j)
      : fI(i), fJ(j), fNodeI(node_i), fNodeJ(node_j) {
      fDiffX = fNodeI.fX - fNodeJ.fX;
      fDiffY = fNodeI.fY - fNodeJ.fY;
      fDiffZ = fNodeI.fZ - fNodeJ.fZ;

      fMag = sqrt(fDiffX*fDiffX + fDiffY*fDiffY + fDiffZ*fDiffZ);

      fDiffC0 = fNodeI.fC0 - fNodeJ.fC0;
      fDiffC1 = fNodeI.fC1 - fNodeJ.fC1;
      fDiffC2 = fNodeI.fC2 - fNodeJ.fC2;

      fDiffOP = {
        fDiffX*fDiffX, fDiffX*fDiffY, fDiffX*fDiffZ,
        fDiffY*fDiffY, fDiffY*fDiffZ, fDiffZ*fDiffZ
      };

      fDiffCOP = {
        fDiffC0*fDiffC0, fDiffC0*fDiffC1, fDiffC0*fDiffC2,
        fDiffC1*fDiffC1, fDiffC1*fDiffC2, fDiffC2*fDiffC2
      };

      fAveCharge = {
        (fNodeI.fC0 + fNodeJ.fC0)/2.,
        (fNodeI.fC1 + fNodeJ.fC1)/2.,
        (fNodeI.fC2 + fNodeJ.fC2)/2.
      };
    };

    size_t fI, fJ;
    SPNode fNodeI, fNodeJ;
    double fDiffX, fDiffY, fDiffZ, fMag, fDiffC0, fDiffC1, fDiffC2;
    std::array<double, 6> fDiffOP, fDiffCOP;
    std::array<double, 3> fAveCharge;
  };
}


class pduneana::PDSPSaveBeamSPs : public art::EDAnalyzer {
public:
  explicit PDSPSaveBeamSPs(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPSaveBeamSPs(PDSPSaveBeamSPs const&) = delete;
  PDSPSaveBeamSPs(PDSPSaveBeamSPs&&) = delete;
  PDSPSaveBeamSPs& operator=(PDSPSaveBeamSPs const&) = delete;
  PDSPSaveBeamSPs& operator=(PDSPSaveBeamSPs&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void endJob() override;
  void beginJob() override;
  void reset();
private:

  //Persistent vectors for writeout 
  std::vector<std::vector<double>> fSpacePointX, fSpacePointY, fSpacePointZ,
                                   fSpacePointCharge0, fSpacePointCharge1,
                                   fSpacePointCharge2;

  TTree *fTree;
  std::vector<double> fSpacePointX_OneDim, fSpacePointY_OneDim, fSpacePointZ_OneDim,
                      fSpacePointCharge0_OneDim, fSpacePointCharge1_OneDim,
                      fSpacePointCharge2_OneDim,
                      fSpacePointRMS0_OneDim, fSpacePointRMS1_OneDim,
                      fSpacePointRMS2_OneDim;

  int fRecoBeamOrigin;
  double fRecoHierarchyBeamFraction;
  int fTrueBeamID, fTrueBeamPDG;
  std::string fTrueBeamEndProcess;
  double fTrueBeamEndX, fTrueBeamEndY, fTrueBeamEndZ;
  int fTrueBeam_nPiPlus, fTrueBeam_nPiMinus, fTrueBeam_nPi0, fTrueBeam_nProton,
      fTrueBeam_nNeutron, fTrueBeam_nNucleus, fTrueBeam_nKaon;

  size_t fNEvents = 0;
  std::map<art::Ptr<recob::SpacePoint>, std::vector<const recob::Hit*>>
      fSpacePointHitMap;

  std::vector<SPEdge> fEdges;

  std::vector<std::vector<int>> fEdge_i, fEdge_j;
  std::vector<std::vector<double>>
      fEdge_dX, fEdge_dY, fEdge_dZ,
      fEdge_dC0, fEdge_dC1, fEdge_dC2,
      fEdge_AveC0, fEdge_AveC1, fEdge_AveC2,
      fEdge_OP0, fEdge_OP1, fEdge_OP2,
      fEdge_OP3, fEdge_OP4, fEdge_OP5,
      fEdge_COP0, fEdge_COP1, fEdge_COP2,
      fEdge_COP3, fEdge_COP4, fEdge_COP5;

  //Configuration 
  std::string fPFParticleTag, fSpacePointTag, fHitTag, fGeneratorTag;
  bool fDebug = false;


  //Utils
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;


  //Methods
  void GetTruthInfo(const simb::MCParticle * particle);
  void GetOriginInfo(//TODO -- Better name
    const art::Event & evt,
    const detinfo::DetectorClocksData & clockData,
    const recob::PFParticle * particle,
    std::pair<double, double> & origin_energies);
  void ResizeSPVectors();
  void FillSPVectors();
  void MakeEdges();
  void GetSpacePoints(
      const art::Event & evt,
      const art::FindManyP<recob::SpacePoint> & space_pts_from_hits,
      const recob::PFParticle * particle);
};


pduneana::PDSPSaveBeamSPs::PDSPSaveBeamSPs(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fPFParticleTag(p.get<std::string>("PFParticleTag", "pandora")),
    fSpacePointTag(p.get<std::string>("SpacePointTag", "reco3d")),
    fHitTag(p.get<std::string>("HitTag", "hitpdune")),
    fGeneratorTag(p.get<std::string>("GeneratorTag", "generator")),
    fDebug(p.get<bool>("Debug", false))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pduneana::PDSPSaveBeamSPs::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("space_points", "");

  fTree->Branch("x", &fSpacePointX_OneDim);
  fTree->Branch("y", &fSpacePointY_OneDim);
  fTree->Branch("z", &fSpacePointZ_OneDim);

  fTree->Branch("charge0", &fSpacePointCharge0_OneDim);
  fTree->Branch("charge1", &fSpacePointCharge1_OneDim);
  fTree->Branch("charge2", &fSpacePointCharge2_OneDim);
  fTree->Branch("rms0", &fSpacePointRMS0_OneDim);
  fTree->Branch("rms1", &fSpacePointRMS1_OneDim);
  fTree->Branch("rms2", &fSpacePointRMS2_OneDim);

  fTree->Branch("edge_i", &fEdge_i);
  fTree->Branch("edge_j", &fEdge_j);
  fTree->Branch("edge_dx", &fEdge_dX);
  fTree->Branch("edge_dy", &fEdge_dY);
  fTree->Branch("edge_dz", &fEdge_dZ);
  fTree->Branch("edge_dc0", &fEdge_dC0);
  fTree->Branch("edge_dc1", &fEdge_dC1);
  fTree->Branch("edge_dc2", &fEdge_dC2);
  fTree->Branch("edge_avec0", &fEdge_AveC0);
  fTree->Branch("edge_avec1", &fEdge_AveC1);
  fTree->Branch("edge_avec2", &fEdge_AveC2);
  fTree->Branch("edge_op0", &fEdge_OP0);
  fTree->Branch("edge_op1", &fEdge_OP1);
  fTree->Branch("edge_op2", &fEdge_OP2);
  fTree->Branch("edge_op3", &fEdge_OP0);
  fTree->Branch("edge_op4", &fEdge_OP4);
  fTree->Branch("edge_op5", &fEdge_OP5);
  fTree->Branch("edge_cop0", &fEdge_COP0);
  fTree->Branch("edge_cop1", &fEdge_COP1);
  fTree->Branch("edge_cop2", &fEdge_COP2);
  fTree->Branch("edge_cop3", &fEdge_COP0);
  fTree->Branch("edge_cop4", &fEdge_COP4);
  fTree->Branch("edge_cop5", &fEdge_COP5);

  fTree->Branch("reco_origin", &fRecoBeamOrigin);
  fTree->Branch("reco_beam_fraction", &fRecoHierarchyBeamFraction);
  fTree->Branch("true_id", &fTrueBeamID);
  fTree->Branch("true_pdg", &fTrueBeamPDG);
  fTree->Branch("true_end_process", &fTrueBeamEndProcess);
  fTree->Branch("true_end_x", &fTrueBeamEndX);
  fTree->Branch("true_end_y", &fTrueBeamEndY);
  fTree->Branch("true_end_z", &fTrueBeamEndZ);

  fTree->Branch("true_n_piplus", &fTrueBeam_nPiPlus);
  fTree->Branch("true_n_piminus", &fTrueBeam_nPiMinus);
  fTree->Branch("true_n_pi0", &fTrueBeam_nPi0);
  fTree->Branch("true_n_proton", &fTrueBeam_nProton);
  fTree->Branch("true_n_neutron", &fTrueBeam_nNeutron);
  fTree->Branch("true_n_nucleus", &fTrueBeam_nNucleus);
  fTree->Branch("true_n_kaon", &fTrueBeam_nKaon);
}

void pduneana::PDSPSaveBeamSPs::reset() {
  fSpacePointX_OneDim.clear();
  fSpacePointY_OneDim.clear();
  fSpacePointZ_OneDim.clear();
  fSpacePointCharge0_OneDim.clear();
  fSpacePointCharge1_OneDim.clear();
  fSpacePointCharge2_OneDim.clear(); 
  fSpacePointRMS0_OneDim.clear();
  fSpacePointRMS1_OneDim.clear();
  fSpacePointRMS2_OneDim.clear(); 

  fRecoBeamOrigin = -999;
  fRecoHierarchyBeamFraction = -999.;
  fTrueBeamID = -999;
  fTrueBeamPDG = -999;
  fTrueBeamEndProcess = "";
  fTrueBeamEndX = -999.;
  fTrueBeamEndY = -999.;
  fTrueBeamEndZ = -999.;

  fTrueBeam_nPiPlus = 0;
  fTrueBeam_nPiMinus = 0;
  fTrueBeam_nPi0 = 0;
  fTrueBeam_nProton = 0;
  fTrueBeam_nNeutron = 0;
  fTrueBeam_nNucleus = 0;
  fTrueBeam_nKaon = 0;

  fEdges.clear();
  fEdge_i.clear();
  fEdge_j.clear();
  fEdge_dX.clear();
  fEdge_dY.clear();
  fEdge_dZ.clear();
  fEdge_dC0.clear();
  fEdge_dC1.clear();
  fEdge_dC2.clear();
  fEdge_AveC0.clear();
  fEdge_AveC1.clear();
  fEdge_AveC2.clear();
  fEdge_OP0.clear();
  fEdge_OP1.clear();
  fEdge_OP2.clear();
  fEdge_OP3.clear();
  fEdge_OP4.clear();
  fEdge_OP5.clear();
  fEdge_COP0.clear();
  fEdge_COP1.clear();
  fEdge_COP2.clear();
  fEdge_COP3.clear();
  fEdge_COP4.clear();
  fEdge_COP5.clear();
}

void pduneana::PDSPSaveBeamSPs::endJob() {
  /*
  std::cout << "Event summaries:" << std::endl;
  for (size_t i = 0; i < fNEvents; ++i) {
    std::cout << "Event " << i << std::endl;
    std::cout << "\tSizes " <<
                 fSpacePointX[i].size() << " " <<
                 fSpacePointY[i].size() << " " <<
                 fSpacePointZ[i].size() << " " <<
                 fSpacePointCharge0[i].size() << " " <<
                 fSpacePointCharge1[i].size() << " " <<
                 fSpacePointCharge2[i].size() << " " <<
                 std::endl; 
  }

  fFile = hep_hpc::hdf5::File("test.hdf5", H5F_ACC_TRUNC);
  */
}

void pduneana::PDSPSaveBeamSPs::GetSpacePoints(
    const art::Event & evt,
    const art::FindManyP<recob::SpacePoint> & space_pts_from_hits,
    const recob::PFParticle * particle) {

  const auto particle_hits
      = pfpUtil.GetPFParticleHits_Ptrs(*particle, evt, fPFParticleTag);
  if (fDebug) 
    std::cout << "Found N hits from particle: " << particle_hits.size() <<
                 std::endl;
  for (const auto & hit : particle_hits) {
    const auto & space_points = space_pts_from_hits.at(hit.key());
    for (const auto & space_point : space_points) {
      fSpacePointHitMap[space_point].push_back(&(*hit)); //Convert to normal ptr?
    }
  }
}

void pduneana::PDSPSaveBeamSPs::MakeEdges() {
  for (size_t i = 0; i < fSpacePointX_OneDim.size(); ++i) {
    SPNode node_i(fSpacePointX_OneDim[i],
                  fSpacePointY_OneDim[i],
                  fSpacePointZ_OneDim[i],
                  fSpacePointCharge0_OneDim[i],
                  fSpacePointCharge1_OneDim[i],
                  fSpacePointCharge2_OneDim[i],
                  fSpacePointRMS0_OneDim[i],
                  fSpacePointRMS1_OneDim[i],
                  fSpacePointRMS2_OneDim[i]
                 );
    for (size_t j = 0; j < fSpacePointX_OneDim.size(); ++j) {
      if (i == j) continue;
      SPNode node_j(fSpacePointX_OneDim[j],
                    fSpacePointY_OneDim[j],
                    fSpacePointZ_OneDim[j],
                    fSpacePointCharge0_OneDim[j],
                    fSpacePointCharge1_OneDim[j],
                    fSpacePointCharge2_OneDim[j],
                    fSpacePointRMS0_OneDim[j],
                    fSpacePointRMS1_OneDim[j],
                    fSpacePointRMS2_OneDim[j]
                   );
      fEdges.push_back(SPEdge(i, j, node_i, node_j));
    }
    std::cout << "N Edges: " << fEdges.size() << std::endl;
    std::sort(fEdges.begin(), fEdges.end(),
              [](SPEdge & a, SPEdge & b) { return a.fMag < b.fMag; });
    std::cout << "Sorted" << std::endl;
    if (fEdges.size() > 1) {
      std::cout << fEdges[0].fMag << " " << fEdges[1].fMag << " " <<
                   fEdges.back().fMag << std::endl;
    }

    fEdge_i.push_back({});
    fEdge_j.push_back({});
    fEdge_dX.push_back({});
    fEdge_dY.push_back({});
    fEdge_dZ.push_back({});

    fEdge_dC0.push_back({});
    fEdge_dC1.push_back({});
    fEdge_dC2.push_back({});

    fEdge_AveC0.push_back({});
    fEdge_AveC1.push_back({});
    fEdge_AveC2.push_back({});

    fEdge_OP0.push_back({});
    fEdge_OP1.push_back({});
    fEdge_OP2.push_back({});
    fEdge_OP3.push_back({});
    fEdge_OP4.push_back({});
    fEdge_OP5.push_back({});

    fEdge_COP0.push_back({});
    fEdge_COP1.push_back({});
    fEdge_COP2.push_back({});
    fEdge_COP3.push_back({});
    fEdge_COP4.push_back({});
    fEdge_COP5.push_back({});

    //Max 50 edges for right now
    size_t max = ((50 <= fEdges.size()) ? 50 : fEdges.size());
    for (size_t j = 0; j < max; ++j) {
      const auto & edge = fEdges[j];
      //fEdge_i.back().push_back(i);
      //fEdge_j.back().push_back(j);
      fEdge_i.back().push_back(edge.fI);
      fEdge_j.back().push_back(edge.fJ);
      fEdge_dX.back().push_back(edge.fDiffX);
      fEdge_dY.back().push_back(edge.fDiffY);
      fEdge_dZ.back().push_back(edge.fDiffZ);

      fEdge_dC0.back().push_back(edge.fDiffC0);
      fEdge_dC1.back().push_back(edge.fDiffC1);
      fEdge_dC2.back().push_back(edge.fDiffC2);

      fEdge_AveC0.back().push_back(edge.fAveCharge[0]);
      fEdge_AveC1.back().push_back(edge.fAveCharge[1]);
      fEdge_AveC2.back().push_back(edge.fAveCharge[2]);

      fEdge_OP0.back().push_back(edge.fDiffOP[0]);
      fEdge_OP1.back().push_back(edge.fDiffOP[1]);
      fEdge_OP2.back().push_back(edge.fDiffOP[2]);
      fEdge_OP3.back().push_back(edge.fDiffOP[3]);
      fEdge_OP4.back().push_back(edge.fDiffOP[4]);
      fEdge_OP5.back().push_back(edge.fDiffOP[5]);

      fEdge_COP0.back().push_back(edge.fDiffCOP[0]);
      fEdge_COP1.back().push_back(edge.fDiffCOP[1]);
      fEdge_COP2.back().push_back(edge.fDiffCOP[2]);
      fEdge_COP3.back().push_back(edge.fDiffCOP[3]);
      fEdge_COP4.back().push_back(edge.fDiffCOP[4]);
      fEdge_COP5.back().push_back(edge.fDiffCOP[5]);
    }

    fEdges.clear();
  }
}

void pduneana::PDSPSaveBeamSPs::FillSPVectors() {
  //Need to call this first to append an empty vec
  ResizeSPVectors();

  for (const auto & sp_hits : fSpacePointHitMap) {
    const auto space_point = sp_hits.first;
    fSpacePointX.back().push_back(space_point->XYZ()[0]);
    fSpacePointY.back().push_back(space_point->XYZ()[1]);
    fSpacePointZ.back().push_back(space_point->XYZ()[2]);

    fSpacePointX_OneDim.push_back(space_point->XYZ()[0]);
    fSpacePointY_OneDim.push_back(space_point->XYZ()[1]);
    fSpacePointZ_OneDim.push_back(space_point->XYZ()[2]);

    bool found_0 = false, found_1 = false, found_2 = false;
    for (const auto * hit : sp_hits.second) {

      if (hit->WireID().Plane == 0) {
        found_0 = true;
        fSpacePointCharge0_OneDim.push_back(hit->Integral());
        fSpacePointRMS0_OneDim.push_back(hit->RMS());
      }
      else if (hit->WireID().Plane == 1) {
        found_1 = true;
        fSpacePointCharge1_OneDim.push_back(hit->Integral());
        fSpacePointRMS1_OneDim.push_back(hit->RMS());
      }
      else if (hit->WireID().Plane == 2) {
        found_2 = true;
        fSpacePointCharge2_OneDim.push_back(hit->Integral());
        fSpacePointRMS2_OneDim.push_back(hit->RMS());
      }
      else {
        throw cet::exception("PDSPSaveBeamSPs_module.cc") <<
              "Error. Found additional plane? " << hit->WireID().Plane;
      }
    }

    if (!found_0) {
      fSpacePointCharge0_OneDim.push_back(0.);
      fSpacePointRMS0_OneDim.push_back(0.);
    }
    if (!found_1) {
      fSpacePointCharge1_OneDim.push_back(0.);
      fSpacePointRMS1_OneDim.push_back(0.);
    }
    if (!found_2) {
      fSpacePointCharge2_OneDim.push_back(0.);
      fSpacePointRMS2_OneDim.push_back(0.);
    }
  }
}

void pduneana::PDSPSaveBeamSPs::ResizeSPVectors() {
  fSpacePointX.push_back({});
  fSpacePointY.push_back({});
  fSpacePointZ.push_back({});
  fSpacePointCharge0.push_back({});
  fSpacePointCharge1.push_back({});
  fSpacePointCharge2.push_back({});
}

void pduneana::PDSPSaveBeamSPs::GetTruthInfo(
    const simb::MCParticle * particle) {
  //Beam Info
  fTrueBeamID = particle->TrackId();
  fTrueBeamPDG = particle->PdgCode();
  fTrueBeamEndProcess = particle->EndProcess();
  fTrueBeamEndX = particle->EndX();
  fTrueBeamEndY = particle->EndY();
  fTrueBeamEndZ = particle->EndZ();

  //Product Info
  const sim::ParticleList & plist = pi_serv->ParticleList();
  for (int i = 0; i < particle->NumberDaughters(); ++i) {
    int id = particle->Daughter(i);
    const auto * part = plist[id];
    if (fTrueBeamEndProcess.find("Inelastic") != std::string::npos) {
      int pdg = part->PdgCode();
      if (pdg == 211) ++fTrueBeam_nPiPlus;
      else if (pdg == -211) ++fTrueBeam_nPiMinus;
      else if (pdg == 111) ++fTrueBeam_nPi0;
      else if (pdg == 2212) ++fTrueBeam_nProton;
      else if (pdg == 2112) ++fTrueBeam_nNeutron;
      else if (pdg > 2212) ++fTrueBeam_nNucleus;
      else if (abs(pdg) == 321) ++fTrueBeam_nKaon;
    }
  }
}

void pduneana::PDSPSaveBeamSPs::GetOriginInfo(//TODO -- Better name
    const art::Event & evt,
    const detinfo::DetectorClocksData & clockData,
    const recob::PFParticle * particle,
    std::pair<double, double> & origin_energies) {

  // Loop over the hits in the particle 
  const auto particle_hits
      = pfpUtil.GetPFParticleHits_Ptrs(*particle, evt, fPFParticleTag);
  for (const auto & hit : particle_hits) {
    const auto & ides = bt_serv->HitToEveTrackIDEs(clockData, *hit);
    for (const auto & ide : ides) {
      int origin = pi_serv->TrackIdToMCTruth_P(ide.trackID)->Origin();
      if (origin == 4) {
        origin_energies.first += ide.energy;
      }
      else if (origin == 2){
        origin_energies.second += ide.energy;
      }
      else {
        std::cout << "Warning: found other origin? " << origin << std::endl;
      }
    }
  }
}

void pduneana::PDSPSaveBeamSPs::analyze(art::Event const& evt) {
  reset();
  auto const clockData
      = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);


  //Get Beam particle if it is there
  std::vector<const recob::PFParticle*> beamParticles
      = pfpUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleTag);
  if (beamParticles.size() == 0) return;
  if (fDebug) std::cout << "N Beam Particles: " << beamParticles.size() << std::endl;
  const auto * particle = beamParticles.at(0);
  if (fDebug) std::cout << "Found particle: " << particle << std::endl;

  //Get the hits and associated space points
  const auto all_hits = evt.getValidHandle<std::vector<recob::Hit>>(fHitTag);
  if (fDebug) {
    std::cout << "Is valid? " << all_hits.isValid() << std::endl;
    std::cout << "N hits in event: " << all_hits->size() << std::endl;
  }
  art::FindManyP<recob::SpacePoint> space_pts_from_hits(
      all_hits, evt, fHitTag);

  if (fDebug) 
    std::cout << "Found " << space_pts_from_hits.size() <<
                 " space points from hits" << std::endl;


  fSpacePointHitMap.clear(); //Clear space point map to free memory
  GetSpacePoints(evt, space_pts_from_hits, particle);
  std::pair<double, double> origin_energies = {0., 0.};
  GetOriginInfo(evt, clockData, particle, origin_energies);
  if (fDebug)
    std::cout << "Energies: " << origin_energies.first << " " <<
                 origin_energies.second << std::endl;

  std::deque<size_t> hierarchy;
  auto pfpVec
      = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  for (size_t daughterID : particle->Daughters()) {
    if (fDebug) std::cout << "Added " << daughterID << std::endl;
    hierarchy.push_back(daughterID);
  }

  while (hierarchy.size() > 0) {
    size_t id = hierarchy[0];
    if (fDebug) std::cout << "Getting SPs for " << id << std::endl;
    const recob::PFParticle * pfp = &(pfpVec->at(id));
    GetSpacePoints(evt, space_pts_from_hits, pfp);
    for (auto d : pfp->Daughters()) {
      if (fDebug) std::cout << "\tAdded " << d << std::endl;
      hierarchy.push_back(d);
    }

    GetOriginInfo(evt, clockData, pfp, origin_energies);
    if (fDebug)
      std::cout << "Energies: " << origin_energies.first << " " <<
                   origin_energies.second << std::endl;

    hierarchy.pop_front();
  }


  ++fNEvents;

  FillSPVectors();
  MakeEdges();

  //Truth info
  if (!evt.isRealData()) {
    auto mcTruths
        = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    auto * true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0], evt);
    if (true_beam_particle == 0x0) //Make sure there is some truth info
      return;

    GetTruthInfo(true_beam_particle);

    auto beam_match
        = truthUtil.GetMCParticleByHits(clockData, *particle, evt,
                                        fPFParticleTag, fHitTag);
    fRecoBeamOrigin
        = pi_serv->TrackIdToMCTruth_P(beam_match.particle->TrackId())->Origin();
    
  }

  double total_energy = origin_energies.first + origin_energies.second;
  fRecoHierarchyBeamFraction = origin_energies.first / total_energy;
  fTree->Fill();
}

DEFINE_ART_MODULE(pduneana::PDSPSaveBeamSPs)
