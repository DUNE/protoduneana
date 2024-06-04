////////////////////////////////////////////////////////////////////////
// Class:       PDSPAnalyzer
// Plugin Type: analyzer (art v3_00_00)
// File:        PDSPAnalyzer_module.cc
// Written by Jake Calcutt (calcuttj@msu.edu -- Slack: @jakecalcutt)
// Reach out for questions/issues/bugs
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

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
#include "protoduneana/Utilities/ProtoDUNEEmptyEventFinder.h"

#include "protoduneana/Utilities/G4ReweightUtils.h"
//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNEDataUtils.h"

//#include "duneprototypes/Protodune/singlephase/DataUtils/ProtoDUNECalibration.h"
#include "protoduneana/Utilities/ProtoDUNECalibration.h"
#include "protoduneana/Utilities/AbsCexReweighter.hh"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/RawData/RDTimeStamp.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"

#include "lardata/ArtDataHelper/MVAReader.h"


#include "geant4reweight/src/ReweightBase/G4ReweighterFactory.hh"
#include "geant4reweight/src/ReweightBase/G4Reweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"
#include "geant4reweight/src/PropBase/G4ReweightParameterMaker.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"
#include "geant4reweight/src/ReweightBase/G4ReweightManager.hh"

#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"

#include "art_root_io/TFileService.h"
#include "TProfile.h"
#include "TFile.h"

// ROOT includes
#include "TTree.h"
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TPolyLine3D.h"
#include "Math/Vector3D.h"
#include "TGraph2D.h"
#include "TROOT.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"


namespace pduneana {

  using std::string;
  using std::vector;
  using std::pair;
  using std::map;
  using std::set;

  using art::Ptr;
  using art::Event;
  using art::ServiceHandle;

  using cheat::BackTrackerService;
  using cheat::ParticleInventoryService;

  using simb::MCParticle;

  using recob::Hit;
  using recob::SpacePoint;


  // Get the angle and the dot product between the vector from the base node to its neighbours
  void GetAngleAndDotProduct(const SpacePoint &baseNode,
    const SpacePoint &n1, const SpacePoint &n2, double &dotProduct, double &angle) {
    TVector3 basePos(baseNode.XYZ());
    TVector3 neighbour1Pos(n1.XYZ());
    TVector3 neighbour2Pos(n2.XYZ());
    TVector3 baseToNeighbour1 = neighbour1Pos - basePos;
    TVector3 baseToNeighbour2 = neighbour2Pos - basePos;
    dotProduct = baseToNeighbour1.Dot(baseToNeighbour2);
    angle = baseToNeighbour1.Angle(baseToNeighbour2);
    return;
  }




  // Use the association between space points and hits to return a charge
  std::map<unsigned int, float> GetSpacePointChargeMap(
    std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit) {

    map<unsigned int, float> ret;

    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      float charge = 0.0;
      for (Ptr<Hit> hit : sp2Hit[spIdx]) {
        charge += hit->Integral();
      }
      ret[spacePoints[spIdx]->ID()] = charge;
    }

    return ret;

  } // function GetSpacePointChargeMap

  std::map<unsigned int, float> GetSpacePointChargeMap(
    art::Event const &evt, const std::string &spLabel)  {

    art::Handle<vector<SpacePoint>> spacePointHandle;
    vector<Ptr<SpacePoint>> spacePoints;
    if (!evt.getByLabel(spLabel, spacePointHandle)) {
      throw art::Exception(art::errors::LogicError)
        << "Could not find spacepoints with module label "
        << spLabel << "!";
    }
    art::fill_ptr_vector(spacePoints, spacePointHandle);
    art::FindManyP<Hit> fmp(spacePointHandle, evt, spLabel);
    vector<vector<Ptr<Hit>>> sp2Hit(spacePoints.size());
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    } // for spacepoint

    return GetSpacePointChargeMap(spacePoints, sp2Hit);

  } // function GetSpacePointChargeMap




  // Get the two nearest neighbours to use for calcuation of angles between them and the node in question
  

  std::map<int,std::pair<int,int>> GetTwoNearestNeighbours(art::Event const &evt,
    const std::vector<const recob::SpacePoint*> &sps) {
    std::map<int,int> closestID;
    std::map<int,int> secondID;
    // Now loop over all of the space points and simultaneously fill our maps
    // Get the space points from the event and make the map
    for(auto sp0 : sps){
      // We want an entry even if it ends up being zero
      int thisSP = sp0->ID();
      int closest = -1;
      int second = -1;
      float closestDist = 99999;
      float secondDist = 99999;
      for(auto sp1 : sps){
        if(thisSP == sp1->ID()) continue;
        // For some reason we have to use arrays
        const double *p0 = sp0->XYZ();
        const double *p1 = sp1->XYZ();
        // Get the distance between the points
        const float dx = p1[0] - p0[0];
        const float dy = p1[1] - p0[1];
        const float dz = p1[2] - p0[2];
        const float dist = sqrt(dx*dx + dy*dy + dz*dz);
        if(dist < closestDist){
          secondDist = closestDist;
          closestDist = dist;
          second = closest;
          closest = sp1->ID();
        }
        else if(dist < secondDist){
          secondDist = dist;
          second = sp1->ID();
        }
      }
      closestID.insert(std::make_pair(thisSP,closest));
      secondID.insert(std::make_pair(thisSP,second));
    }
    map<int,pair<int,int>> finalMap;
    for(unsigned int m = 0; m < closestID.size(); ++m){
      finalMap[m] = std::make_pair(closestID[m],secondID[m]);
    }
    return finalMap;
  }

/*
  std::map<int,std::pair<int,int>> GetTwoNearestNeighbours_label(art::Event const &evt,
    const std::string &spLabel) {
    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
    if(evt.getByLabel(spLabel,spacePointHandle)){
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }
    return GetTwoNearestNeighbours(evt,allSpacePoints);
  }

  
  std::map<int,std::pair<int,int>> GetTwoNearestNeighbours_vec(art::Event const &evt,
    const std::map<unsigned int, art::Ptr<recob::SpacePoint>> &sps)  {

    vector<Ptr<SpacePoint>> vec;
    for(auto m : sps){
      vec.push_back(m.second);
    }
    return GetTwoNearestNeighbours(evt,vec);
  }

*/




  std::vector<std::map<int,unsigned int>> GetNeighboursForRadii(art::Event const &evt,
    const std::vector<float>& rangeCuts, const std::vector<float>& chargeRangeCuts, const std::vector<const recob::SpacePoint*> &sps) {
    std::vector<std::map<int,unsigned int>> result;
    // Initialise a map for each range cut value
    for(unsigned int m = 0; m < rangeCuts.size(); ++m){
      result.push_back(map<int,unsigned int>());
    }
    // Initialize a map for each charge range cut value
    for(unsigned int m = 0; m < chargeRangeCuts.size(); ++m){
      result.push_back(map<int,unsigned int>());
    }
 
    // Comput the charge map
    // NOTE: space point label is hardcoded
    std::map<unsigned int,float> chargeMap = GetSpacePointChargeMap(evt, "pandora");

 
    for(auto sp0 : sps){
      // We want an entry even if it ends up being zero
      for(auto &m : result){
        m[sp0->ID()] = 0;
      }
      for(auto sp1 : sps){
        if(sp0->ID() == sp1->ID()) continue;
        // For some reason we have to use arrays
        const double *p0 = sp0->XYZ();
        const double *p1 = sp1->XYZ();
        // Get the distance between the points
        const float dx = p1[0] - p0[0];
        const float dy = p1[1] - p0[1];
        const float dz = p1[2] - p0[2];
        const float dist = sqrt(dx*dx + dy*dy + dz*dz);
        // Fill the maps if we satify the criteria
        for(unsigned int r = 0; r < rangeCuts.size(); ++r){
          if(dist < rangeCuts[r]){
            result[r][sp0->ID()] = result[r][sp0->ID()] + 1;
          }
        }
        // Fill the maps with the charge within that radius
        for(unsigned int r = 0; r < chargeRangeCuts.size(); ++r){
          if(dist < chargeRangeCuts[r]){
            float charge = chargeMap.at(sp1->ID()); 
            // Instead of adding 1 I should add the charge of SP1!! (to be confirmed by Leigh)
            result[r+rangeCuts.size()][sp0->ID()] = result[r+rangeCuts.size()][sp0->ID()] + charge;
          }
        }
 


      }
    }
    return result;
  }

/*
  // Sometimes we might want to know the number of neighbours within various radii
  std::vector<std::map<int,unsigned int>> GetNeighboursForRadii_label(art::Event const &evt, const std::vector<float>& rangeCuts, const std::vector<float>& chargeRangeCuts, const std::string &spLabel) {
    // Get the space points from the event and make the map
    art::Handle<vector<SpacePoint>> spacePointHandle;
    vector<Ptr<SpacePoint>> allSpacePoints;
    if(evt.getByLabel(spLabel,spacePointHandle)){
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }
    return GetNeighboursForRadii(evt,rangeCuts, chargeRangeCuts, allSpacePoints);
  }

  std::vector<std::map<int,unsigned int>> GetNeighboursForRadii_vec(art::Event const &evt,
    const std::vector<float>& rangeCuts, const std::vector<float>& chargeRangeCuts, const std::map<unsigned int,art::Ptr<recob::SpacePoint>> &sps) {
    std::vector<art::Ptr<recob::SpacePoint>> vec;
    for(auto m : sps){
      vec.push_back(m.second);
    }
    return GetNeighboursForRadii(evt, rangeCuts, chargeRangeCuts, vec);
  }

*/





  class PDSPAnalyzer;


  //Used to fit lines to sets of reconstructed hits in order to get
  //directions
  double distance2(double x, double y, double z, double * p);
  void line(double t, double * p, double & x, double & y, double & z);
  void SumDistance2(int &, double *, double & sum, double * par, int);
  TVector3 FitLine(const std::vector<TVector3> & input);

  //Used to create truth thin slices, though this is sort of 
  //deprecated because IDEs get ate up
  bool sort_IDEs( const sim::IDE * i1, const sim::IDE * i2){
    return( i1->z < i2->z );
  }
  std::map<int, std::vector<const sim::IDE*>> slice_IDEs(
      std::vector<const sim::IDE*> ides, double the_z0, double the_pitch,
      double true_endZ){

    std::map<int, std::vector<const sim::IDE*>> results;

    for (size_t i = 0; i < ides.size(); ++i) {
      int slice_num = std::floor(
          (ides[i]->z - (the_z0 - the_pitch/2.)) / the_pitch);
      /*
      std::cout << "IDE: " << i << " ID: " << ides[i]->trackID << " Edep: "
                << ides[i]->energy << " (X,Y,Z): " << "(" << ides[i]->x << ","
                << ides[i]->y<<","<<ides[i]->z << ") Z0: " << the_z0
                << " Slice: " << slice_num << std::endl;*/
      results[slice_num].push_back(ides[i]);
    }

    return results;
  }

  
  //To make auto-finding files easier for the dEdX template file
  TFile * OpenFile(const std::string filename) {
    TFile * theFile = 0x0;
    mf::LogInfo("pduneana::OpenFile") << "Searching for " << filename;
    if (cet::file_exists(filename)) {
      mf::LogInfo("pduneana::OpenFile") << "File exists. Opening " << filename;
      theFile = new TFile(filename.c_str());
      if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
        delete theFile;
        theFile = 0x0;
        throw cet::exception("PDSPAnalyzer_module.cc") << "Could not open " << filename;
      }
    }
    else {
      mf::LogInfo("pduneana::OpenFile") << "File does not exist here. Searching FW_SEARCH_PATH";
      cet::search_path sp{"FW_SEARCH_PATH"};
      std::string found_filename;
      auto found = sp.find_file(filename, found_filename);
      if (!found) {
        throw cet::exception("PDSPAnalyzer_module.cc") << "Could not find " << filename;
      }

      mf::LogInfo("pduneana::OpenFile") << "Found file " << found_filename;
      theFile = new TFile(found_filename.c_str());
      if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) {
        delete theFile;
        theFile = 0x0;
        throw cet::exception("PDSPAnalyzer_module.cc") << "Could not open " << found_filename;
      }
    }
    return theFile;
  }

  //Struct to store CNN output
  struct cnnOutput2D{
    cnnOutput2D();

    double track;
    double em;
    double michel;
    double none;
    size_t nHits;
    double track_weight_by_charge;
    double em_weight_by_charge;
    double michel_weight_by_charge;
    double none_weight_by_charge;
  };

  //Struct to handle calorimetry info more easily
  struct calo_point{

    calo_point();
    calo_point(size_t w, double in_tick, double p, double dqdx, double dedx, double dq,
               double cali_dqdx, double cali_dedx, double r, size_t index, double input_wire_z, int t,
               double efield, double input_x, double input_y, double input_z)
        : wire(w), tick(in_tick), pitch(p), dQdX(dqdx), dEdX(dedx), dQ(dq),
          calibrated_dQdX(cali_dqdx), calibrated_dEdX(cali_dedx),
          res_range(r), hit_index(index), wire_z(input_wire_z), tpc(t),
          EField(efield), x(input_x), y(input_y), z(input_z) {};

    size_t wire;
    double tick;
    double pitch;
    double dQdX;
    double dEdX;
    double dQ;
    double calibrated_dQdX;
    double calibrated_dEdX;
    double res_range;
    size_t hit_index;
    double wire_z;
    int tpc;
    double EField;
    double x, y, z;

    std::vector<double> IDE_electrons, IDE_energies;
    std::vector<int> IDE_IDs;
    std::vector<int> IDE_origins;

    void SetIDE_electrons(std::vector<double> input) {IDE_electrons = input;};
    void SetIDE_energies(std::vector<double> input) {IDE_energies = input;};
    void SetIDE_IDs(std::vector<int> input) {IDE_IDs = input;};
    void SetIDE_origins(std::vector<int> input) {IDE_origins = input;};
  };

  //Util to get CNN output
  cnnOutput2D GetCNNOutputFromPFParticle(
      const recob::PFParticle & part, const art::Event & evt,
      const anab::MVAReader<recob::Hit,4> & CNN_results,
      protoana::ProtoDUNEPFParticleUtils & pfpUtil,
      std::string fPFParticleTag) {

    cnnOutput2D output;
    const std::vector< art::Ptr< recob::Hit > > daughterPFP_hits
        = pfpUtil.GetPFParticleHits_Ptrs(part, evt, fPFParticleTag);

    double tot_charge = 0.;
    for (size_t h = 0; h < daughterPFP_hits.size(); ++h) {
      std::array<float,4> cnn_out = CNN_results.getOutput( daughterPFP_hits[h] );
      double hitcharge = daughterPFP_hits[h]->Integral();
      double track_score = cnn_out[ CNN_results.getIndex("track") ];
      double em_score = cnn_out[ CNN_results.getIndex("em") ];
      double michel_score = cnn_out[ CNN_results.getIndex("michel") ];
      double none_score = cnn_out[ CNN_results.getIndex("none") ];
      output.track  += track_score;
      output.em     += em_score;
      output.michel += michel_score;
      output.none   += none_score;
      output.track_weight_by_charge  += hitcharge*track_score;
      output.em_weight_by_charge     += hitcharge*em_score;
      output.michel_weight_by_charge += hitcharge*michel_score;
      output.none_weight_by_charge   += hitcharge*none_score;
      tot_charge += hitcharge;
    }
    output.nHits = daughterPFP_hits.size();
    if (tot_charge != 0) {
      output.track_weight_by_charge  /= tot_charge;
      output.em_weight_by_charge     /= tot_charge;
      output.michel_weight_by_charge /= tot_charge;
      output.none_weight_by_charge   /= tot_charge;
    }
    else {
      output.track_weight_by_charge  = -999.;
      output.em_weight_by_charge     = -999.;
      output.michel_weight_by_charge = -999.;
      output.none_weight_by_charge   = -999.;
    }

    return output;
  }


  //Another CNN output util
  cnnOutput2D GetCNNOutputFromPFParticleFromPlane( const recob::PFParticle & part, const art::Event & evt, const anab::MVAReader<recob::Hit,4> & CNN_results,  protoana::ProtoDUNEPFParticleUtils & pfpUtil, std::string fPFParticleTag, size_t planeID ){

    cnnOutput2D output;
    const std::vector< art::Ptr< recob::Hit > > daughterPFP_hits = pfpUtil.GetPFParticleHitsFromPlane_Ptrs( part, evt, fPFParticleTag, planeID );

    double tot_charge = 0.;
    for (size_t h = 0; h < daughterPFP_hits.size(); ++h) {
      std::array<float,4> cnn_out = CNN_results.getOutput( daughterPFP_hits[h] );
      double hitcharge = daughterPFP_hits[h]->Integral();
      double track_score = cnn_out[ CNN_results.getIndex("track") ];
      double em_score = cnn_out[ CNN_results.getIndex("em") ];
      double michel_score = cnn_out[ CNN_results.getIndex("michel") ];
      double none_score = cnn_out[ CNN_results.getIndex("none") ];
      output.track  += track_score;
      output.em     += em_score;
      output.michel += michel_score;
      output.none   += none_score;
      output.track_weight_by_charge  += hitcharge*track_score;
      output.em_weight_by_charge     += hitcharge*em_score;
      output.michel_weight_by_charge += hitcharge*michel_score;
      output.none_weight_by_charge   += hitcharge*none_score;
      tot_charge += hitcharge;
    }
    output.nHits = daughterPFP_hits.size();
    if (tot_charge != 0) {
      output.track_weight_by_charge  /= tot_charge;
      output.em_weight_by_charge     /= tot_charge;
      output.michel_weight_by_charge /= tot_charge;
      output.none_weight_by_charge   /= tot_charge;
    }
    else {
      output.track_weight_by_charge  = -999.;
      output.em_weight_by_charge     = -999.;
      output.michel_weight_by_charge = -999.;
      output.none_weight_by_charge   = -999.;
    }

    return output;
  }
}

//Useful Geant4Reweight util methods
using protoana::G4ReweightUtils::CreateRWTraj;
using protoana::G4ReweightUtils::BuildHierarchy;
using protoana::G4ReweightUtils::CreateNRWTrajs;
using protoana::G4ReweightUtils::GetNTrajWeightFromSetPars;

//Constructor for cnn struct
pduneana::cnnOutput2D::cnnOutput2D() : track(0), em(0), michel(0), none(0), nHits(0), track_weight_by_charge(0), em_weight_by_charge(0), michel_weight_by_charge(0), none_weight_by_charge(0) { }

class pduneana::PDSPAnalyzer : public art::EDAnalyzer {
public:
  explicit PDSPAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPAnalyzer(PDSPAnalyzer const&) = delete;
  PDSPAnalyzer(PDSPAnalyzer&&) = delete;
  PDSPAnalyzer& operator=(PDSPAnalyzer const&) = delete;
  PDSPAnalyzer& operator=(PDSPAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset();
  //Deprecated -- To be removed
  double lateralDist( TVector3 & n, TVector3 & x0, TVector3 & p );

private:

  void BeamTrackInfo(const art::Event & evt,
                     const recob::Track * thisTrack,
                     detinfo::DetectorClocksData const& clockData);
  void BeamShowerInfo(const art::Event & evt, const recob::Shower* thisShower);
  void BeamPFPInfo(const art::Event & evt,
                   const recob::PFParticle* particle,
                   anab::MVAReader<recob::Hit,4> * hitResults);
  void BeamForcedTrackInfo(const art::Event & evt,
                           const recob::PFParticle * particle);
  void CheckEff(const art::Event & evt, int rightTrackID);
  void TrueBeamInfo(const art::Event & evt,
                    const simb::MCParticle* true_beam_particle,
                    detinfo::DetectorClocksData const& clockData,
                    const sim::ParticleList & plist,
                    std::map<int, std::vector<int>> & trueToPFPs,
                    anab::MVAReader<recob::Hit,4> * hitResults);
  void BeamInstInfo(const art::Event & evt);
  void DaughterPFPInfo(const art::Event & evt,
                       const recob::PFParticle* particle,
                       detinfo::DetectorClocksData const& clockData,
                       anab::MVAReader<recob::Hit,4> * hitResults);
  void DaughterMatchInfo(const art::Event & evt,
                         const recob::PFParticle * daughterPFP,
                         detinfo::DetectorClocksData const& clockData);
  void BeamMatchInfo(const art::Event & evt, const recob::PFParticle * particle,
                     const simb::MCParticle * true_beam_particle,
                     detinfo::DetectorClocksData const& clockData);
  void BeamForcedTrackInfo();
  void GetG4RWCoeffs(std::vector<double> & weights, std::vector<double> & coeffs);
  double GetG4RWExpCoeffs(std::vector<double> & weights, std::vector<double> & coeffs);
  void G4RWGridWeights(
      std::vector<std::vector<G4ReweightTraj *>> & hierarchy,
      std::vector<fhicl::ParameterSet> & pars,
      std::vector<std::vector<double>> & weights,
      G4MultiReweighter * multi_rw);
  std::vector<int> PrimaryHierarchy(
      int ID, const sim::ParticleList & plist, bool verbose=false);
  // Declare member data here.
  const art::InputTag fTrackModuleLabel;

  TTree *fTree;
  // Run information
  int run;
  int subrun;
  int event;
  int MC;

  /////////////////////////////////////////////
  //Truth level info of the primary beam particle
  //that generated the event
  int true_beam_PDG;
  double true_beam_mass;
  int true_beam_ID;
  std::vector<int> true_beam_hierarchy;
  int true_beam_hierarchy_size;
  bool reco_beam_true_byHits_in_primary_hierarchy;
  std::string true_beam_endProcess;
  double true_beam_endX;
  double true_beam_endY;
  double true_beam_endZ;
  double true_beam_endX_SCE;
  double true_beam_endY_SCE;
  double true_beam_endZ_SCE;
  double true_beam_startX;
  double true_beam_startY;
  double true_beam_startZ;

  double true_beam_startDirX;
  double true_beam_startDirY;
  double true_beam_startDirZ;

  double true_beam_startPx;
  double true_beam_startPy;
  double true_beam_startPz;
  double true_beam_startP;

  double true_beam_endPx;
  double true_beam_endPy;
  double true_beam_endPz;
  double true_beam_endP, true_beam_endP2, true_beam_last_len;

  int  true_beam_nElasticScatters;
  int  true_beam_nHits;
  std::vector< double > true_beam_elastic_costheta, true_beam_elastic_X,
                        true_beam_elastic_Y, true_beam_elastic_Z,
                        true_beam_elastic_deltaE, true_beam_elastic_IDE_edep;

  double true_beam_IDE_totalDep;
  std::vector< std::string > true_beam_processes;


  std::vector< std::vector< int > > true_beam_reco_byHits_PFP_ID, true_beam_reco_byHits_PFP_nHits,
                                    true_beam_reco_byHits_allTrack_ID;
  //////////////////////////////////////////////////////

  //////////////////////////////////////////////////////
  //Truth level info of the daughter MCParticles coming out of the
  //true primary particle
  std::vector< int > true_beam_daughter_PDG;
  std::vector< int > true_beam_daughter_ID;
  std::vector< double > true_beam_daughter_len;
  std::vector< std::string > true_beam_daughter_Process, true_beam_daughter_endProcess;

  std::vector< double > true_beam_daughter_startX, true_beam_daughter_startY, true_beam_daughter_startZ;
  std::vector< double > true_beam_daughter_startP, true_beam_daughter_startPx, true_beam_daughter_startPy, true_beam_daughter_startPz;
  std::vector< double > true_beam_daughter_endX, true_beam_daughter_endY, true_beam_daughter_endZ;
  std::vector< int >    true_beam_daughter_nHits;


  //going from true to reco byHits
  std::vector< std::vector< int > > true_beam_daughter_reco_byHits_PFP_ID, true_beam_daughter_reco_byHits_PFP_nHits,
                                    true_beam_daughter_reco_byHits_allTrack_ID, true_beam_daughter_reco_byHits_allShower_ID;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_PFP_trackScore;                           
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allTrack_startX, true_beam_daughter_reco_byHits_allTrack_startY, true_beam_daughter_reco_byHits_allTrack_startZ;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allTrack_endX, true_beam_daughter_reco_byHits_allTrack_endY, true_beam_daughter_reco_byHits_allTrack_endZ;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allTrack_len;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allShower_startX, true_beam_daughter_reco_byHits_allShower_startY, true_beam_daughter_reco_byHits_allShower_startZ;
  std::vector< std::vector< double > > true_beam_daughter_reco_byHits_allShower_len;
  //////////////////////////////////////////////////////


  //Decay products from pi0s
  std::vector< int > true_beam_Pi0_decay_PDG, true_beam_Pi0_decay_ID, true_beam_Pi0_decay_parID;
  std::vector< double > true_beam_Pi0_decay_startP, true_beam_Pi0_decay_startPx, true_beam_Pi0_decay_startPy, true_beam_Pi0_decay_startPz;
  std::vector< double > true_beam_Pi0_decay_startX, true_beam_Pi0_decay_startY, true_beam_Pi0_decay_startZ;
  std::vector< int > true_beam_Pi0_decay_nHits;
  std::vector< std::vector< int > > true_beam_Pi0_decay_reco_byHits_PFP_ID, true_beam_Pi0_decay_reco_byHits_PFP_nHits,
                                    true_beam_Pi0_decay_reco_byHits_allTrack_ID, true_beam_Pi0_decay_reco_byHits_allShower_ID;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_PFP_trackScore;                           
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allTrack_startX, true_beam_Pi0_decay_reco_byHits_allTrack_startY, true_beam_Pi0_decay_reco_byHits_allTrack_startZ;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allTrack_endX, true_beam_Pi0_decay_reco_byHits_allTrack_endY, true_beam_Pi0_decay_reco_byHits_allTrack_endZ;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allTrack_len;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allShower_startX, true_beam_Pi0_decay_reco_byHits_allShower_startY, true_beam_Pi0_decay_reco_byHits_allShower_startZ;
  std::vector< std::vector< double > > true_beam_Pi0_decay_reco_byHits_allShower_len;
  //also reco nhits
  std::vector< double > true_beam_Pi0_decay_len;

  std::vector< int > true_beam_grand_daughter_PDG, true_beam_grand_daughter_ID, true_beam_grand_daughter_parID;
  std::vector< int > true_beam_grand_daughter_nHits;
  std::vector< std::string > true_beam_grand_daughter_Process, true_beam_grand_daughter_endProcess;

  //How many of each true particle came out of the true primary beam particle?
  int true_daughter_nPiPlus, true_daughter_nPiMinus, true_daughter_nPi0;
  int true_daughter_nProton, true_daughter_nNeutron, true_daughter_nNucleus;

  //Matched to vertex/slice?
  //
  int reco_beam_vertex_slice;
  ////////////////////////


  //Reconstructed track info
  //EDIT: STANDARDIZE
  double reco_beam_startX, reco_beam_startY, reco_beam_startZ;
  double reco_beam_endX, reco_beam_endY, reco_beam_endZ;
  int test_branch;
  double true_beam_len;
  double reco_beam_len, reco_beam_alt_len;
  int for_truncation_method;
  double reco_beam_alt_len_allTrack;
  double reco_beam_vertex_michel_score;
  int reco_beam_vertex_nHits;
  double reco_beam_vertex_michel_score_allTrack;
  int reco_beam_vertex_nHits_allTrack;
  double reco_beam_vertex_michel_score_weight_by_charge;
  double reco_beam_vertex_michel_score_weight_by_charge_allTrack;

  //position from SCE corrected calo
  double reco_beam_calo_startX, reco_beam_calo_startY, reco_beam_calo_startZ;
  double reco_beam_calo_endX, reco_beam_calo_endY, reco_beam_calo_endZ;
  double reco_beam_calo_startX_allTrack, reco_beam_calo_startY_allTrack, reco_beam_calo_startZ_allTrack;
  double reco_beam_calo_endX_allTrack, reco_beam_calo_endY_allTrack, reco_beam_calo_endZ_allTrack;
  std::vector<double> reco_beam_calo_X, reco_beam_calo_Y, reco_beam_calo_Z;
  std::vector<double> reco_beam_calo_X_allTrack, reco_beam_calo_Y_allTrack, reco_beam_calo_Z_allTrack;
  std::vector<double> reco_beam_calo_startDirX, reco_beam_calo_endDirX;
  std::vector<double> reco_beam_calo_startDirY, reco_beam_calo_endDirY;
  std::vector<double> reco_beam_calo_startDirZ, reco_beam_calo_endDirZ;

  double reco_beam_trackDirX, reco_beam_trackDirY, reco_beam_trackDirZ;
  double reco_beam_trackEndDirX, reco_beam_trackEndDirY, reco_beam_trackEndDirZ;

  std::vector<double> reco_beam_dEdX_SCE, reco_beam_dQdX_SCE, reco_beam_EField_SCE, reco_beam_resRange_SCE, reco_beam_TrkPitch_SCE;
  std::vector<double> reco_beam_TrkPitch_SCE_allTrack;
  std::vector<double> reco_beam_calibrated_dEdX_SCE, reco_beam_calibrated_dQdX_SCE, reco_beam_dQ;

  std::vector<std::vector<int>> reco_beam_hit_IDE_IDs, reco_beam_hit_IDE_origins;
  std::vector<std::vector<double>> reco_beam_hit_IDE_electrons, reco_beam_hit_IDE_energies;
  std::vector<double> reco_beam_hit_IDE_cosmic_electrons, reco_beam_hit_IDE_beam_electrons;
  std::vector<double> reco_beam_hit_IDE_cosmic_energies, reco_beam_hit_IDE_beam_energies;

  std::vector<double> reco_beam_dEdX_NoSCE, reco_beam_dQdX_NoSCE, reco_beam_resRange_NoSCE, reco_beam_TrkPitch_NoSCE;
  std::vector<double> reco_beam_calibrated_dEdX_NoSCE;

  std::vector<double> reco_beam_calo_wire, reco_beam_calo_tick, reco_beam_calo_wire_z;
  std::vector<double> reco_beam_calo_wire_allTrack;
  std::vector<double> reco_beam_calo_wire_NoSCE, reco_beam_dQ_NoSCE, reco_beam_calo_wire_z_NoSCE;
  std::vector<int> reco_beam_calo_TPC, reco_beam_calo_TPC_NoSCE;

  int reco_beam_trackID;
  int n_beam_slices, n_beam_particles;
  std::vector<int> beam_track_IDs;
  std::vector<double> beam_particle_scores;
  bool reco_beam_flipped;

  //fix
  bool reco_beam_passes_beam_cuts;              

  int reco_beam_type;
  double reco_beam_Chi2_proton, reco_beam_Chi2_muon;
  int    reco_beam_Chi2_ndof;

  ////////////////////////

  //For all track info 
  std::vector<double> reco_track_startX, reco_track_startY, reco_track_startZ,
                      reco_track_endX, reco_track_endY, reco_track_endZ,
                      reco_track_michel_score, reco_track_michel_score_weight_by_charge;
  std::vector<int> reco_track_ID, reco_track_nHits;


  //GeantReweight stuff
  // -- Maybe think of new naming scheme?
  //std::vector<double> g4rw_primary_plus_sigma_weight;
  //std::vector<double> g4rw_primary_minus_sigma_weight;
  //std::vector<std::string> g4rw_primary_var;

  //std::vector<double> g4rw_alt_primary_plus_sigma_weight;
  //std::vector<double> g4rw_alt_primary_minus_sigma_weight;
  //std::vector<double> g4rw_full_primary_plus_sigma_weight;
  //std::vector<double> g4rw_full_primary_minus_sigma_weight;

  //std::vector<std::vector<double>> g4rw_full_grid_weights,
  //                                 g4rw_full_grid_coeffs;
  std::vector<std::vector<double>> g4rw_primary_grid_weights,
                                   g4rw_primary_grid_coeffs;
  //std::vector<double> g4rw_primary_grid_pair_weights;

  std::vector<std::vector<double>> g4rw_piplus_traj_ps,
                                   g4rw_piplus_traj_lens;
  std::vector<int> g4rw_piplus_traj_npiplus,
                   g4rw_piplus_traj_npiminus,
                   g4rw_piplus_traj_npi0;
  std::vector<std::vector<double>> g4rw_full_grid_piplus_weights,
                                   g4rw_full_grid_piplus_coeffs;
  std::vector<std::vector<double>> g4rw_full_grid_piplus_weights_fake_data,
                                   g4rw_full_grid_piplus_coeffs_fake_data;
  std::vector<std::vector<double>> g4rw_full_grid_piminus_weights;
  std::vector<std::vector<double>> g4rw_full_grid_proton_weights,
                                   g4rw_full_grid_proton_coeffs;
  std::vector<std::vector<double>> g4rw_full_grid_neutron_weights,
                                   g4rw_full_grid_neutron_coeffs;
  std::vector<std::vector<double>> g4rw_full_grid_kplus_weights,
                                   g4rw_full_grid_kplus_coeffs;

  std::vector<std::vector<double>> g4rw_full_grid_abscex_weights,
                                   g4rw_full_grid_abscex_coeffs,
                                   g4rw_primary_grid_abscex_weights,
                                   g4rw_primary_grid_abscex_coeffs,
                                   g4rw_downstream_grid_abscex_weights,
                                   g4rw_downstream_grid_abscex_coeffs;

  std::vector<std::vector<double>> g4rw_downstream_grid_piplus_weights,
                                   g4rw_downstream_grid_piplus_coeffs;

  std::vector<std::vector<double>> g4rw_full_fine_piplus_weights,
                                   g4rw_full_fine_piplus_coeffs;

  std::vector<double> g4rw_primary_grid_exp_fit_chi2,
                      g4rw_full_grid_piplus_exp_fit_chi2,
                      g4rw_full_grid_piplus_exp_fit_chi2_fake_data,
                      g4rw_full_grid_proton_exp_fit_chi2,
                      g4rw_full_grid_neutron_exp_fit_chi2,
                      g4rw_full_grid_kplus_exp_fit_chi2,
                      g4rw_downstream_grid_piplus_exp_fit_chi2,
                      g4rw_full_fine_piplus_exp_fit_chi2;

  //EDIT: STANDARDIZE
  //EndProcess --> endProcess ?
  std::string reco_beam_true_byE_endProcess, reco_beam_true_byHits_endProcess; //What process ended the reco beam particle
  std::string reco_beam_true_byE_process, reco_beam_true_byHits_process;    //What process created the reco beam particle
  int reco_beam_true_byE_PDG, reco_beam_true_byHits_PDG;
  int reco_beam_true_byE_ID, reco_beam_true_byHits_ID;
  bool reco_beam_true_byE_matched, reco_beam_true_byHits_matched; //Does the true particle contributing most to the
                                           //reconstructed beam track coincide with the actual
                                           //beam particle that generated the event
  int reco_beam_true_byE_origin, reco_beam_true_byHits_origin; //What is the origin of the reconstructed beam track?
  //EDIT: STANDARDIZE
  //End_P --> endP, etc.
  double reco_beam_true_byE_endPx,   reco_beam_true_byHits_endPx;
  double reco_beam_true_byE_endPy,   reco_beam_true_byHits_endPy;
  double reco_beam_true_byE_endPz,   reco_beam_true_byHits_endPz;
  double reco_beam_true_byE_endE,    reco_beam_true_byHits_endE;
  double reco_beam_true_byE_endP,    reco_beam_true_byHits_endP;
                          
  double reco_beam_true_byE_startPx, reco_beam_true_byHits_startPx;
  double reco_beam_true_byE_startPy, reco_beam_true_byHits_startPy;
  double reco_beam_true_byE_startPz, reco_beam_true_byHits_startPz;
  double reco_beam_true_byE_startE,  reco_beam_true_byHits_startE;
  double reco_beam_true_byE_startP,  reco_beam_true_byHits_startP;
  //also throw in byE
  double reco_beam_true_byHits_purity;             
  //////////////////////////

  std::vector< double > reco_beam_incidentEnergies;
  double reco_beam_interactingEnergy;
  std::vector< double > reco_beam_incidentEnergies_allTrack;
  double reco_beam_interactingEnergy_allTrack;
  std::vector< double > true_beam_incidentEnergies;
  std::vector< int >    true_beam_slices, true_beam_slices_found;
  std::vector< double > true_beam_slices_deltaE;
  double true_beam_interactingEnergy;
  double em_energy;
  std::vector<double> true_beam_traj_X;
  std::vector<double> true_beam_traj_Y;
  std::vector<double> true_beam_traj_Z;
  std::vector<double> true_beam_traj_Px;
  std::vector<double> true_beam_traj_Py;
  std::vector<double> true_beam_traj_Pz;
  std::vector<double> true_beam_traj_KE;
  std::vector<double> true_beam_traj_X_SCE;
  std::vector<double> true_beam_traj_Y_SCE;
  std::vector<double> true_beam_traj_Z_SCE;
  bool true_beam_is_scraper;

  int    reco_beam_PFP_ID;
  int    reco_beam_PFP_nHits;
  double reco_beam_PFP_trackScore, reco_beam_PFP_trackScore_weight_by_charge;
  double reco_beam_PFP_emScore, reco_beam_PFP_emScore_weight_by_charge;
  double reco_beam_PFP_michelScore, reco_beam_PFP_michelScore_weight_by_charge;
  double reco_beam_PFP_trackScore_collection, reco_beam_PFP_trackScore_collection_weight_by_charge;
  double reco_beam_PFP_emScore_collection, reco_beam_PFP_emScore_collection_weight_by_charge;
  double reco_beam_PFP_michelScore_collection, reco_beam_PFP_michelScore_collection_weight_by_charge;

  int    reco_beam_allTrack_ID;
  bool   reco_beam_allTrack_beam_cuts, reco_beam_allTrack_flipped;
  double reco_beam_allTrack_len;
  double reco_beam_allTrack_startX, reco_beam_allTrack_startY, reco_beam_allTrack_startZ;
  double reco_beam_allTrack_endX, reco_beam_allTrack_endY, reco_beam_allTrack_endZ;
  double reco_beam_allTrack_trackDirX, reco_beam_allTrack_trackDirY, reco_beam_allTrack_trackDirZ;
  double reco_beam_allTrack_trackEndDirX, reco_beam_allTrack_trackEndDirY, reco_beam_allTrack_trackEndDirZ;
  std::vector< double > reco_beam_allTrack_resRange;
  std::vector< double > reco_beam_allTrack_calibrated_dEdX;
  double reco_beam_allTrack_Chi2_proton;
  int    reco_beam_allTrack_Chi2_ndof;

  /////////////////////////////////////////////////////
  //Info from the BI
  /////////////////////////////////////////////////////
  double beam_inst_P;
  bool beam_inst_valid;
  int beam_inst_trigger;
  std::vector<double> beam_inst_TOF;
  std::vector< int > beam_inst_PDG_candidates, beam_inst_TOF_Chan;
  double beam_inst_X, beam_inst_Y, beam_inst_Z;
  double beam_inst_dirX, beam_inst_dirY, beam_inst_dirZ;
  int beam_inst_nFibersP1, beam_inst_nFibersP2, beam_inst_nFibersP3;
  int beam_inst_nTracks, beam_inst_nMomenta;
  int beam_inst_C0, beam_inst_C1;
  double beam_inst_C0_pressure, beam_inst_C1_pressure;
  ////////////////////////////////////////////////////

  bool reco_reconstructable_beam_event;

  //Alternative Reco values
  //EDIT: track_score --> trkScore, etc.
  std::vector<int> reco_daughter_PFP_ID;
  std::vector<int> reco_daughter_pandora_type;
  std::vector<int> reco_daughter_PFP_nHits,
                   reco_daughter_PFP_nHits_collection;
  std::vector< double > reco_daughter_PFP_trackScore;
  std::vector< double > reco_daughter_PFP_emScore;
  std::vector< double > reco_daughter_PFP_michelScore;
  std::vector< double > reco_daughter_PFP_trackScore_collection;
  std::vector< double > reco_daughter_PFP_emScore_collection;
  std::vector< double > reco_daughter_PFP_michelScore_collection;

  std::vector< double > reco_daughter_PFP_trackScore_weight_by_charge;
  std::vector< double > reco_daughter_PFP_emScore_weight_by_charge;
  std::vector< double > reco_daughter_PFP_michelScore_weight_by_charge;
  std::vector< double > reco_daughter_PFP_trackScore_collection_weight_by_charge;
  std::vector< double > reco_daughter_PFP_emScore_collection_weight_by_charge;
  std::vector< double > reco_daughter_PFP_michelScore_collection_weight_by_charge;


  //EDIT: reco_daughter_PFP_true_byY_XXX
  std::vector< int > reco_daughter_PFP_true_byHits_PDG;
  std::vector< int > reco_daughter_PFP_true_byHits_ID;
  std::vector< int > reco_daughter_PFP_true_byHits_origin;
  std::vector< int > reco_daughter_PFP_true_byHits_parID;
  std::vector< int > reco_daughter_PFP_true_byHits_parPDG;
  std::vector< std::string > reco_daughter_PFP_true_byHits_process;
  std::vector< double > reco_daughter_PFP_true_byHits_purity;///EDIT: quality
  std::vector< size_t > reco_daughter_PFP_true_byHits_sharedHits, reco_daughter_PFP_true_byHits_emHits;
  std::vector< double > reco_daughter_PFP_true_byHits_completeness;

  std::vector< double > reco_daughter_PFP_true_byHits_len;
  std::vector< double > reco_daughter_PFP_true_byHits_startX;
  std::vector< double > reco_daughter_PFP_true_byHits_startY;
  std::vector< double > reco_daughter_PFP_true_byHits_startZ;
  std::vector< double > reco_daughter_PFP_true_byHits_endX;
  std::vector< double > reco_daughter_PFP_true_byHits_endY;
  std::vector< double > reco_daughter_PFP_true_byHits_endZ;

  std::vector< double > reco_daughter_PFP_true_byHits_startPx;
  std::vector< double > reco_daughter_PFP_true_byHits_startPy;
  std::vector< double > reco_daughter_PFP_true_byHits_startPz;
  std::vector< double > reco_daughter_PFP_true_byHits_startE;
  std::vector< double > reco_daughter_PFP_true_byHits_startP;

  std::vector< std::string > reco_daughter_PFP_true_byHits_endProcess;

  std::vector< int > reco_daughter_PFP_true_byE_PDG;
  std::vector< double > reco_daughter_PFP_true_byE_len;
  std::vector< double > reco_daughter_PFP_true_byE_completeness;
  std::vector< double > reco_daughter_PFP_true_byE_purity;



  //////////////////////////////////////

  //EDIT: reco_daughter_allTrack_XXX
  std::vector< int > reco_daughter_allTrack_ID;
  std::vector< double > reco_daughter_allTrack_Theta;
  std::vector< double > reco_daughter_allTrack_Phi;
  std::vector<double> reco_daughter_allTrack_startDirX,
                      reco_daughter_allTrack_startDirY,
                      reco_daughter_allTrack_startDirZ;
  std::vector< std::vector< double > > reco_daughter_allTrack_dQdX_SCE, reco_daughter_allTrack_dEdX_SCE, reco_daughter_allTrack_resRange_SCE;
  std::vector< std::vector< double > > reco_daughter_allTrack_calibrated_dEdX_SCE, reco_daughter_allTrack_calibrated_dQdX_SCE, reco_daughter_allTrack_EField_SCE, reco_daughter_allTrack_calo_X, reco_daughter_allTrack_calo_Y, reco_daughter_allTrack_calo_Z;
  std::vector< double > reco_daughter_allTrack_Chi2_proton, reco_daughter_allTrack_Chi2_muon, reco_daughter_allTrack_Chi2_pion;
  std::vector< int >    reco_daughter_allTrack_Chi2_ndof, reco_daughter_allTrack_Chi2_ndof_muon, reco_daughter_allTrack_Chi2_ndof_pion;

  //New: calorimetry + chi2 for planes 0 and 1
  std::vector<std::vector<double>>
      reco_daughter_allTrack_calibrated_dEdX_SCE_plane0,
      reco_daughter_allTrack_calibrated_dEdX_SCE_plane1;

  std::vector<std::vector<double>>
      reco_daughter_allTrack_resRange_plane0,
      reco_daughter_allTrack_resRange_plane1;

  std::vector<double> reco_daughter_allTrack_Chi2_proton_plane0,
                      reco_daughter_allTrack_Chi2_proton_plane1;

  std::vector<int> reco_daughter_allTrack_Chi2_ndof_plane0,
                   reco_daughter_allTrack_Chi2_ndof_plane1;
  //////////////////////////////////////////////

  std::vector< double > reco_daughter_allTrack_startX, reco_daughter_allTrack_endX;
  std::vector< double > reco_daughter_allTrack_startY, reco_daughter_allTrack_endY;
  std::vector< double > reco_daughter_allTrack_startZ, reco_daughter_allTrack_endZ;
  std::vector< double > reco_daughter_allTrack_len, reco_daughter_allTrack_alt_len;

  std::vector<double> reco_daughter_allTrack_vertex_michel_score;
  std::vector<int> reco_daughter_allTrack_vertex_nHits;
  //

  std::vector<int>    reco_daughter_allShower_ID;
  std::vector<double> reco_daughter_allShower_len,
                      reco_daughter_allShower_startX,
                      reco_daughter_allShower_startY,
                      reco_daughter_allShower_startZ,
                      reco_daughter_allShower_dirX,
                      reco_daughter_allShower_dirY,
                      reco_daughter_allShower_dirZ,
                      reco_daughter_allShower_energy,
                      reco_daughter_allShower_calibrated_energy;


  std::vector<double> reco_daughter_allTrack_momByRange_proton;
  std::vector<double> reco_daughter_allTrack_momByRange_muon;
  double reco_beam_momByRange_proton;
  double reco_beam_momByRange_muon;

  std::vector<double> reco_daughter_allTrack_momByRange_alt_proton;
  std::vector<double> reco_daughter_allTrack_momByRange_alt_muon;
  double reco_beam_momByRange_alt_proton;
  double reco_beam_momByRange_alt_muon;

  //New hits info for SparseNet
  std::vector< double > reco_beam_spacePts_X, reco_beam_spacePts_Y, reco_beam_spacePts_Z;
  std::vector< std::vector< double > > reco_daughter_spacePts_X, reco_daughter_spacePts_Y, reco_daughter_spacePts_Z;
  std::vector< std::vector< double > > reco_daughter_shower_spacePts_X, reco_daughter_shower_spacePts_Y, reco_daughter_shower_spacePts_Z;
 
  std::vector< std::vector< double > > sparsenet_features_charge;
  std::vector< std::vector< double > > sparsenet_features_angle;
  std::vector< std::vector< double > > sparsenet_features_dot_product;
  std::vector< std::vector< double > > sparsenet_features_neighboring_nodes_3;
  std::vector< std::vector< double > > sparsenet_features_neighboring_nodes_10;
  std::vector< std::vector< double > > sparsenet_features_neighboring_nodes_30;
  std::vector< std::vector< double > > sparsenet_features_charge_distance_3;
  std::vector< std::vector< double > > sparsenet_features_charge_distance_10;
  std::vector< std::vector< double > > sparsenet_features_charge_distance_30;
  



  ////New section -- mechanical class members
  std::map< int, TProfile* > templates;

  //FCL pars
  std::string fCalorimetryTagSCE;
  std::string fPandora2CaloSCE;
  std::string fCalorimetryTagNoSCE;
  std::string fTrackerTag;
  std::string fHitTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fBeamModuleLabel;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  protoana::ProtoDUNEEmptyEventFinder fEmptyEventFinder;
  double fBeamPIDMomentum;
  std::string dEdX_template_name;
  TFile * dEdX_template_file;
  bool fVerbose;    
  fhicl::ParameterSet BeamPars;
  fhicl::ParameterSet BeamCuts;
  protoana::ProtoDUNEBeamCuts beam_cuts;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNEShowerUtils showerUtil;
  protoana::ProtoDUNETruthUtils truthUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  art::ServiceHandle<geo::Geometry> geom;
  geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  protoana::ProtoDUNECalibration calibration_SCE;
  protoana::ProtoDUNECalibration calibration_NoSCE;
  bool fSaveHits;
  bool fSaveHitIDEInfo;
  bool fSaveRecoBeamAllTrack;
  bool fSkipMVA;
  bool fTrueToReco;
  bool fDoReweight;
  bool fSaveG4RWWeights;
  bool fDoProtReweight;
  bool fGetTrackMichel;
  bool fMCHasBI;
  bool fRecalibrate;
  bool fGetCalibratedShowerEnergy;
  bool fGetUncalibratedShowerEnergy;
  bool fCheckSlicesForBeam;
  bool fCheckTruncation;
  double fBeamInstPFix = 1.;
  // SparseNet params 
  /// Module label for input space points
  std::string fSpacePointLabel;
  /// Radii for calculating number of neighbours for any number of cut values
  std::vector<float> fNeighbourRadii;
  /// Radii for calculating charge over distance
  std::vector<float> fChargeRadii;

  bool fSCE;

  double fZ0, fPitch;

  //Geant4Reweight stuff
  TFile * FracsFile/*, * PiMinusFracsFile*/;
  TFile * ProtFracsFile;
  std::vector<fhicl::ParameterSet> ParSet, FakeDataParSet, ProtParSet,
                                   fAbsCexPars, FineParSet;//, PiMinusParSet;
  std::vector<double> fGridPoints;
  std::pair<double, double> fGridPair;
  //G4ReweightParameterMaker ParMaker, FakeDataParameterMaker, ProtParMaker;//, PiMinusParMaker;
  G4ReweightParameterMaker ParMaker;
  AbsCexReweighter * fAbsCex_reweighter;
  G4MultiReweighter * MultiRW, * ProtMultiRW, * PiMinusMultiRW,
                    * KPlusMultiRW, * NeutronMultiRW, * FakeDataMultiRW,
                    * FineMultiRW, * fAbsCexMultiRW;
  G4ReweightManager * RWManager;
};


pduneana::PDSPAnalyzer::PDSPAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  ,
  fTrackModuleLabel(p.get< art::InputTag >("TrackModuleLabel")),
  fCalorimetryTagSCE(p.get<std::string>("CalorimetryTagSCE")),
  fPandora2CaloSCE(p.get<std::string>("Pandora2CaloSCE")),
  fCalorimetryTagNoSCE(p.get<std::string>("CalorimetryTagNoSCE")),
  fTrackerTag(p.get<std::string>("TrackerTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fBeamModuleLabel(p.get<std::string>("BeamModuleLabel")),
  fBeamlineUtils(p.get< fhicl::ParameterSet >("BeamlineUtils")),
  fEmptyEventFinder(p.get< fhicl::ParameterSet >("EmptyEventFinder")),
  fBeamPIDMomentum(p.get<double>("BeamPIDMomentum", 1.)),
  dEdX_template_name(p.get<std::string>("dEdX_template_name")),
  //dEdX_template_file( dEdX_template_name.c_str(), "OPEN" ),
  fVerbose(p.get<bool>("Verbose")),
  BeamPars(p.get<fhicl::ParameterSet>("BeamPars")),
  BeamCuts(p.get<fhicl::ParameterSet>("BeamCuts")),
  calibration_SCE(p.get<fhicl::ParameterSet>("CalibrationParsSCE")),
  calibration_NoSCE(p.get<fhicl::ParameterSet>("CalibrationParsNoSCE")),
  fSaveHits( p.get<bool>( "SaveHits" ) ),
  fSaveHitIDEInfo(p.get<bool>( "SaveHitIDEInfo", false)),
  fSaveRecoBeamAllTrack(p.get<bool>("SaveRecoBeamAllTrack", false)),
  fSkipMVA( p.get<bool>( "SkipMVA" ) ),
  fTrueToReco( p.get<bool>( "TrueToReco" ) ),
  fDoReweight(p.get<bool>("DoReweight")),
  fSaveG4RWWeights(p.get<bool>("SaveG4RWWeights", true)),
  fDoProtReweight(p.get<bool>("DoProtReweight")),
  fGetTrackMichel(p.get<bool>("GetTrackMichel")),
  fRecalibrate(p.get<bool>("Recalibrate", true)),
  fGetCalibratedShowerEnergy(p.get<bool>("GetCalibratedShowerEnergy", false)),
  fGetUncalibratedShowerEnergy(p.get<bool>("GetUncalibratedShowerEnergy", true)),
  fCheckSlicesForBeam(p.get<bool>("CheckSlicesForBeam", false)),
  fCheckTruncation(p.get<bool>("CheckTruncation", false)),

  fBeamInstPFix(p.get<double>("BeamInstPFix", 1.)),

  // SparseNet params
  fSpacePointLabel (p.get<std::string>    ("SpacePointLabel")), 
  fNeighbourRadii (p.get<std::vector<float>>("NeighbourRadii")),
  fChargeRadii (p.get<std::vector<float>>("ChargeRadii")),


  fSCE(p.get<bool>("SCE", true)){

  dEdX_template_file = OpenFile(dEdX_template_name);
  templates[ 211 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_pi"  );
  templates[ 321 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_ka"  );
  templates[ 13 ]   = (TProfile*)dEdX_template_file->Get( "dedx_range_mu"  );
  templates[ 2212 ] = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );

  beam_cuts = protoana::ProtoDUNEBeamCuts( BeamCuts );

  if (fRecalibrate) {
    calibration_SCE = p.get<fhicl::ParameterSet>("CalibrationParsSCE");
    calibration_NoSCE = p.get<fhicl::ParameterSet>("CalibrationParsNoSCE");
  }


  if (fDoReweight) {
    //Pion Reweighter
    auto material = p.get<fhicl::ParameterSet>("Material");
    RWManager = new G4ReweightManager({material});
    FracsFile = OpenFile(p.get< std::string >("FracsFile"));
    ParSet = p.get<std::vector<fhicl::ParameterSet>>("ParameterSet");
    FineParSet = p.get<std::vector<fhicl::ParameterSet>>("FineParameterSet");
    //ParMaker = G4ReweightParameterMaker(ParSet);
    MultiRW = new G4MultiReweighter(211, *FracsFile, ParSet,
                                    material,
                                    RWManager);
    FineMultiRW = new G4MultiReweighter(211, *FracsFile, FineParSet,
                                    material,
                                    RWManager);
    FakeDataParSet = p.get<std::vector<fhicl::ParameterSet>>("FakeDataParameterSet");
    FakeDataMultiRW = new G4MultiReweighter(211, *FracsFile, FakeDataParSet,
                                            material,
                                            RWManager);
    //Piminus Reweighter
    //PiMinusFracsFile = OpenFile(p.get< std::string >("PiMinusFracsFile"));
    //PiMinusMultiRW = new G4MultiReweighter(-211, *PiMinusFracsFile, ParSet,
    //                                material,
    //                                RWManager);
    
    //Proton Reweighter
    ProtFracsFile = OpenFile(p.get< std::string >("ProtFracsFile"));
    ProtParSet = p.get<std::vector<fhicl::ParameterSet>>("ProtParameterSet");
    //ProtParMaker = G4ReweightParameterMaker(ProtParSet);
    ProtMultiRW = new G4MultiReweighter(
        2212, *ProtFracsFile, ProtParSet,
        material,
        RWManager);
    

    //K Plus Reweighter
    // -- Just use proton fracs and parameter set
    // -- because it's just reaction
    KPlusMultiRW = new G4MultiReweighter(
        321, *ProtFracsFile, ProtParSet,
        material,
        RWManager);

    //Neutron Reweighter
    // -- Just use proton fracs and parameter set
    // -- because it's just reaction
    NeutronMultiRW = new G4MultiReweighter(
        2112, *ProtFracsFile, ProtParSet,
        material,
        RWManager);
    
    TFile * abscex_file = TFile::Open(p.get<std::string>("AbsCexFile").c_str());
    fAbsCexPars
        = p.get<std::vector<fhicl::ParameterSet>>("AbsCexParameterSet");
    ParMaker = G4ReweightParameterMaker(fAbsCexPars, {"abs", "cex", "other"});

    fAbsCex_reweighter = new AbsCexReweighter(
        abscex_file, ParMaker.GetFSHists(),
        material,
        RWManager,
        ParMaker.GetElasticHist());
    fAbsCexMultiRW = new G4MultiReweighter(
        211, *abscex_file, fAbsCexPars, ParMaker, material, fAbsCex_reweighter);

    //Setting up grid points for mutliple variations
    double start = p.get<double>("ParameterGridStart");
    double delta = p.get<double>("ParameterGridDelta");
    int n = p.get<int>("ParameterGridN"); 
    for (int i = 0; i < n; ++i) {
      fGridPoints.push_back(start);
      start += delta;
    }
    fGridPair = p.get<std::pair<double, double>>("GridPair");
  }

  constexpr geo::PlaneID planeID{0, 1, 2};
  fZ0 = wireReadout.Wire( geo::WireID(planeID, 0) ).GetCenter().Z();
  fPitch = wireReadout.Plane(planeID).WirePitch();

}

void pduneana::PDSPAnalyzer::analyze(art::Event const & evt) {

  //reset containers
  reset();


  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  std::cout<<"########## EvtNo."<<event<<std::endl;

  if( !evt.isRealData() ) MC = 1;
  else MC = 0;

  // Is this a reconstructable beam event?
  reco_reconstructable_beam_event = !fEmptyEventFinder.IsEmptyEvent(evt);

  // Get various utilities
  protoana::ProtoDUNEPFParticleUtils                    pfpUtil;
  auto pfpVec = evt.getValidHandle< std::vector< recob::PFParticle > >( fPFParticleTag );

  if (fGetTrackMichel && !fSkipMVA) {
    for (const recob::PFParticle & pfp : (*pfpVec)) {
      const recob::Track* tempTrack = pfpUtil.GetPFParticleTrack(pfp, evt,
                                                                 fPFParticleTag,
                                                                 fTrackerTag);
      if (tempTrack) {
        auto const start = tempTrack->Start();
        double startX = start.X();
        double startY = start.Y();
        double startZ = start.Z();

        auto const end = tempTrack->End();
        double endX = end.X();
        double endY = end.Y();
        double endZ = end.Z();

        int end_tpc = geom->FindTPCAtPosition(end).TPC;
        int start_tpc = geom->FindTPCAtPosition(start).TPC;

        if (!((end_tpc == 1 || end_tpc == 5) &&
              (start_tpc == 1 || start_tpc == 5)))
          continue;

        std::pair<double, int> vertex_michel_score =
            trackUtil.GetVertexMichelScore(*tempTrack, evt, fTrackerTag,
                                           fHitTag);
        std::pair<double, double> vertex_michel_score_weight_by_charge =
            trackUtil.GetVertexMichelScore_weight_by_charge(*tempTrack, evt, fTrackerTag,
                                           fHitTag);

        reco_track_michel_score.push_back(vertex_michel_score.first);
        reco_track_nHits.push_back(vertex_michel_score.second);
        reco_track_michel_score_weight_by_charge.push_back( vertex_michel_score_weight_by_charge.second != 0 ? vertex_michel_score_weight_by_charge.first/vertex_michel_score_weight_by_charge.second : -999. );
        reco_track_ID.push_back(tempTrack->ID());
        reco_track_startX.push_back(startX);
        reco_track_startY.push_back(startY);
        reco_track_startZ.push_back(startZ);
        reco_track_endX.push_back(endX);
        reco_track_endY.push_back(endY);
        reco_track_endZ.push_back(endZ);
      }
    }
  }

  const sim::ParticleList & plist = pi_serv->ParticleList();

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp =  art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);
  trkf::TrackMomentumCalculator track_p_calc;
  ////////////////////////////////////////


  //double z0 = wireReadout.Wire( geo::WireID(0, 1, 2, 0) ).GetCenter().Z();
  //double pitch = wireReadout.WirePitch( 2, 1, 0);

  if (fVerbose) {
    std::cout << "Z0: " << fZ0 << std::endl;
    std::cout << "Pitch: " << fPitch << std::endl;

    double z0_APA2 = wireReadout.Wire(geo::WireID(0, 5, 2, 0)).GetCenter().Z();
    std::cout << "APA 2 Z0: " << z0_APA2 << std::endl;
  }

  // This gets the true beam particle that generated the event
  const simb::MCParticle* true_beam_particle = 0x0;
  if( !evt.isRealData() ){
    auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
    true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);
    if( !true_beam_particle ){
      MF_LOG_WARNING("PDSPAnalyzer") << "No true beam particle" << std::endl;
      return;
    }
    if (fVerbose) {
      std::cout << "Got " << (*mcTruths)[0].NParticles() <<
                   " particles in mcTruth" << std::endl;
      for (int i = 0; i < (*mcTruths)[0].NParticles(); ++i) {
        simb::MCParticle part = (*mcTruths)[0].GetParticle(i);
        std::cout << part.Process() << " " << part.TrackId() << " " <<
                     part.PdgCode() << std::endl;

      }
    }
  }
  ////////////////////////////

  //Get the beam instrumentation from the event
  BeamInstInfo(evt);


  // Helper to get hits and the 4 associated CNN outputs
  // CNN Outputs: EM, Track, Michel, Empty
  // outputNames: track, em, none, michel
  anab::MVAReader<recob::Hit,4> * hitResults = 0x0;
  if (!fSkipMVA) { 
    hitResults = new anab::MVAReader<recob::Hit, 4>(evt, "emtrkmichelid:emtrkmichel" );
  }

  /*
  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fTrackerTag);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,fTrackerTag);

  auto recoShowers = evt.getValidHandle< std::vector< recob::Shower > >(fShowerTag);
  art::FindManyP<recob::Hit> findHitsFromShowers(recoShowers,evt,fShowerTag);
  */

  std::map< int, std::vector< int > > trueToPFPs;
  if( fTrueToReco ){
    trueToPFPs = truthUtil.GetMapMCToPFPs_ByHits( clockData, evt, fPFParticleTag, fHitTag );
  }

  // Processing beam slices
  n_beam_slices = 0;
  n_beam_particles = 0;
  const std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap
      = pfpUtil.GetPFParticleSliceMap(evt, fPFParticleTag);
  std::vector<std::vector<const recob::PFParticle*>> beam_slices;
  for (auto slice : sliceMap) {
    for (auto particle : slice.second) {
      bool added = false;
      if (pfpUtil.IsBeamParticle(*particle,evt, fPFParticleTag)) {
        if (!added) {
          beam_slices.push_back(slice.second);
          ++n_beam_slices;
          added = true;
        }
        //std::cout << "Slice: " << slice.first << " N Part: " <<
        //             slice.second.size() << std::endl;
        for (const auto * p : slice.second) {
          const recob::Track* track
              = pfpUtil.GetPFParticleTrack(*p,evt,fPFParticleTag,fTrackerTag);
          ++n_beam_particles;
          beam_particle_scores.push_back(pfpUtil.GetBeamCosmicScore(*p, evt, fPFParticleTag));

          if (track) {
            //std::cout << "Slice: " << slice.first << " ID: " << track->ID() <<
            //             std::endl;
            beam_track_IDs.push_back(track->ID());
          }
          else {
            beam_track_IDs.push_back(-999);
          }
        }

      }
      if (!added) {continue;}
      //else {
      //  std::cout << "Not beam particle" << std::endl;
      //}
    }
  }
  //std::cout << "Got " << beam_slices.size() <<" beam slices" << std::endl;

  int rightTrackID = -1;
  ///Gets the beam pfparticle
  if(beam_slices.size() == 0){
    std::cout << "We found no beam particles for this event... moving on" << std::endl;
    //return;
  }
  else {
    std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt,fPFParticleTag);
    if (fVerbose)
      std::cout << "Found " << beamParticles.size() << " particles" << std::endl;

    //////////////////////////////////////////////////////////////////
    int ii = 0; // index of beam particle candidate
    int iiloop = 0;
    double cut_v = 9999.; // cut value
    for (std::vector<const recob::PFParticle*> beam_slice : beam_slices) {
      const recob::PFParticle* particle = beam_slice.at(0);
      const TVector3 vtx = pfpUtil.GetPFParticleVertex(*particle,evt,fPFParticleTag,fTrackerTag);

      if (fVerbose)
        std::cout << "#Start position of beam particle NO." << iiloop <<
                     ": (" << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() <<
                     ")" << std::endl;
      
      // add selection
      double fom = abs(vtx.Z() - 30); // figure of merit (to be compared to the cut value)
      if (fom < cut_v) {
        cut_v = fom;
        ii = iiloop;
      }
      ++iiloop;
    }
    // Get the reconstructed PFParticle tagged as beam by Pandora
    const recob::PFParticle* particle 
        = ((fCheckSlicesForBeam && ii >= 0) ?
           beam_slices.at(ii).at(0) :
           beamParticles.at(0));

    //Get info from the PFP object identified as beam by pandora
    BeamPFPInfo(evt, particle, hitResults);
    //If MC, attempt to match to some MCParticle
    if( !evt.isRealData() ){
      BeamMatchInfo(evt, particle, true_beam_particle, clockData);
    }


    // Determine if the beam particle is track-like or shower-like
    const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*particle,evt,fPFParticleTag,fTrackerTag);
    const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*particle,evt,fPFParticleTag,fShowerTag);
    if( thisTrack ){
      //Get reconstructed beam track info
      BeamTrackInfo(evt, thisTrack, clockData);
      rightTrackID = thisTrack->ID();
    }
    else if( thisShower ){
      //Get reconstructed beam shower info 
      BeamShowerInfo(evt, thisShower);
    }

    //Get info from all PFPs associated as daughters to the primary particle
    DaughterPFPInfo(evt, particle, clockData, hitResults);

    //Gets info from forcing pandora to fit a track to the beam PFP regardless
    //of the BDT score
    BeamForcedTrackInfo(evt, particle);
    //To do: BeamForcedShowerInfo?

  }

  if (fCheckTruncation) CheckEff(evt, rightTrackID);
  //If MC, attempt to match to some MCParticle
  if( !evt.isRealData() ){
    TrueBeamInfo(evt, true_beam_particle, clockData, plist, trueToPFPs, hitResults);
  }

  //New geant4reweight stuff
  //To do: put in its own function
  /*
  if (!evt.isRealData() && fDoReweight) {
    if (fVerbose) std::cout << "Doing reweight" << std::endl;

    //Doing reweighting if the primary is a piplus
    if (true_beam_PDG == 211) {

      std::vector<G4ReweightTraj *> trajs = CreateNRWTrajs(
          *true_beam_particle, plist,
          fGeometryService, event);
      
      bool added = false;
      for (size_t i = 0; i < trajs.size(); ++i) {
        if (trajs[i]->GetNSteps() > 0) {
          for (size_t j = 0; j < ParSet.size(); ++j) {
            std::pair<double, double> pm_weights =
                MultiRW->GetPlusMinusSigmaParWeight((*trajs[i]), j);

            if (!added) {
              g4rw_alt_primary_plus_sigma_weight.push_back(pm_weights.first);
              g4rw_alt_primary_minus_sigma_weight.push_back(pm_weights.second);
            }
            else {
              g4rw_alt_primary_plus_sigma_weight[j] *= pm_weights.first;
              g4rw_alt_primary_minus_sigma_weight[j] *= pm_weights.second;
            }
          }
          added = true;
        }
      }

      //Weighting according to the full heirarchy
      std::vector<std::vector<G4ReweightTraj *>> new_full_created
          = BuildHierarchy(true_beam_ID, 211, plist, fGeometryService,
                           event, "LAr", false);
      if (fVerbose) {
        std::cout << "Created " << new_full_created.size() << " reweightable pi+"
                  << std::endl;
      }

      bool new_added = false;
      for (size_t i = 0; i < new_full_created.size(); ++i) {
        std::vector<G4ReweightTraj *> temp_trajs = new_full_created[i];
        if (fVerbose) std::cout << i << " n trajs: " << temp_trajs.size() << std::endl;
        for (size_t j = 0; j < temp_trajs.size(); ++j) {
          G4ReweightTraj * this_traj = temp_trajs[j];
          if (this_traj->GetNSteps() > 0) {
            for (size_t k = 0; k < ParSet.size(); ++k) {
              std::pair<double, double> pm_weights =
                  MultiRW->GetPlusMinusSigmaParWeight((*this_traj), k);

              if (!new_added) {
                g4rw_full_primary_plus_sigma_weight.push_back(pm_weights.first);
                g4rw_full_primary_minus_sigma_weight.push_back(pm_weights.second);
              }
              else {
                g4rw_full_primary_plus_sigma_weight[k] *= pm_weights.first;
                g4rw_full_primary_minus_sigma_weight[k] *= pm_weights.second;
              }
            }
            new_added = true;
          }
        }
      }

      G4RWGridWeights(new_full_created, ParSet, g4rw_full_grid_weights,
                      MultiRW);
      for (auto weights : g4rw_full_grid_weights) {
        g4rw_full_grid_coeffs.push_back(std::vector<double>());
        GetG4RWCoeffs(weights, g4rw_full_grid_coeffs.back());
      }
    }
  }*/

  //New style of weighting to get full hierarchy (i.e. from full event not 
  //just from the primary + downstreams) -- currently for piplus and proton
  if (!evt.isRealData() && fDoReweight) {
    std::vector<std::vector<G4ReweightTraj *>> piplus_hierarchy 
        = BuildHierarchy(true_beam_ID, 211, plist, *geom,
                         event, "LAr", false);

    //std::cout << "Primary hierarchy" << std::endl; 
    int ipart = 0;
    for (auto & traj_vec : piplus_hierarchy) {
      g4rw_piplus_traj_ps.push_back(std::vector<double>());
      g4rw_piplus_traj_lens.push_back(std::vector<double>());
      int itraj = 0;
      for (auto * traj : traj_vec) {
        //std::cout << "Part " << ipart << " Traj " << itraj << std::endl;
        for (size_t istep = 0; istep < traj->GetNSteps(); ++istep) {
          auto * step = traj->GetStep(istep);
          /*std::cout << "\t" << step->GetStepLength() << " " <<
                       step->GetFullPreStepP() <<
                       " " << " "  << std::endl;*/
          g4rw_piplus_traj_ps.back().push_back(step->GetFullPreStepP());
          g4rw_piplus_traj_lens.back().push_back(step->GetStepLength());
        }
        /*std::cout << "\t" << traj->HasChild(211).size() << " " <<
                     traj->HasChild(-211).size() << " " <<
                     traj->HasChild(111).size() << std::endl;*/
        ++itraj;
      }
      //Children should only be at the last traj
      g4rw_piplus_traj_npiplus.push_back(traj_vec.back()->HasChild(211).size());
      g4rw_piplus_traj_npiminus.push_back(traj_vec.back()->HasChild(-211).size());
      g4rw_piplus_traj_npi0.push_back(traj_vec.back()->HasChild(111).size());
      ++ipart;
    }

    G4RWGridWeights(piplus_hierarchy, ParSet, g4rw_full_grid_piplus_weights,
                    MultiRW);
    for (auto weights : g4rw_full_grid_piplus_weights) {
      g4rw_full_grid_piplus_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_full_grid_piplus_coeffs.back());


    }

    //Fine Weights
    G4RWGridWeights(piplus_hierarchy, FineParSet, g4rw_full_fine_piplus_weights,
                    FineMultiRW);
    for (auto weights : g4rw_full_fine_piplus_weights) {
      g4rw_full_fine_piplus_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_full_fine_piplus_coeffs.back());

    }


    G4RWGridWeights(piplus_hierarchy, FakeDataParSet,
                    g4rw_full_grid_piplus_weights_fake_data,
                    FakeDataMultiRW);
    for (auto weights : g4rw_full_grid_piplus_weights_fake_data) {
      g4rw_full_grid_piplus_coeffs_fake_data.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_full_grid_piplus_coeffs_fake_data.back());

    }

    //AbsCex Weights
    
    std::cout << "Getting abscex weigths" << std::endl;
    G4RWGridWeights(piplus_hierarchy, fAbsCexPars,
                    g4rw_full_grid_abscex_weights, fAbsCexMultiRW);
    for (auto weights : g4rw_full_grid_abscex_weights) {
      g4rw_full_grid_abscex_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_full_grid_abscex_coeffs.back());
    }

    //Proton
    std::vector<std::vector<G4ReweightTraj *>> proton_hierarchy 
        = BuildHierarchy(true_beam_ID, 2212, plist, *geom,
                         event, "LAr", false);
    G4RWGridWeights(proton_hierarchy, ProtParSet, g4rw_full_grid_proton_weights,
                    ProtMultiRW);
    for (auto weights : g4rw_full_grid_proton_weights) {
      g4rw_full_grid_proton_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_full_grid_proton_coeffs.back());

    }

    std::vector<std::vector<G4ReweightTraj *>> neutron_hierarchy 
        = BuildHierarchy(true_beam_ID, 2112, plist, *geom,
                         event, "LAr", false);
    G4RWGridWeights(neutron_hierarchy, ProtParSet,
                    g4rw_full_grid_neutron_weights, NeutronMultiRW);
    for (auto weights : g4rw_full_grid_neutron_weights) {
      g4rw_full_grid_neutron_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_full_grid_neutron_coeffs.back());

    }

    std::vector<std::vector<G4ReweightTraj *>> kplus_hierarchy 
        = BuildHierarchy(true_beam_ID, 321, plist, *geom,
                         event, "LAr", false);
    G4RWGridWeights(kplus_hierarchy, ProtParSet, g4rw_full_grid_kplus_weights,
                    KPlusMultiRW);
    for (auto weights : g4rw_full_grid_kplus_weights) {
      g4rw_full_grid_kplus_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_full_grid_kplus_coeffs.back());

    }

    std::vector<std::vector<G4ReweightTraj *>> downstream_piplus_hierarchy 
        = BuildHierarchy(true_beam_ID, 211, plist, *geom,
                         event, "LAr", true);
    G4RWGridWeights(downstream_piplus_hierarchy, ParSet,
                    g4rw_downstream_grid_piplus_weights, MultiRW);
    for (auto weights : g4rw_downstream_grid_piplus_weights) {
      g4rw_downstream_grid_piplus_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_downstream_grid_piplus_coeffs.back());
    }
    G4RWGridWeights(downstream_piplus_hierarchy, fAbsCexPars,
                    g4rw_downstream_grid_abscex_weights, fAbsCexMultiRW);
    for (auto weights : g4rw_downstream_grid_abscex_weights) {
      g4rw_downstream_grid_abscex_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_downstream_grid_abscex_coeffs.back());
    }


    std::cout << "Making trajs " << std::endl;
    std::vector<G4ReweightTraj *> trajs = CreateNRWTrajs(
        *true_beam_particle, plist,
        *geom, event);
    std::cout << trajs.size() << std::endl;

    std::vector<std::vector<G4ReweightTraj *>> temp_hierarchy = {trajs};
    G4MultiReweighter * primary_rw = 0x0;
    std::cout << true_beam_PDG << std::endl;
    if (true_beam_PDG == 211) {
      primary_rw = MultiRW;
    }
    else if (true_beam_PDG == 2212) {
      primary_rw = ProtMultiRW;
    }
    else if (true_beam_PDG == 321) {
      primary_rw = KPlusMultiRW;
    }
    else {
      primary_rw = MultiRW;
    }

    G4RWGridWeights(
        temp_hierarchy, ((true_beam_PDG == 211 || true_beam_PDG == -13) ?
                         ParSet : ProtParSet),
        g4rw_primary_grid_weights, primary_rw);
    std::cout << "Weights; " << g4rw_primary_grid_weights.size() << std::endl;
    for (auto weights : g4rw_primary_grid_weights) {
      g4rw_primary_grid_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_primary_grid_coeffs.back());
    }

    G4RWGridWeights(
        temp_hierarchy, fAbsCexPars,
        g4rw_primary_grid_abscex_weights, fAbsCexMultiRW);
    for (auto weights : g4rw_primary_grid_abscex_weights) {
      g4rw_primary_grid_abscex_coeffs.push_back(std::vector<double>());
      GetG4RWCoeffs(weights, g4rw_primary_grid_abscex_coeffs.back());
    }

    //FineWeights
    //

  }

  fTree->Fill();
}

void pduneana::PDSPAnalyzer::beginJob() {

  gROOT->SetBatch(1);

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("beamana","beam analysis tree");

  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("event", &event);
  fTree->Branch("MC", &MC);



  ///Reconstructed info
  fTree->Branch("reco_reconstructable_beam_event",&reco_reconstructable_beam_event);
  fTree->Branch("reco_beam_type", &reco_beam_type);
  fTree->Branch("reco_beam_startX", &reco_beam_startX);
  fTree->Branch("reco_beam_startY", &reco_beam_startY);
  fTree->Branch("reco_beam_startZ", &reco_beam_startZ);
  fTree->Branch("reco_beam_endX", &reco_beam_endX);
  fTree->Branch("reco_beam_endY", &reco_beam_endY);
  fTree->Branch("reco_beam_endZ", &reco_beam_endZ);
  fTree->Branch("true_beam_len", &true_beam_len);
  fTree->Branch("reco_beam_len", &reco_beam_len);
  fTree->Branch("test_branch", &test_branch);
  fTree->Branch("reco_beam_alt_len", &reco_beam_alt_len);
  fTree->Branch("for_truncation_method", &for_truncation_method);
  fTree->Branch("reco_beam_calo_startX", &reco_beam_calo_startX);
  fTree->Branch("reco_beam_calo_startY", &reco_beam_calo_startY);
  fTree->Branch("reco_beam_calo_startZ", &reco_beam_calo_startZ);
  fTree->Branch("reco_beam_calo_endX", &reco_beam_calo_endX);
  fTree->Branch("reco_beam_calo_endY", &reco_beam_calo_endY);
  fTree->Branch("reco_beam_calo_endZ", &reco_beam_calo_endZ);
  fTree->Branch("reco_beam_calo_startDirX", &reco_beam_calo_startDirX);
  fTree->Branch("reco_beam_calo_startDirY", &reco_beam_calo_startDirY);
  fTree->Branch("reco_beam_calo_startDirZ", &reco_beam_calo_startDirZ);
  fTree->Branch("reco_beam_calo_endDirX", &reco_beam_calo_endDirX);
  fTree->Branch("reco_beam_calo_endDirY", &reco_beam_calo_endDirY);
  fTree->Branch("reco_beam_calo_endDirZ", &reco_beam_calo_endDirZ);

  fTree->Branch("reco_beam_trackDirX", &reco_beam_trackDirX);
  fTree->Branch("reco_beam_trackDirY", &reco_beam_trackDirY);
  fTree->Branch("reco_beam_trackDirZ", &reco_beam_trackDirZ);
  fTree->Branch("reco_beam_trackEndDirX", &reco_beam_trackEndDirX);
  fTree->Branch("reco_beam_trackEndDirY", &reco_beam_trackEndDirY);
  fTree->Branch("reco_beam_trackEndDirZ", &reco_beam_trackEndDirZ);
  fTree->Branch("reco_beam_vertex_nHits", &reco_beam_vertex_nHits);
  fTree->Branch("reco_beam_vertex_michel_score", &reco_beam_vertex_michel_score);
  fTree->Branch("reco_beam_vertex_michel_score_weight_by_charge", &reco_beam_vertex_michel_score_weight_by_charge);
  fTree->Branch("reco_beam_trackID", &reco_beam_trackID);
  fTree->Branch("n_beam_slices", &n_beam_slices);
  fTree->Branch("n_beam_particles", &n_beam_particles);
  fTree->Branch("beam_track_IDs", &beam_track_IDs);
  fTree->Branch("beam_particle_scores", &beam_particle_scores);

  fTree->Branch("reco_beam_dQdX_SCE", &reco_beam_dQdX_SCE);
  fTree->Branch("reco_beam_EField_SCE", &reco_beam_EField_SCE);
  fTree->Branch("reco_beam_calo_X", &reco_beam_calo_X);
  fTree->Branch("reco_beam_calo_Y", &reco_beam_calo_Y);
  fTree->Branch("reco_beam_calo_Z", &reco_beam_calo_Z);
  fTree->Branch("reco_beam_dQ", &reco_beam_dQ);
  fTree->Branch("reco_beam_dEdX_SCE", &reco_beam_dEdX_SCE);
  fTree->Branch("reco_beam_calibrated_dEdX_SCE", &reco_beam_calibrated_dEdX_SCE);
  fTree->Branch("reco_beam_calibrated_dQdX_SCE", &reco_beam_calibrated_dQdX_SCE);
  fTree->Branch("reco_beam_resRange_SCE", &reco_beam_resRange_SCE);
  fTree->Branch("reco_beam_TrkPitch_SCE", &reco_beam_TrkPitch_SCE);

  if (fSaveHitIDEInfo) {
    fTree->Branch("reco_beam_hit_IDE_IDs", &reco_beam_hit_IDE_IDs);
    fTree->Branch("reco_beam_hit_IDE_electrons", &reco_beam_hit_IDE_electrons);
    fTree->Branch("reco_beam_hit_IDE_energies", &reco_beam_hit_IDE_energies);
    fTree->Branch("reco_beam_hit_IDE_origins", &reco_beam_hit_IDE_origins);
    fTree->Branch("reco_beam_hit_IDE_cosmic_electrons", &reco_beam_hit_IDE_cosmic_electrons);
    fTree->Branch("reco_beam_hit_IDE_beam_electrons", &reco_beam_hit_IDE_beam_electrons);
    fTree->Branch("reco_beam_hit_IDE_cosmic_energies", &reco_beam_hit_IDE_cosmic_energies);
    fTree->Branch("reco_beam_hit_IDE_beam_energies", &reco_beam_hit_IDE_beam_energies);
  }

  fTree->Branch("reco_beam_dQdX_NoSCE", &reco_beam_dQdX_NoSCE);
  fTree->Branch("reco_beam_dQ_NoSCE", &reco_beam_dQ_NoSCE);
  fTree->Branch("reco_beam_dEdX_NoSCE", &reco_beam_dEdX_NoSCE);
  fTree->Branch("reco_beam_calibrated_dEdX_NoSCE", &reco_beam_calibrated_dEdX_NoSCE);
  fTree->Branch("reco_beam_resRange_NoSCE", &reco_beam_resRange_NoSCE);
  fTree->Branch("reco_beam_TrkPitch_NoSCE", &reco_beam_TrkPitch_NoSCE);

  fTree->Branch("reco_beam_calo_wire", &reco_beam_calo_wire);
  fTree->Branch("reco_beam_calo_wire_z", &reco_beam_calo_wire_z);
  fTree->Branch("reco_beam_calo_wire_NoSCE", &reco_beam_calo_wire_NoSCE);
  fTree->Branch("reco_beam_calo_wire_z_NoSCE", &reco_beam_calo_wire_z_NoSCE);
  fTree->Branch("reco_beam_calo_tick", &reco_beam_calo_tick);
  fTree->Branch("reco_beam_calo_TPC", &reco_beam_calo_TPC);
  fTree->Branch("reco_beam_calo_TPC_NoSCE", &reco_beam_calo_TPC_NoSCE);

  fTree->Branch("reco_beam_flipped", &reco_beam_flipped);
  fTree->Branch("reco_beam_passes_beam_cuts", &reco_beam_passes_beam_cuts);

  fTree->Branch("reco_beam_PFP_ID", &reco_beam_PFP_ID);
  fTree->Branch("reco_beam_PFP_nHits", &reco_beam_PFP_nHits);
  fTree->Branch("reco_beam_PFP_trackScore", &reco_beam_PFP_trackScore);
  fTree->Branch("reco_beam_PFP_emScore", &reco_beam_PFP_emScore);
  fTree->Branch("reco_beam_PFP_michelScore", &reco_beam_PFP_michelScore);
  fTree->Branch("reco_beam_PFP_trackScore_collection", &reco_beam_PFP_trackScore_collection);
  fTree->Branch("reco_beam_PFP_emScore_collection", &reco_beam_PFP_emScore_collection);
  fTree->Branch("reco_beam_PFP_michelScore_collection", &reco_beam_PFP_michelScore_collection);
  fTree->Branch("reco_beam_PFP_trackScore_weight_by_charge", &reco_beam_PFP_trackScore_weight_by_charge);
  fTree->Branch("reco_beam_PFP_emScore_weight_by_charge", &reco_beam_PFP_emScore_weight_by_charge);
  fTree->Branch("reco_beam_PFP_michelScore_weight_by_charge", &reco_beam_PFP_michelScore_weight_by_charge);
  fTree->Branch("reco_beam_PFP_trackScore_collection_weight_by_charge", &reco_beam_PFP_trackScore_collection_weight_by_charge);
  fTree->Branch("reco_beam_PFP_emScore_collection_weight_by_charge", &reco_beam_PFP_emScore_collection_weight_by_charge);
  fTree->Branch("reco_beam_PFP_michelScore_collection_weight_by_charge", &reco_beam_PFP_michelScore_collection_weight_by_charge);

  if (fSaveRecoBeamAllTrack) {
    fTree->Branch("reco_beam_allTrack_ID",              &reco_beam_allTrack_ID);
    fTree->Branch("reco_beam_allTrack_beam_cuts",       &reco_beam_allTrack_beam_cuts);
    fTree->Branch("reco_beam_allTrack_flipped",         &reco_beam_allTrack_flipped);
    fTree->Branch("reco_beam_allTrack_len",             &reco_beam_allTrack_len);
    fTree->Branch("reco_beam_allTrack_startX",          &reco_beam_allTrack_startX);
    fTree->Branch("reco_beam_allTrack_startY",          &reco_beam_allTrack_startY);
    fTree->Branch("reco_beam_allTrack_startZ",          &reco_beam_allTrack_startZ);
    fTree->Branch("reco_beam_allTrack_endX",            &reco_beam_allTrack_endX);
    fTree->Branch("reco_beam_allTrack_endY",            &reco_beam_allTrack_endY);
    fTree->Branch("reco_beam_allTrack_endZ",            &reco_beam_allTrack_endZ);
    fTree->Branch("reco_beam_allTrack_trackDirX",       &reco_beam_allTrack_trackDirX);
    fTree->Branch("reco_beam_allTrack_trackDirY",       &reco_beam_allTrack_trackDirY);
    fTree->Branch("reco_beam_allTrack_trackDirZ",       &reco_beam_allTrack_trackDirZ);
    fTree->Branch("reco_beam_allTrack_trackEndDirX",    &reco_beam_allTrack_trackEndDirX);
    fTree->Branch("reco_beam_allTrack_trackEndDirY",    &reco_beam_allTrack_trackEndDirY);
    fTree->Branch("reco_beam_allTrack_trackEndDirZ",    &reco_beam_allTrack_trackEndDirZ);
    fTree->Branch("reco_beam_allTrack_resRange",        &reco_beam_allTrack_resRange);
    fTree->Branch("reco_beam_allTrack_calibrated_dEdX", &reco_beam_allTrack_calibrated_dEdX);
    fTree->Branch("reco_beam_allTrack_Chi2_proton",     &reco_beam_allTrack_Chi2_proton);
    fTree->Branch("reco_beam_allTrack_Chi2_ndof",       &reco_beam_allTrack_Chi2_ndof);
    fTree->Branch("reco_beam_alt_len_allTrack", &reco_beam_alt_len_allTrack);
    fTree->Branch("reco_beam_calo_startX_allTrack", &reco_beam_calo_startX_allTrack);
    fTree->Branch("reco_beam_calo_startY_allTrack", &reco_beam_calo_startY_allTrack);
    fTree->Branch("reco_beam_calo_startZ_allTrack", &reco_beam_calo_startZ_allTrack);
    fTree->Branch("reco_beam_calo_endX_allTrack", &reco_beam_calo_endX_allTrack);
    fTree->Branch("reco_beam_calo_endY_allTrack", &reco_beam_calo_endY_allTrack);
    fTree->Branch("reco_beam_calo_endZ_allTrack", &reco_beam_calo_endZ_allTrack);
    fTree->Branch("reco_beam_interactingEnergy_allTrack", &reco_beam_interactingEnergy_allTrack);
    fTree->Branch("reco_beam_calo_wire_allTrack", &reco_beam_calo_wire_allTrack);
    fTree->Branch("reco_beam_calo_X_allTrack", &reco_beam_calo_X_allTrack);
    fTree->Branch("reco_beam_calo_Y_allTrack", &reco_beam_calo_Y_allTrack);
    fTree->Branch("reco_beam_calo_Z_allTrack", &reco_beam_calo_Z_allTrack);
    fTree->Branch("reco_beam_TrkPitch_SCE_allTrack", &reco_beam_TrkPitch_SCE_allTrack);
    fTree->Branch("reco_beam_vertex_michel_score_weight_by_charge_allTrack", &reco_beam_vertex_michel_score_weight_by_charge_allTrack);
    fTree->Branch("reco_beam_vertex_nHits_allTrack", &reco_beam_vertex_nHits_allTrack);
    fTree->Branch("reco_beam_vertex_michel_score_allTrack", &reco_beam_vertex_michel_score_allTrack);
  }

  fTree->Branch("reco_track_startX", &reco_track_startX);
  fTree->Branch("reco_track_startY", &reco_track_startY);
  fTree->Branch("reco_track_startZ", &reco_track_startZ);
  fTree->Branch("reco_track_endX", &reco_track_endX);
  fTree->Branch("reco_track_endY", &reco_track_endY);
  fTree->Branch("reco_track_endZ", &reco_track_endZ);
  fTree->Branch("reco_track_michel_score", &reco_track_michel_score);
  fTree->Branch("reco_track_michel_score_weight_by_charge", &reco_track_michel_score_weight_by_charge);
  fTree->Branch("reco_track_ID", &reco_track_ID);
  fTree->Branch("reco_track_nHits", &reco_track_nHits);

  //Alternative reco
  fTree->Branch("reco_daughter_PFP_true_byHits_PDG", &reco_daughter_PFP_true_byHits_PDG);
  fTree->Branch("reco_daughter_PFP_true_byHits_ID", &reco_daughter_PFP_true_byHits_ID);
  fTree->Branch("reco_daughter_PFP_true_byHits_origin", &reco_daughter_PFP_true_byHits_origin);
  fTree->Branch("reco_daughter_PFP_true_byHits_parID", &reco_daughter_PFP_true_byHits_parID);
  fTree->Branch("reco_daughter_PFP_true_byHits_parPDG", &reco_daughter_PFP_true_byHits_parPDG);
  fTree->Branch("reco_daughter_PFP_true_byHits_process", &reco_daughter_PFP_true_byHits_process);
  fTree->Branch("reco_daughter_PFP_true_byHits_sharedHits", &reco_daughter_PFP_true_byHits_sharedHits);
  fTree->Branch("reco_daughter_PFP_true_byHits_emHits", &reco_daughter_PFP_true_byHits_emHits);

  fTree->Branch("reco_daughter_PFP_true_byHits_len", &reco_daughter_PFP_true_byHits_len);
  fTree->Branch("reco_daughter_PFP_true_byHits_startX", &reco_daughter_PFP_true_byHits_startX);
  fTree->Branch("reco_daughter_PFP_true_byHits_startY", &reco_daughter_PFP_true_byHits_startY);
  fTree->Branch("reco_daughter_PFP_true_byHits_startZ", &reco_daughter_PFP_true_byHits_startZ);
  fTree->Branch("reco_daughter_PFP_true_byHits_endX", &reco_daughter_PFP_true_byHits_endX);
  fTree->Branch("reco_daughter_PFP_true_byHits_endY", &reco_daughter_PFP_true_byHits_endY);
  fTree->Branch("reco_daughter_PFP_true_byHits_endZ", &reco_daughter_PFP_true_byHits_endZ);

  fTree->Branch("reco_daughter_PFP_true_byHits_startPx", &reco_daughter_PFP_true_byHits_startPx);
  fTree->Branch("reco_daughter_PFP_true_byHits_startPy", &reco_daughter_PFP_true_byHits_startPy);
  fTree->Branch("reco_daughter_PFP_true_byHits_startPz", &reco_daughter_PFP_true_byHits_startPz);
  fTree->Branch("reco_daughter_PFP_true_byHits_startP", &reco_daughter_PFP_true_byHits_startP);
  fTree->Branch("reco_daughter_PFP_true_byHits_startE", &reco_daughter_PFP_true_byHits_startE);
  fTree->Branch("reco_daughter_PFP_true_byHits_endProcess", &reco_daughter_PFP_true_byHits_endProcess);
  fTree->Branch("reco_daughter_PFP_true_byHits_purity", &reco_daughter_PFP_true_byHits_purity);
  fTree->Branch("reco_daughter_PFP_true_byHits_completeness", &reco_daughter_PFP_true_byHits_completeness);
  fTree->Branch("reco_daughter_PFP_true_byE_PDG", &reco_daughter_PFP_true_byE_PDG);
  fTree->Branch("reco_daughter_PFP_true_byE_len", &reco_daughter_PFP_true_byE_len);
  fTree->Branch("reco_daughter_PFP_true_byE_completeness", &reco_daughter_PFP_true_byE_completeness);
  fTree->Branch("reco_daughter_PFP_true_byE_purity", &reco_daughter_PFP_true_byE_purity);

  fTree->Branch("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID);
  fTree->Branch("reco_daughter_allTrack_dQdX_SCE", &reco_daughter_allTrack_dQdX_SCE);
  fTree->Branch("reco_daughter_allTrack_calibrated_dQdX_SCE", &reco_daughter_allTrack_calibrated_dQdX_SCE);
  fTree->Branch("reco_daughter_allTrack_EField_SCE", &reco_daughter_allTrack_EField_SCE);
  fTree->Branch("reco_daughter_allTrack_dEdX_SCE", &reco_daughter_allTrack_dEdX_SCE);
  fTree->Branch("reco_daughter_allTrack_resRange_SCE", &reco_daughter_allTrack_resRange_SCE);
  fTree->Branch("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE);

  fTree->Branch("reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton);
  fTree->Branch("reco_daughter_allTrack_Chi2_pion", &reco_daughter_allTrack_Chi2_pion);
  fTree->Branch("reco_daughter_allTrack_Chi2_muon", &reco_daughter_allTrack_Chi2_muon);
  fTree->Branch("reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof);
  fTree->Branch("reco_daughter_allTrack_Chi2_ndof_pion", &reco_daughter_allTrack_Chi2_ndof_pion);
  fTree->Branch("reco_daughter_allTrack_Chi2_ndof_muon", &reco_daughter_allTrack_Chi2_ndof_muon);


  ///Calorimetry/chi2 planes 0 and 1
  fTree->Branch("reco_daughter_allTrack_Chi2_proton_plane0",
                &reco_daughter_allTrack_Chi2_proton_plane0);
  fTree->Branch("reco_daughter_allTrack_Chi2_proton_plane1",
                &reco_daughter_allTrack_Chi2_proton_plane1);

  fTree->Branch("reco_daughter_allTrack_Chi2_ndof_plane0",
                &reco_daughter_allTrack_Chi2_ndof_plane0);
  fTree->Branch("reco_daughter_allTrack_Chi2_ndof_plane1",
                &reco_daughter_allTrack_Chi2_ndof_plane1);

  fTree->Branch("reco_daughter_allTrack_calibrated_dEdX_SCE_plane0",
                &reco_daughter_allTrack_calibrated_dEdX_SCE_plane0);
  fTree->Branch("reco_daughter_allTrack_calibrated_dEdX_SCE_plane1",
                &reco_daughter_allTrack_calibrated_dEdX_SCE_plane1);
  fTree->Branch("reco_daughter_allTrack_resRange_plane0",
                &reco_daughter_allTrack_resRange_plane0);
  fTree->Branch("reco_daughter_allTrack_resRange_plane1",
                &reco_daughter_allTrack_resRange_plane1);
  ///////////////////////////////////

  fTree->Branch("reco_daughter_allTrack_Theta", &reco_daughter_allTrack_Theta);
  fTree->Branch("reco_daughter_allTrack_Phi", &reco_daughter_allTrack_Phi);
  fTree->Branch("reco_daughter_allTrack_startDirX", &reco_daughter_allTrack_startDirX);
  fTree->Branch("reco_daughter_allTrack_startDirY", &reco_daughter_allTrack_startDirY);
  fTree->Branch("reco_daughter_allTrack_startDirZ", &reco_daughter_allTrack_startDirZ);

  fTree->Branch("reco_daughter_allTrack_len", &reco_daughter_allTrack_len);
  fTree->Branch("reco_daughter_allTrack_alt_len", &reco_daughter_allTrack_alt_len);
  fTree->Branch("reco_daughter_allTrack_startX", &reco_daughter_allTrack_startX);
  fTree->Branch("reco_daughter_allTrack_startY", &reco_daughter_allTrack_startY);
  fTree->Branch("reco_daughter_allTrack_startZ", &reco_daughter_allTrack_startZ);
  fTree->Branch("reco_daughter_allTrack_endX", &reco_daughter_allTrack_endX);
  fTree->Branch("reco_daughter_allTrack_endY", &reco_daughter_allTrack_endY);
  fTree->Branch("reco_daughter_allTrack_endZ", &reco_daughter_allTrack_endZ);
  fTree->Branch("reco_daughter_allTrack_calo_X", &reco_daughter_allTrack_calo_X);
  fTree->Branch("reco_daughter_allTrack_calo_Y", &reco_daughter_allTrack_calo_Y);
  fTree->Branch("reco_daughter_allTrack_calo_Z", &reco_daughter_allTrack_calo_Z);

  fTree->Branch("reco_daughter_allTrack_vertex_michel_score",
                &reco_daughter_allTrack_vertex_michel_score);
  fTree->Branch("reco_daughter_allTrack_vertex_nHits",
                &reco_daughter_allTrack_vertex_nHits);
  //////

  fTree->Branch("reco_daughter_allShower_ID", &reco_daughter_allShower_ID);
  fTree->Branch("reco_daughter_allShower_len", &reco_daughter_allShower_len);
  fTree->Branch("reco_daughter_allShower_startX", &reco_daughter_allShower_startX);
  fTree->Branch("reco_daughter_allShower_startY", &reco_daughter_allShower_startY);
  fTree->Branch("reco_daughter_allShower_startZ", &reco_daughter_allShower_startZ);
  fTree->Branch("reco_daughter_allShower_dirX", &reco_daughter_allShower_dirX);
  fTree->Branch("reco_daughter_allShower_dirY", &reco_daughter_allShower_dirY);
  fTree->Branch("reco_daughter_allShower_dirZ", &reco_daughter_allShower_dirZ);
  fTree->Branch("reco_daughter_allShower_energy", &reco_daughter_allShower_energy);
  fTree->Branch("reco_daughter_allShower_calibrated_energy", &reco_daughter_allShower_calibrated_energy);



  fTree->Branch("reco_daughter_PFP_ID", &reco_daughter_PFP_ID);
  fTree->Branch("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);
  fTree->Branch("reco_daughter_PFP_nHits_collection",
                &reco_daughter_PFP_nHits_collection);
  fTree->Branch("reco_daughter_PFP_trackScore", &reco_daughter_PFP_trackScore);
  fTree->Branch("reco_daughter_PFP_emScore", &reco_daughter_PFP_emScore);
  fTree->Branch("reco_daughter_PFP_michelScore", &reco_daughter_PFP_michelScore);
  fTree->Branch("reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection);
  fTree->Branch("reco_daughter_PFP_emScore_collection", &reco_daughter_PFP_emScore_collection);
  fTree->Branch("reco_daughter_PFP_michelScore_collection", &reco_daughter_PFP_michelScore_collection);

  fTree->Branch("reco_daughter_PFP_trackScore_weight_by_charge", &reco_daughter_PFP_trackScore_weight_by_charge);
  fTree->Branch("reco_daughter_PFP_emScore_weight_by_charge", &reco_daughter_PFP_emScore_weight_by_charge);
  fTree->Branch("reco_daughter_PFP_michelScore_weight_by_charge", &reco_daughter_PFP_michelScore_weight_by_charge);
  fTree->Branch("reco_daughter_PFP_trackScore_collection_weight_by_charge", &reco_daughter_PFP_trackScore_collection_weight_by_charge);
  fTree->Branch("reco_daughter_PFP_emScore_collection_weight_by_charge", &reco_daughter_PFP_emScore_collection_weight_by_charge);
  fTree->Branch("reco_daughter_PFP_michelScore_collection_weight_by_charge", &reco_daughter_PFP_michelScore_collection_weight_by_charge);

  fTree->Branch("reco_daughter_pandora_type", &reco_daughter_pandora_type);

  fTree->Branch("true_beam_PDG", &true_beam_PDG);
  fTree->Branch("true_beam_mass", &true_beam_mass);
  fTree->Branch("true_beam_ID", &true_beam_ID);
  fTree->Branch("reco_beam_true_byHits_in_primary_hierarchy",
                &reco_beam_true_byHits_in_primary_hierarchy);
  fTree->Branch("true_beam_endProcess", &true_beam_endProcess);
  fTree->Branch("true_beam_endX", &true_beam_endX);
  fTree->Branch("true_beam_endY", &true_beam_endY);
  fTree->Branch("true_beam_endZ", &true_beam_endZ);
  fTree->Branch("true_beam_endX_SCE", &true_beam_endX_SCE);
  fTree->Branch("true_beam_endY_SCE", &true_beam_endY_SCE);
  fTree->Branch("true_beam_endZ_SCE", &true_beam_endZ_SCE);
  fTree->Branch("true_beam_startX", &true_beam_startX);
  fTree->Branch("true_beam_startY", &true_beam_startY);
  fTree->Branch("true_beam_startZ", &true_beam_startZ);

  fTree->Branch("true_beam_startPx", &true_beam_startPx);
  fTree->Branch("true_beam_startPy", &true_beam_startPy);
  fTree->Branch("true_beam_startPz", &true_beam_startPz);
  fTree->Branch("true_beam_startP", &true_beam_startP);

  fTree->Branch("true_beam_endPx", &true_beam_endPx);
  fTree->Branch("true_beam_endPy", &true_beam_endPy);
  fTree->Branch("true_beam_endPz", &true_beam_endPz);
  fTree->Branch("true_beam_endP", &true_beam_endP);
  fTree->Branch("true_beam_endP2", &true_beam_endP2);
  fTree->Branch("true_beam_last_len", &true_beam_last_len);

  fTree->Branch("true_beam_startDirX", &true_beam_startDirX);
  fTree->Branch("true_beam_startDirY", &true_beam_startDirY);
  fTree->Branch("true_beam_startDirZ", &true_beam_startDirZ);

  fTree->Branch("true_beam_nElasticScatters", &true_beam_nElasticScatters);
  fTree->Branch("true_beam_elastic_costheta", &true_beam_elastic_costheta);
  fTree->Branch("true_beam_elastic_X", &true_beam_elastic_X);
  fTree->Branch("true_beam_elastic_Y", &true_beam_elastic_Y);
  fTree->Branch("true_beam_elastic_Z", &true_beam_elastic_Z);
  fTree->Branch("true_beam_elastic_deltaE", &true_beam_elastic_deltaE);
  fTree->Branch("true_beam_elastic_IDE_edep", &true_beam_elastic_IDE_edep);
  fTree->Branch("true_beam_IDE_totalDep",    &true_beam_IDE_totalDep);

  fTree->Branch("true_beam_nHits", &true_beam_nHits);
  fTree->Branch("true_beam_reco_byHits_PFP_ID", &true_beam_reco_byHits_PFP_ID);
  fTree->Branch("true_beam_reco_byHits_PFP_nHits", &true_beam_reco_byHits_PFP_nHits);
  fTree->Branch("true_beam_reco_byHits_allTrack_ID", &true_beam_reco_byHits_allTrack_ID);

  fTree->Branch("true_daughter_nPi0", &true_daughter_nPi0);
  fTree->Branch("true_daughter_nPiPlus", &true_daughter_nPiPlus);
  fTree->Branch("true_daughter_nProton", &true_daughter_nProton);
  fTree->Branch("true_daughter_nNeutron", &true_daughter_nNeutron);
  fTree->Branch("true_daughter_nPiMinus", &true_daughter_nPiMinus);
  fTree->Branch("true_daughter_nNucleus", &true_daughter_nNucleus);

  fTree->Branch("reco_beam_vertex_slice", &reco_beam_vertex_slice);

  fTree->Branch("true_beam_daughter_PDG", &true_beam_daughter_PDG);
  fTree->Branch("true_beam_daughter_ID", &true_beam_daughter_ID);
  fTree->Branch("true_beam_daughter_len", &true_beam_daughter_len);
  fTree->Branch("true_beam_daughter_startX", &true_beam_daughter_startX);
  fTree->Branch("true_beam_daughter_startY", &true_beam_daughter_startY);
  fTree->Branch("true_beam_daughter_startZ", &true_beam_daughter_startZ);
  fTree->Branch("true_beam_daughter_startPx", &true_beam_daughter_startPx);
  fTree->Branch("true_beam_daughter_startPy", &true_beam_daughter_startPy);
  fTree->Branch("true_beam_daughter_startPz", &true_beam_daughter_startPz);
  fTree->Branch("true_beam_daughter_startP", &true_beam_daughter_startP);
  fTree->Branch("true_beam_daughter_endX", &true_beam_daughter_endX);
  fTree->Branch("true_beam_daughter_endY", &true_beam_daughter_endY);
  fTree->Branch("true_beam_daughter_endZ", &true_beam_daughter_endZ);
  fTree->Branch("true_beam_daughter_Process", &true_beam_daughter_Process);
  fTree->Branch("true_beam_daughter_endProcess", &true_beam_daughter_endProcess);
  fTree->Branch("true_beam_daughter_nHits", &true_beam_daughter_nHits);

  fTree->Branch("true_beam_daughter_reco_byHits_PFP_ID", &true_beam_daughter_reco_byHits_PFP_ID);
  fTree->Branch("true_beam_daughter_reco_byHits_PFP_nHits", &true_beam_daughter_reco_byHits_PFP_nHits);
  fTree->Branch("true_beam_daughter_reco_byHits_PFP_trackScore", &true_beam_daughter_reco_byHits_PFP_trackScore);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_ID", &true_beam_daughter_reco_byHits_allTrack_ID);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_startX", &true_beam_daughter_reco_byHits_allTrack_startX);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_startY", &true_beam_daughter_reco_byHits_allTrack_startY);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_startZ", &true_beam_daughter_reco_byHits_allTrack_startZ);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_endX", &true_beam_daughter_reco_byHits_allTrack_endX);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_endY", &true_beam_daughter_reco_byHits_allTrack_endY);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_endZ", &true_beam_daughter_reco_byHits_allTrack_endZ);
  fTree->Branch("true_beam_daughter_reco_byHits_allTrack_len", &true_beam_daughter_reco_byHits_allTrack_len);

  fTree->Branch("true_beam_daughter_reco_byHits_allShower_ID", &true_beam_daughter_reco_byHits_allShower_ID);
  fTree->Branch("true_beam_daughter_reco_byHits_allShower_startX", &true_beam_daughter_reco_byHits_allShower_startX);
  fTree->Branch("true_beam_daughter_reco_byHits_allShower_startY", &true_beam_daughter_reco_byHits_allShower_startY);
  fTree->Branch("true_beam_daughter_reco_byHits_allShower_startZ", &true_beam_daughter_reco_byHits_allShower_startZ);
  fTree->Branch("true_beam_daughter_reco_byHits_allShower_len", &true_beam_daughter_reco_byHits_allShower_len);

  fTree->Branch("true_beam_Pi0_decay_ID", &true_beam_Pi0_decay_ID);
  fTree->Branch("true_beam_Pi0_decay_parID", &true_beam_Pi0_decay_parID);
  fTree->Branch("true_beam_Pi0_decay_PDG", &true_beam_Pi0_decay_PDG);
  fTree->Branch("true_beam_Pi0_decay_startP", &true_beam_Pi0_decay_startP);
  fTree->Branch("true_beam_Pi0_decay_startPx", &true_beam_Pi0_decay_startPx);
  fTree->Branch("true_beam_Pi0_decay_startPy", &true_beam_Pi0_decay_startPy);
  fTree->Branch("true_beam_Pi0_decay_startPz", &true_beam_Pi0_decay_startPz);
  fTree->Branch("true_beam_Pi0_decay_startX", &true_beam_Pi0_decay_startX);
  fTree->Branch("true_beam_Pi0_decay_startY", &true_beam_Pi0_decay_startY);
  fTree->Branch("true_beam_Pi0_decay_startZ", &true_beam_Pi0_decay_startZ);

  fTree->Branch("true_beam_Pi0_decay_len", &true_beam_Pi0_decay_len);
  fTree->Branch("true_beam_Pi0_decay_nHits", &true_beam_Pi0_decay_nHits);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_PFP_ID", &true_beam_Pi0_decay_reco_byHits_PFP_ID);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_PFP_nHits", &true_beam_Pi0_decay_reco_byHits_PFP_nHits);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_PFP_trackScore", &true_beam_Pi0_decay_reco_byHits_PFP_trackScore);

  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_ID", &true_beam_Pi0_decay_reco_byHits_allTrack_ID);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_startX", &true_beam_Pi0_decay_reco_byHits_allTrack_startX);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_startY", &true_beam_Pi0_decay_reco_byHits_allTrack_startY);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_startZ", &true_beam_Pi0_decay_reco_byHits_allTrack_startZ);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_endX", &true_beam_Pi0_decay_reco_byHits_allTrack_endX);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_endY", &true_beam_Pi0_decay_reco_byHits_allTrack_endY);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_endZ", &true_beam_Pi0_decay_reco_byHits_allTrack_endZ);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allTrack_len", &true_beam_Pi0_decay_reco_byHits_allTrack_len);

  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_ID", &true_beam_Pi0_decay_reco_byHits_allShower_ID);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_startX", &true_beam_Pi0_decay_reco_byHits_allShower_startX);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_startY", &true_beam_Pi0_decay_reco_byHits_allShower_startY);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_startZ", &true_beam_Pi0_decay_reco_byHits_allShower_startZ);
  fTree->Branch("true_beam_Pi0_decay_reco_byHits_allShower_len", &true_beam_Pi0_decay_reco_byHits_allShower_len);

  fTree->Branch("true_beam_grand_daughter_ID", &true_beam_grand_daughter_ID);
  fTree->Branch("true_beam_grand_daughter_parID", &true_beam_grand_daughter_parID);
  fTree->Branch("true_beam_grand_daughter_PDG", &true_beam_grand_daughter_PDG);
  fTree->Branch("true_beam_grand_daughter_nHits", &true_beam_grand_daughter_nHits);
  fTree->Branch("true_beam_grand_daughter_Process", &true_beam_grand_daughter_Process);
  fTree->Branch("true_beam_grand_daughter_endProcess", &true_beam_grand_daughter_endProcess);

  ////Matching reco to truth
  fTree->Branch("reco_beam_true_byE_endProcess", &reco_beam_true_byE_endProcess);
  fTree->Branch("reco_beam_true_byE_process", &reco_beam_true_byE_process);
  fTree->Branch("reco_beam_true_byE_origin", &reco_beam_true_byE_origin);
  fTree->Branch("reco_beam_true_byE_PDG", &reco_beam_true_byE_PDG);
  fTree->Branch("reco_beam_true_byE_ID", &reco_beam_true_byE_ID);

  fTree->Branch("reco_beam_true_byHits_endProcess", &reco_beam_true_byHits_endProcess);
  fTree->Branch("reco_beam_true_byHits_process", &reco_beam_true_byHits_process);
  fTree->Branch("reco_beam_true_byHits_origin", &reco_beam_true_byHits_origin);
  fTree->Branch("reco_beam_true_byHits_PDG", &reco_beam_true_byHits_PDG);
  fTree->Branch("reco_beam_true_byHits_ID", &reco_beam_true_byHits_ID);

  fTree->Branch("reco_beam_true_byE_matched", &reco_beam_true_byE_matched);
  fTree->Branch("reco_beam_true_byHits_matched", &reco_beam_true_byHits_matched);
  fTree->Branch("reco_beam_true_byHits_purity", &reco_beam_true_byHits_purity);

  fTree->Branch("true_beam_processes", &true_beam_processes);

  fTree->Branch("beam_inst_P", &beam_inst_P);
  fTree->Branch("beam_inst_C0", &beam_inst_C0);
  fTree->Branch("beam_inst_C1", &beam_inst_C1);
  fTree->Branch("beam_inst_C0_pressure", &beam_inst_C0_pressure);
  fTree->Branch("beam_inst_C1_pressure", &beam_inst_C1_pressure);
  fTree->Branch("beam_inst_TOF", &beam_inst_TOF);
  fTree->Branch("beam_inst_TOF_Chan", &beam_inst_TOF_Chan);
  fTree->Branch("beam_inst_X", &beam_inst_X);
  fTree->Branch("beam_inst_Y", &beam_inst_Y);
  fTree->Branch("beam_inst_Z", &beam_inst_Z);
  fTree->Branch("beam_inst_dirX", &beam_inst_dirX);
  fTree->Branch("beam_inst_dirY", &beam_inst_dirY);
  fTree->Branch("beam_inst_dirZ", &beam_inst_dirZ);

  fTree->Branch("beam_inst_nFibersP1", &beam_inst_nFibersP1);
  fTree->Branch("beam_inst_nFibersP2", &beam_inst_nFibersP2);
  fTree->Branch("beam_inst_nFibersP3", &beam_inst_nFibersP3);
  fTree->Branch("beam_inst_PDG_candidates", &beam_inst_PDG_candidates);
  fTree->Branch("beam_inst_nTracks", &beam_inst_nTracks);
  fTree->Branch("beam_inst_nMomenta", &beam_inst_nMomenta);
  fTree->Branch("beam_inst_valid", &beam_inst_valid);
  fTree->Branch("beam_inst_trigger", &beam_inst_trigger);


  fTree->Branch("reco_beam_Chi2_proton", &reco_beam_Chi2_proton);
  fTree->Branch("reco_beam_Chi2_muon", &reco_beam_Chi2_muon);
  fTree->Branch("reco_beam_Chi2_ndof", &reco_beam_Chi2_ndof);


  fTree->Branch("reco_daughter_allTrack_momByRange_proton", &reco_daughter_allTrack_momByRange_proton);
  fTree->Branch("reco_daughter_allTrack_momByRange_muon", &reco_daughter_allTrack_momByRange_muon);
  fTree->Branch("reco_beam_momByRange_proton", &reco_beam_momByRange_proton);
  fTree->Branch("reco_beam_momByRange_muon", &reco_beam_momByRange_muon);

  fTree->Branch("reco_daughter_allTrack_momByRange_alt_proton", &reco_daughter_allTrack_momByRange_alt_proton);
  fTree->Branch("reco_daughter_allTrack_momByRange_alt_muon", &reco_daughter_allTrack_momByRange_alt_muon);
  fTree->Branch("reco_beam_momByRange_alt_proton", &reco_beam_momByRange_alt_proton);
  fTree->Branch("reco_beam_momByRange_alt_muon", &reco_beam_momByRange_alt_muon);

  fTree->Branch("reco_beam_true_byE_endPx", &reco_beam_true_byE_endPx);
  fTree->Branch("reco_beam_true_byE_endPy", &reco_beam_true_byE_endPy);
  fTree->Branch("reco_beam_true_byE_endPz", &reco_beam_true_byE_endPz);
  fTree->Branch("reco_beam_true_byE_endE", &reco_beam_true_byE_endE);
  fTree->Branch("reco_beam_true_byE_endP", &reco_beam_true_byE_endP);

  fTree->Branch("reco_beam_true_byE_startPx", &reco_beam_true_byE_startPx);
  fTree->Branch("reco_beam_true_byE_startPy", &reco_beam_true_byE_startPy);
  fTree->Branch("reco_beam_true_byE_startPz", &reco_beam_true_byE_startPz);
  fTree->Branch("reco_beam_true_byE_startE", &reco_beam_true_byE_startE);
  fTree->Branch("reco_beam_true_byE_startP", &reco_beam_true_byE_startP);


  fTree->Branch("reco_beam_true_byHits_endPx", &reco_beam_true_byHits_endPx);
  fTree->Branch("reco_beam_true_byHits_endPy", &reco_beam_true_byHits_endPy);
  fTree->Branch("reco_beam_true_byHits_endPz", &reco_beam_true_byHits_endPz);
  fTree->Branch("reco_beam_true_byHits_endE", &reco_beam_true_byHits_endE);
  fTree->Branch("reco_beam_true_byHits_endP", &reco_beam_true_byHits_endP);

  fTree->Branch("reco_beam_true_byHits_startPx", &reco_beam_true_byHits_startPx);
  fTree->Branch("reco_beam_true_byHits_startPy", &reco_beam_true_byHits_startPy);
  fTree->Branch("reco_beam_true_byHits_startPz", &reco_beam_true_byHits_startPz);
  fTree->Branch("reco_beam_true_byHits_startE", &reco_beam_true_byHits_startE);
  fTree->Branch("reco_beam_true_byHits_startP", &reco_beam_true_byHits_startP);

  fTree->Branch("reco_beam_incidentEnergies", &reco_beam_incidentEnergies);
  fTree->Branch("reco_beam_interactingEnergy", &reco_beam_interactingEnergy);
  //fTree->Branch("reco_beam_incidentEnergies_allTrack", &reco_beam_incidentEnergies_allTrack);
  fTree->Branch("true_beam_incidentEnergies", &true_beam_incidentEnergies);
  fTree->Branch("true_beam_interactingEnergy", &true_beam_interactingEnergy);
  fTree->Branch("true_beam_slices", &true_beam_slices);
  fTree->Branch("true_beam_slices_found", &true_beam_slices_found);
  fTree->Branch("true_beam_slices_deltaE", &true_beam_slices_deltaE);
  fTree->Branch("em_energy", &em_energy);
  fTree->Branch("true_beam_traj_X", &true_beam_traj_X);
  fTree->Branch("true_beam_traj_Y", &true_beam_traj_Y);
  fTree->Branch("true_beam_traj_Z", &true_beam_traj_Z);
  fTree->Branch("true_beam_traj_Px", &true_beam_traj_Px);
  fTree->Branch("true_beam_traj_Py", &true_beam_traj_Py);
  fTree->Branch("true_beam_traj_Pz", &true_beam_traj_Pz);
  fTree->Branch("true_beam_traj_KE", &true_beam_traj_KE);
  fTree->Branch("true_beam_traj_X_SCE", &true_beam_traj_X_SCE);
  fTree->Branch("true_beam_traj_Y_SCE", &true_beam_traj_Y_SCE);
  fTree->Branch("true_beam_traj_Z_SCE", &true_beam_traj_Z_SCE);
  fTree->Branch("true_beam_is_scraper", &true_beam_is_scraper);

  //fTree->Branch("g4rw_primary_plus_sigma_weight", &g4rw_primary_plus_sigma_weight);
  //fTree->Branch("g4rw_primary_minus_sigma_weight", &g4rw_primary_minus_sigma_weight);
  //fTree->Branch("g4rw_primary_var", &g4rw_primary_var);

  /*
  fTree->Branch("g4rw_alt_primary_plus_sigma_weight",
                &g4rw_alt_primary_plus_sigma_weight);
  fTree->Branch("g4rw_alt_primary_minus_sigma_weight",
                &g4rw_alt_primary_minus_sigma_weight);

  fTree->Branch("g4rw_full_primary_plus_sigma_weight",
                &g4rw_full_primary_plus_sigma_weight);
  fTree->Branch("g4rw_full_primary_minus_sigma_weight",
                &g4rw_full_primary_minus_sigma_weight);*/
  if (fSaveG4RWWeights) {
    //fTree->Branch("g4rw_full_grid_weights", &g4rw_full_grid_weights);
    fTree->Branch("g4rw_full_grid_piplus_weights",  &g4rw_full_grid_piplus_weights);
    fTree->Branch("g4rw_full_fine_piplus_weights",  &g4rw_full_fine_piplus_weights);
    fTree->Branch("g4rw_full_grid_piplus_weights_fake_data",  &g4rw_full_grid_piplus_weights_fake_data);
    fTree->Branch("g4rw_downstream_grid_piplus_weights",  &g4rw_downstream_grid_piplus_weights);
    fTree->Branch("g4rw_full_grid_piminus_weights", &g4rw_full_grid_piminus_weights);
    fTree->Branch("g4rw_full_grid_proton_weights",  &g4rw_full_grid_proton_weights);
    fTree->Branch("g4rw_full_grid_neutron_weights",  &g4rw_full_grid_neutron_weights);
    fTree->Branch("g4rw_full_grid_kplus_weights",  &g4rw_full_grid_kplus_weights);
    fTree->Branch("g4rw_primary_grid_weights", &g4rw_primary_grid_weights);
  }

  //fTree->Branch("g4rw_full_grid_coeffs", &g4rw_full_grid_coeffs);
  fTree->Branch("g4rw_full_grid_piplus_coeffs",  &g4rw_full_grid_piplus_coeffs);
  fTree->Branch("g4rw_full_grid_abscex_coeffs",  &g4rw_full_grid_abscex_coeffs);
  fTree->Branch("g4rw_full_grid_abscex_weights",  &g4rw_full_grid_abscex_weights);
  fTree->Branch("g4rw_primary_grid_abscex_coeffs",  &g4rw_primary_grid_abscex_coeffs);
  fTree->Branch("g4rw_primary_grid_abscex_weights",  &g4rw_primary_grid_abscex_weights);
  fTree->Branch("g4rw_downstream_grid_abscex_coeffs",  &g4rw_downstream_grid_abscex_coeffs);
  fTree->Branch("g4rw_downstream_grid_abscex_weights",  &g4rw_downstream_grid_abscex_weights);

  fTree->Branch("g4rw_full_fine_piplus_coeffs",  &g4rw_full_fine_piplus_coeffs);
  fTree->Branch("g4rw_full_grid_piplus_coeffs_fake_data",  &g4rw_full_grid_piplus_coeffs_fake_data);
  fTree->Branch("g4rw_downstream_grid_piplus_coeffs",  &g4rw_downstream_grid_piplus_coeffs);
  fTree->Branch("g4rw_full_grid_proton_coeffs", &g4rw_full_grid_proton_coeffs);
  fTree->Branch("g4rw_full_grid_neutron_coeffs", &g4rw_full_grid_neutron_coeffs);
  fTree->Branch("g4rw_full_grid_kplus_coeffs", &g4rw_full_grid_kplus_coeffs);
  fTree->Branch("g4rw_primary_grid_coeffs", &g4rw_primary_grid_coeffs);


  fTree->Branch("g4rw_full_grid_piplus_exp_fit_chi2",  &g4rw_full_grid_piplus_exp_fit_chi2);
  fTree->Branch("g4rw_full_fine_piplus_exp_fit_chi2",  &g4rw_full_fine_piplus_exp_fit_chi2);
  fTree->Branch("g4rw_full_grid_piplus_exp_fit_chi2_fake_data",  &g4rw_full_grid_piplus_exp_fit_chi2_fake_data);
  fTree->Branch("g4rw_downstream_grid_piplus_exp_fit_chi2",  &g4rw_downstream_grid_piplus_exp_fit_chi2);
  fTree->Branch("g4rw_full_grid_proton_exp_fit_chi2", &g4rw_full_grid_proton_exp_fit_chi2);
  fTree->Branch("g4rw_full_grid_neutron_exp_fit_chi2", &g4rw_full_grid_neutron_exp_fit_chi2);
  fTree->Branch("g4rw_full_grid_kplus_exp_fit_chi2", &g4rw_full_grid_kplus_exp_fit_chi2);
  fTree->Branch("g4rw_primary_grid_exp_fit_chi2", &g4rw_primary_grid_exp_fit_chi2);
  //fTree->Branch("g4rw_primary_grid_pair_weights", &g4rw_primary_grid_pair_weights);

  fTree->Branch("g4rw_piplus_traj_ps", &g4rw_piplus_traj_ps);
  fTree->Branch("g4rw_piplus_traj_lens", &g4rw_piplus_traj_lens);
  fTree->Branch("g4rw_piplus_traj_npiplus", &g4rw_piplus_traj_npiplus);
  fTree->Branch("g4rw_piplus_traj_npiminus", &g4rw_piplus_traj_npiminus);
  fTree->Branch("g4rw_piplus_traj_npi0", &g4rw_piplus_traj_npi0);
  if( fSaveHits ){
    fTree->Branch( "reco_beam_spacePts_X", &reco_beam_spacePts_X );
    fTree->Branch( "reco_beam_spacePts_Y", &reco_beam_spacePts_Y );
    fTree->Branch( "reco_beam_spacePts_Z", &reco_beam_spacePts_Z );

    fTree->Branch( "reco_daughter_spacePts_X", &reco_daughter_spacePts_X );
    fTree->Branch( "reco_daughter_spacePts_Y", &reco_daughter_spacePts_Y );
    fTree->Branch( "reco_daughter_spacePts_Z", &reco_daughter_spacePts_Z );

    fTree->Branch( "reco_daughter_shower_spacePts_X", &reco_daughter_shower_spacePts_X );
    fTree->Branch( "reco_daughter_shower_spacePts_Y", &reco_daughter_shower_spacePts_Y );
    fTree->Branch( "reco_daughter_shower_spacePts_Z", &reco_daughter_shower_spacePts_Z );
    
    fTree->Branch( "sparsenet_features_charge", &sparsenet_features_charge );
    fTree->Branch( "sparsenet_features_angle", &sparsenet_features_angle );
    fTree->Branch( "sparsenet_features_dot_product", &sparsenet_features_dot_product );
    fTree->Branch( "sparsenet_features_neighboring_nodes_3", &sparsenet_features_neighboring_nodes_3 );
    fTree->Branch( "sparsenet_features_neighboring_nodes_10", &sparsenet_features_neighboring_nodes_10 );
    fTree->Branch( "sparsenet_features_neighboring_nodes_30", &sparsenet_features_neighboring_nodes_30 );
    fTree->Branch( "sparsenet_features_charge_distance_3", &sparsenet_features_charge_distance_3 );
    fTree->Branch( "sparsenet_features_charge_distance_10", &sparsenet_features_charge_distance_10 );
    fTree->Branch( "sparsenet_features_charge_distance_30", &sparsenet_features_charge_distance_30 );
 
  }

}

void pduneana::PDSPAnalyzer::endJob()
{
  dEdX_template_file->Close();
}

double pduneana::PDSPAnalyzer::lateralDist(TVector3 &n, TVector3 &x0, TVector3 &p){
  TVector3 x = ( (p - x0)*n )*n;
  return (x - (p - x0)).Mag();
}

void pduneana::PDSPAnalyzer::reset()
{
  reco_reconstructable_beam_event = false;
  reco_beam_startX = -999;
  reco_beam_startY = -999;
  reco_beam_startZ = -999;
  reco_beam_endX = -999;
  reco_beam_endY = -999;
  reco_beam_endZ = -999;
  reco_beam_flipped = false;
  reco_beam_trackEndDirX = -999;
  reco_beam_trackEndDirY = -999;
  reco_beam_trackEndDirZ = -999;
  reco_beam_trackDirX = -999;
  reco_beam_trackDirY = -999;
  reco_beam_trackDirZ = -999;

  true_beam_len = -999.;
  reco_beam_len = -999;
  test_branch = -1;
  reco_beam_alt_len = -999;
  for_truncation_method = -1;
  reco_beam_alt_len_allTrack = -999;
  reco_beam_calo_startX = -999;
  reco_beam_calo_startY = -999;
  reco_beam_calo_startZ = -999;
  reco_beam_calo_endX = -999;
  reco_beam_calo_endY = -999;
  reco_beam_calo_endZ = -999;
  reco_beam_calo_startX_allTrack = -999;
  reco_beam_calo_startY_allTrack = -999;
  reco_beam_calo_startZ_allTrack = -999;
  reco_beam_calo_endX_allTrack = -999;
  reco_beam_calo_endY_allTrack = -999;
  reco_beam_calo_endZ_allTrack = -999;
  reco_beam_calo_startDirX.clear();
  reco_beam_calo_startDirY.clear();
  reco_beam_calo_startDirZ.clear();
  reco_beam_calo_endDirX.clear();
  reco_beam_calo_endDirY.clear();
  reco_beam_calo_endDirZ.clear();

  reco_track_startX.clear();
  reco_track_startY.clear();
  reco_track_startZ.clear();
  reco_track_endX.clear();
  reco_track_endY.clear();
  reco_track_endZ.clear();
  reco_track_michel_score.clear();
  reco_track_michel_score_weight_by_charge.clear();
  reco_track_ID.clear();
  reco_track_nHits.clear();

  reco_beam_type = -999;
  reco_beam_passes_beam_cuts = false;

  reco_beam_vertex_slice = std::numeric_limits<int>::max();

  true_daughter_nPi0 = 0;
  true_daughter_nPiPlus = 0;
  true_daughter_nPiMinus = 0;
  true_daughter_nProton = 0;
  true_daughter_nNeutron = 0;
  true_daughter_nNucleus = 0;

  reco_beam_true_byE_PDG = -999;
  reco_beam_true_byE_ID = -999;
  reco_beam_true_byHits_PDG = -999;
  reco_beam_true_byHits_ID = -999;

  true_beam_PDG = -999;
  true_beam_mass = -999.;
  true_beam_ID = -999;
  true_beam_hierarchy.clear();
  true_beam_hierarchy_size = -999;
  reco_beam_true_byHits_in_primary_hierarchy = false;
  true_beam_endProcess ="";
  true_beam_endX = -999.;
  true_beam_endY = -999.;
  true_beam_endZ = -999.;
  true_beam_endX_SCE = -999.;
  true_beam_endY_SCE = -999.;
  true_beam_endZ_SCE = -999.;
  true_beam_startX = -999.;
  true_beam_startY = -999.;
  true_beam_startZ = -999.;

  true_beam_startPx   = -999.;
  true_beam_startPy   = -999.;
  true_beam_startPz   = -999.;
  true_beam_startP    = -999.;

  true_beam_endPx   = -999.;
  true_beam_endPy   = -999.;
  true_beam_endPz   = -999.;
  true_beam_endP    = -999.;
  true_beam_endP2    = -999.;
  true_beam_last_len    = -999.;

  true_beam_startDirX = -999.;
  true_beam_startDirY = -999.;
  true_beam_startDirZ = -999.;
  true_beam_nHits = -999;


  true_beam_processes.clear();
  true_beam_nElasticScatters = 0;
  true_beam_elastic_costheta.clear();
  true_beam_elastic_X.clear();
  true_beam_elastic_Y.clear();
  true_beam_elastic_Z.clear();
  true_beam_elastic_deltaE.clear();
  true_beam_elastic_IDE_edep.clear();
  true_beam_IDE_totalDep = -999.;

  true_beam_reco_byHits_PFP_ID.clear();
  true_beam_reco_byHits_PFP_nHits.clear();
  true_beam_reco_byHits_allTrack_ID.clear();

  reco_beam_true_byE_endProcess ="";
  reco_beam_true_byE_process ="";
  reco_beam_true_byE_origin = -999;

  reco_beam_true_byE_endPx = -999.;
  reco_beam_true_byE_endPy = -999.;
  reco_beam_true_byE_endPz = -999.;
  reco_beam_true_byE_endE = -999.;
  reco_beam_true_byE_endP = -999.;

  reco_beam_true_byE_startPx = -999.;
  reco_beam_true_byE_startPy = -999.;
  reco_beam_true_byE_startPz = -999.;
  reco_beam_true_byE_startE = -999.;
  reco_beam_true_byE_startP = -999.;

  reco_beam_true_byHits_endProcess ="";
  reco_beam_true_byHits_process ="";
  reco_beam_true_byHits_origin = -999;

  reco_beam_true_byHits_endPx = -999.;
  reco_beam_true_byHits_endPy = -999.;
  reco_beam_true_byHits_endPz = -999.;
  reco_beam_true_byHits_endE = -999.;
  reco_beam_true_byHits_endP = -999.;

  reco_beam_true_byHits_startPx = -999.;
  reco_beam_true_byHits_startPy = -999.;
  reco_beam_true_byHits_startPz = -999.;
  reco_beam_true_byHits_startE = -999.;
  reco_beam_true_byHits_startP = -999.;

  reco_beam_true_byE_matched = false;
  reco_beam_true_byHits_matched = false;
  reco_beam_true_byHits_purity = -999.;


  //reco_daughter_true_byE_isPrimary = false;
  reco_beam_Chi2_proton = -999.;
  reco_beam_Chi2_muon = -999.;


  beam_inst_P = -999.;
  beam_inst_C0 = -999;
  beam_inst_C1 = -999;
  beam_inst_C0_pressure = -999.;
  beam_inst_C1_pressure = -999.;
  beam_inst_X = -999.;
  beam_inst_Y = -999.;
  beam_inst_Z = -999.;
  beam_inst_dirX = -999.;
  beam_inst_dirY = -999.;
  beam_inst_dirZ = -999.;
  beam_inst_nFibersP1 = -999;
  beam_inst_nFibersP2 = -999;
  beam_inst_nFibersP3 = -999;
  beam_inst_PDG_candidates.clear();
  beam_inst_TOF.clear();
  beam_inst_TOF_Chan.clear();
  beam_inst_nTracks = -999;
  beam_inst_nMomenta = -999;
  beam_inst_valid = true;
  beam_inst_trigger = -999;

  reco_beam_Chi2_ndof = -999;

  reco_daughter_allTrack_momByRange_proton.clear();
  reco_daughter_allTrack_momByRange_muon.clear();
  reco_beam_momByRange_proton = -999.;
  reco_beam_momByRange_muon = -999.;

  reco_daughter_allTrack_momByRange_alt_proton.clear();
  reco_daughter_allTrack_momByRange_alt_muon.clear();
  reco_beam_momByRange_alt_proton = -999.;
  reco_beam_momByRange_alt_muon = -999.;

  true_beam_daughter_PDG.clear();
  true_beam_daughter_len.clear();
  true_beam_daughter_startX.clear();
  true_beam_daughter_startY.clear();
  true_beam_daughter_startZ.clear();
  true_beam_daughter_startPx.clear();
  true_beam_daughter_startPy.clear();
  true_beam_daughter_startPz.clear();
  true_beam_daughter_startP.clear();
  true_beam_daughter_endX.clear();
  true_beam_daughter_endY.clear();
  true_beam_daughter_endZ.clear();
  true_beam_daughter_Process.clear();
  true_beam_daughter_endProcess.clear();
  true_beam_daughter_nHits.clear();

  true_beam_daughter_reco_byHits_PFP_ID.clear();
  true_beam_daughter_reco_byHits_PFP_nHits.clear();
  true_beam_daughter_reco_byHits_PFP_trackScore.clear();

  true_beam_daughter_reco_byHits_allTrack_ID.clear();
  true_beam_daughter_reco_byHits_allTrack_startX.clear();
  true_beam_daughter_reco_byHits_allTrack_startY.clear();
  true_beam_daughter_reco_byHits_allTrack_startZ.clear();
  true_beam_daughter_reco_byHits_allTrack_endX.clear();
  true_beam_daughter_reco_byHits_allTrack_endY.clear();
  true_beam_daughter_reco_byHits_allTrack_endZ.clear();
  true_beam_daughter_reco_byHits_allTrack_len.clear();

  true_beam_daughter_reco_byHits_allShower_ID.clear();
  true_beam_daughter_reco_byHits_allShower_startX.clear();
  true_beam_daughter_reco_byHits_allShower_startY.clear();
  true_beam_daughter_reco_byHits_allShower_startZ.clear();
  true_beam_daughter_reco_byHits_allShower_len.clear();

  true_beam_Pi0_decay_ID.clear();
  true_beam_Pi0_decay_parID.clear();
  true_beam_Pi0_decay_startP.clear();
  true_beam_Pi0_decay_startPx.clear();
  true_beam_Pi0_decay_startPy.clear();
  true_beam_Pi0_decay_startPz.clear();
  true_beam_Pi0_decay_startX.clear();
  true_beam_Pi0_decay_startY.clear();
  true_beam_Pi0_decay_startZ.clear();
  true_beam_Pi0_decay_PDG.clear();
  true_beam_Pi0_decay_len.clear();
  true_beam_Pi0_decay_nHits.clear();
  true_beam_Pi0_decay_reco_byHits_PFP_ID.clear();
  true_beam_Pi0_decay_reco_byHits_PFP_nHits.clear();
  true_beam_Pi0_decay_reco_byHits_PFP_trackScore.clear();

  true_beam_Pi0_decay_reco_byHits_allTrack_ID.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_startX.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_startY.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_startZ.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_endX.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_endY.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_endZ.clear();
  true_beam_Pi0_decay_reco_byHits_allTrack_len.clear();

  true_beam_Pi0_decay_reco_byHits_allShower_ID.clear();
  true_beam_Pi0_decay_reco_byHits_allShower_startX.clear();
  true_beam_Pi0_decay_reco_byHits_allShower_startY.clear();
  true_beam_Pi0_decay_reco_byHits_allShower_startZ.clear();
  true_beam_Pi0_decay_reco_byHits_allShower_len.clear();

  true_beam_grand_daughter_ID.clear();
  true_beam_grand_daughter_parID.clear();
  true_beam_grand_daughter_PDG.clear();
  true_beam_grand_daughter_nHits.clear();
  true_beam_grand_daughter_Process.clear();
  true_beam_grand_daughter_endProcess.clear();
  true_beam_daughter_ID.clear();

  reco_daughter_PFP_ID.clear();
  reco_daughter_pandora_type.clear();

  reco_daughter_PFP_nHits.clear();
  reco_daughter_PFP_nHits_collection.clear();
  reco_daughter_PFP_trackScore.clear();
  reco_daughter_PFP_emScore.clear();
  reco_daughter_PFP_michelScore.clear();
  reco_daughter_PFP_trackScore_collection.clear();
  reco_daughter_PFP_emScore_collection.clear();
  reco_daughter_PFP_michelScore_collection.clear();

  reco_daughter_PFP_trackScore_weight_by_charge.clear();
  reco_daughter_PFP_emScore_weight_by_charge.clear();
  reco_daughter_PFP_michelScore_weight_by_charge.clear();
  reco_daughter_PFP_trackScore_collection_weight_by_charge.clear();
  reco_daughter_PFP_emScore_collection_weight_by_charge.clear();
  reco_daughter_PFP_michelScore_collection_weight_by_charge.clear();

  reco_beam_PFP_ID = -999;
  reco_beam_PFP_nHits = -999;
  reco_beam_PFP_trackScore = -999;
  reco_beam_PFP_emScore = -999;
  reco_beam_PFP_michelScore = -999;
  reco_beam_PFP_trackScore_collection = -999;
  reco_beam_PFP_emScore_collection = -999;
  reco_beam_PFP_michelScore_collection = -999;
  reco_beam_PFP_trackScore_weight_by_charge = -999;
  reco_beam_PFP_emScore_weight_by_charge = -999;
  reco_beam_PFP_michelScore_weight_by_charge = -999;
  reco_beam_PFP_trackScore_collection_weight_by_charge = -999;
  reco_beam_PFP_emScore_collection_weight_by_charge = -999;
  reco_beam_PFP_michelScore_collection_weight_by_charge = -999;

  reco_beam_allTrack_ID = -999;
  reco_beam_allTrack_beam_cuts = -999;
  reco_beam_allTrack_flipped = -999;
  reco_beam_allTrack_len = -999;
  reco_beam_allTrack_startX = -999;
  reco_beam_allTrack_startY = -999;
  reco_beam_allTrack_startZ = -999;
  reco_beam_allTrack_endX = -999;
  reco_beam_allTrack_endY = -999;
  reco_beam_allTrack_endZ = -999;
  reco_beam_allTrack_trackDirX = -999;
  reco_beam_allTrack_trackDirY = -999;
  reco_beam_allTrack_trackDirZ = -999;
  reco_beam_allTrack_trackEndDirX = -999;
  reco_beam_allTrack_trackEndDirY = -999;
  reco_beam_allTrack_trackEndDirZ = -999;
  reco_beam_allTrack_resRange.clear();
  reco_beam_allTrack_calibrated_dEdX.clear();
  reco_beam_allTrack_Chi2_proton = -999;
  reco_beam_allTrack_Chi2_ndof = -999;

  reco_beam_dQdX_NoSCE.clear();
  reco_beam_dEdX_NoSCE.clear();
  reco_beam_calibrated_dEdX_NoSCE.clear();
  reco_beam_resRange_NoSCE.clear();
  reco_beam_TrkPitch_NoSCE.clear();

  reco_beam_dQdX_SCE.clear();
  reco_beam_EField_SCE.clear();
  reco_beam_dQ.clear();
  reco_beam_dQ_NoSCE.clear();
  reco_beam_dEdX_SCE.clear();
  reco_beam_calibrated_dEdX_SCE.clear();
  reco_beam_calibrated_dQdX_SCE.clear();
  reco_beam_vertex_nHits = -999;
  reco_beam_vertex_michel_score = -999.;
  reco_beam_vertex_nHits_allTrack = -999;
  reco_beam_vertex_michel_score_allTrack = -999.;
  reco_beam_vertex_michel_score_weight_by_charge = -999.;
  reco_beam_vertex_michel_score_weight_by_charge_allTrack = -999.;

  reco_beam_resRange_SCE.clear();
  reco_beam_TrkPitch_SCE.clear();
  reco_beam_TrkPitch_SCE_allTrack.clear();
  reco_beam_calo_wire.clear();
  reco_beam_calo_wire_allTrack.clear();
  reco_beam_calo_wire_z.clear();
  reco_beam_calo_X.clear();
  reco_beam_calo_Y.clear();
  reco_beam_calo_Z.clear();
  reco_beam_calo_X_allTrack.clear();
  reco_beam_calo_Y_allTrack.clear();
  reco_beam_calo_Z_allTrack.clear();
  reco_beam_calo_wire_NoSCE.clear();
  reco_beam_calo_wire_z_NoSCE.clear();
  reco_beam_calo_tick.clear();
  reco_beam_calo_TPC.clear();
  reco_beam_calo_TPC_NoSCE.clear();
  reco_beam_hit_IDE_IDs.clear();
  reco_beam_hit_IDE_electrons.clear();
  reco_beam_hit_IDE_energies.clear();
  reco_beam_hit_IDE_origins.clear();
  reco_beam_hit_IDE_cosmic_electrons.clear();
  reco_beam_hit_IDE_beam_electrons.clear();
  reco_beam_hit_IDE_cosmic_energies.clear();
  reco_beam_hit_IDE_beam_energies.clear();

  reco_beam_trackID = -999;

  n_beam_slices =  -999;
  n_beam_particles = -999;
  beam_track_IDs.clear();
  beam_particle_scores.clear();

  reco_beam_incidentEnergies.clear();
  reco_beam_interactingEnergy = -999.;
  reco_beam_incidentEnergies_allTrack.clear();
  reco_beam_interactingEnergy_allTrack = -999.;
  true_beam_incidentEnergies.clear();
  true_beam_slices.clear();
  true_beam_slices_found.clear();
  true_beam_slices_deltaE.clear();
  true_beam_interactingEnergy = -999.;
  em_energy = 0.;
  true_beam_traj_X.clear();
  true_beam_traj_Y.clear();
  true_beam_traj_Z.clear();
  true_beam_traj_Px.clear();
  true_beam_traj_Py.clear();
  true_beam_traj_Pz.clear();
  true_beam_traj_KE.clear();
  true_beam_traj_X_SCE.clear();
  true_beam_traj_Y_SCE.clear();
  true_beam_traj_Z_SCE.clear();
  true_beam_is_scraper = false;

  //Alternative Reco
  reco_daughter_PFP_true_byHits_PDG.clear();
  reco_daughter_PFP_true_byHits_ID.clear();
  reco_daughter_PFP_true_byHits_origin.clear();
  reco_daughter_PFP_true_byHits_parID.clear();
  reco_daughter_PFP_true_byHits_parPDG.clear();
  reco_daughter_PFP_true_byHits_process.clear();
  reco_daughter_PFP_true_byHits_sharedHits.clear();
  reco_daughter_PFP_true_byHits_emHits.clear();

  reco_daughter_PFP_true_byHits_len.clear();
  reco_daughter_PFP_true_byHits_startX.clear();
  reco_daughter_PFP_true_byHits_startY.clear();
  reco_daughter_PFP_true_byHits_startZ.clear();
  reco_daughter_PFP_true_byHits_endX.clear();
  reco_daughter_PFP_true_byHits_endY.clear();
  reco_daughter_PFP_true_byHits_endZ.clear();

  reco_daughter_PFP_true_byHits_startPx.clear();
  reco_daughter_PFP_true_byHits_startPy.clear();
  reco_daughter_PFP_true_byHits_startPz.clear();
  reco_daughter_PFP_true_byHits_startP.clear();
  reco_daughter_PFP_true_byHits_startE.clear();
  reco_daughter_PFP_true_byHits_endProcess.clear();
  reco_daughter_PFP_true_byHits_purity.clear();
  reco_daughter_PFP_true_byHits_completeness.clear();

  reco_daughter_PFP_true_byE_PDG.clear();
  reco_daughter_PFP_true_byE_len.clear();
  reco_daughter_PFP_true_byE_completeness.clear();
  reco_daughter_PFP_true_byE_purity.clear();



  reco_daughter_allTrack_ID.clear();
  reco_daughter_allTrack_dEdX_SCE.clear();
  reco_daughter_allTrack_dQdX_SCE.clear();
  reco_daughter_allTrack_calibrated_dQdX_SCE.clear();
  reco_daughter_allTrack_EField_SCE.clear();
  reco_daughter_allTrack_resRange_SCE.clear();


  //Calorimetry + chi2 for planes 0 and 1
  reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.clear();
  reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.clear();

  reco_daughter_allTrack_resRange_plane0.clear();
  reco_daughter_allTrack_resRange_plane1.clear();

  reco_daughter_allTrack_Chi2_proton_plane0.clear();
  reco_daughter_allTrack_Chi2_proton_plane1.clear();

  reco_daughter_allTrack_Chi2_ndof_plane0.clear();
  reco_daughter_allTrack_Chi2_ndof_plane1.clear();
  ///////////////////////////////////////////


  //reco_daughter_allTrack_calibrated_dEdX.clear();
  reco_daughter_allTrack_calibrated_dEdX_SCE.clear();

  reco_daughter_allTrack_Chi2_proton.clear();
  reco_daughter_allTrack_Chi2_muon.clear();
  reco_daughter_allTrack_Chi2_pion.clear();
  reco_daughter_allTrack_Chi2_ndof.clear();
  reco_daughter_allTrack_Chi2_ndof_muon.clear();
  reco_daughter_allTrack_Chi2_ndof_pion.clear();

  reco_daughter_allTrack_Theta.clear();
  reco_daughter_allTrack_startDirX.clear();
  reco_daughter_allTrack_startDirY.clear();
  reco_daughter_allTrack_startDirZ.clear();
  reco_daughter_allTrack_Phi.clear();
  reco_daughter_allTrack_len.clear();
  reco_daughter_allTrack_alt_len.clear();
  reco_daughter_allTrack_startX.clear();
  reco_daughter_allTrack_startY.clear();
  reco_daughter_allTrack_startZ.clear();
  reco_daughter_allTrack_endX.clear();
  reco_daughter_allTrack_endY.clear();
  reco_daughter_allTrack_endZ.clear();
  reco_daughter_allTrack_calo_X.clear();
  reco_daughter_allTrack_calo_Y.clear();
  reco_daughter_allTrack_calo_Z.clear();
  reco_daughter_allTrack_vertex_michel_score.clear();
  reco_daughter_allTrack_vertex_nHits.clear();

  reco_daughter_allShower_ID.clear();
  reco_daughter_allShower_len.clear();
  reco_daughter_allShower_startX.clear();
  reco_daughter_allShower_startY.clear();
  reco_daughter_allShower_startZ.clear();

  reco_daughter_allShower_dirX.clear();
  reco_daughter_allShower_dirY.clear();
  reco_daughter_allShower_dirZ.clear();
  reco_daughter_allShower_energy.clear();
  reco_daughter_allShower_calibrated_energy.clear();
  ///////


  //New Hits info
  reco_beam_spacePts_X.clear();
  reco_beam_spacePts_Y.clear();
  reco_beam_spacePts_Z.clear();

  reco_daughter_spacePts_X.clear();
  reco_daughter_spacePts_Y.clear();
  reco_daughter_spacePts_Z.clear();

  reco_daughter_shower_spacePts_X.clear();
  reco_daughter_shower_spacePts_Y.clear();
  reco_daughter_shower_spacePts_Z.clear();
 
  // SparseNet
  sparsenet_features_charge.clear();
  sparsenet_features_angle.clear();
  sparsenet_features_dot_product.clear();
  sparsenet_features_neighboring_nodes_3.clear();
  sparsenet_features_neighboring_nodes_10.clear();
  sparsenet_features_neighboring_nodes_30.clear();
  sparsenet_features_charge_distance_3.clear();
  sparsenet_features_charge_distance_10.clear();
  sparsenet_features_charge_distance_30.clear();
  //

  //g4rw_primary_plus_sigma_weight.clear();
  //g4rw_primary_minus_sigma_weight.clear();
  //g4rw_primary_var.clear();
  //g4rw_alt_primary_plus_sigma_weight.clear();
  //g4rw_alt_primary_minus_sigma_weight.clear();
  //g4rw_full_primary_plus_sigma_weight.clear();
  //g4rw_full_primary_minus_sigma_weight.clear();
  //g4rw_full_grid_weights.clear();
  //g4rw_full_grid_coeffs.clear();
  g4rw_full_grid_piplus_weights.clear();
  g4rw_full_grid_piplus_coeffs.clear();
  g4rw_full_grid_abscex_weights.clear();
  g4rw_full_grid_abscex_coeffs.clear();
  g4rw_primary_grid_abscex_weights.clear();
  g4rw_primary_grid_abscex_coeffs.clear();
  g4rw_downstream_grid_abscex_weights.clear();
  g4rw_downstream_grid_abscex_coeffs.clear();

  g4rw_full_fine_piplus_weights.clear();
  g4rw_full_fine_piplus_coeffs.clear();
  g4rw_full_grid_piplus_weights_fake_data.clear();
  g4rw_full_grid_piplus_coeffs_fake_data.clear();
  g4rw_downstream_grid_piplus_weights.clear();
  g4rw_downstream_grid_piplus_coeffs.clear();
  g4rw_full_grid_piminus_weights.clear();
  g4rw_full_grid_proton_weights.clear();
  g4rw_full_grid_proton_coeffs.clear();
  g4rw_full_grid_neutron_weights.clear();
  g4rw_full_grid_neutron_coeffs.clear();
  g4rw_full_grid_kplus_weights.clear();
  g4rw_full_grid_kplus_coeffs.clear();
  g4rw_primary_grid_weights.clear();
  g4rw_primary_grid_coeffs.clear();
  //g4rw_primary_grid_pair_weights.clear();
  //
  g4rw_primary_grid_exp_fit_chi2.clear();
  g4rw_full_grid_piplus_exp_fit_chi2.clear();
  g4rw_full_grid_piplus_exp_fit_chi2_fake_data.clear();
  g4rw_full_grid_proton_exp_fit_chi2.clear();
  g4rw_full_grid_neutron_exp_fit_chi2.clear();
  g4rw_full_grid_kplus_exp_fit_chi2.clear();
  g4rw_downstream_grid_piplus_exp_fit_chi2.clear();
  g4rw_full_fine_piplus_exp_fit_chi2.clear();

  g4rw_piplus_traj_ps.clear();
  g4rw_piplus_traj_lens.clear();
  g4rw_piplus_traj_npiplus.clear();
  g4rw_piplus_traj_npiminus.clear();
  g4rw_piplus_traj_npi0.clear();
}


void pduneana::PDSPAnalyzer::BeamPFPInfo(
    const art::Event & evt, const recob::PFParticle* particle,
    anab::MVAReader<recob::Hit,4> * hitResults) {

  //Get CNN output for the beam
  reco_beam_PFP_ID = particle->Self();
  const std::vector< art::Ptr< recob::Hit > > beamPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( *particle, evt, fPFParticleTag );
  reco_beam_PFP_nHits = beamPFP_hits.size();
  if (!fSkipMVA) {
    cnnOutput2D cnn = GetCNNOutputFromPFParticle( *particle, evt, *hitResults, pfpUtil, fPFParticleTag );
    reco_beam_PFP_trackScore = (cnn.nHits > 0 ?
                                (cnn.track / cnn.nHits) : -999.);
    reco_beam_PFP_emScore = (cnn.nHits > 0 ?
                             (cnn.em / cnn.nHits) : -999.);
    reco_beam_PFP_michelScore = (cnn.nHits > 0 ?
                                 (cnn.michel / cnn.nHits) : -999.);
    reco_beam_PFP_trackScore_weight_by_charge = cnn.track_weight_by_charge;
    reco_beam_PFP_emScore_weight_by_charge = cnn.em_weight_by_charge;
    reco_beam_PFP_michelScore_weight_by_charge = cnn.michel_weight_by_charge;

    cnnOutput2D cnn_collection = GetCNNOutputFromPFParticleFromPlane( *particle, evt, *hitResults, pfpUtil, fPFParticleTag, 2 );
    reco_beam_PFP_trackScore_collection = (cnn_collection.nHits > 0 ?
                                           cnn_collection.track / cnn_collection.nHits : -999.);
    reco_beam_PFP_emScore_collection = (cnn_collection.nHits > 0 ?
                                        cnn_collection.em / cnn_collection.nHits : -999.);
    reco_beam_PFP_michelScore_collection = (cnn_collection.nHits > 0 ?
                                            cnn_collection.michel / cnn_collection.nHits : -999.);
    reco_beam_PFP_trackScore_collection_weight_by_charge = cnn_collection.track_weight_by_charge;
    reco_beam_PFP_emScore_collection_weight_by_charge = cnn_collection.em_weight_by_charge;
    reco_beam_PFP_michelScore_collection_weight_by_charge = cnn_collection.michel_weight_by_charge;
    
  }

  const auto space_pts = pfpUtil.GetPFParticleSpacePoints(*particle, evt, fPFParticleTag);
  for (const auto * sp : space_pts) {
    reco_beam_spacePts_X.push_back(sp->XYZ()[0]);
    reco_beam_spacePts_Y.push_back(sp->XYZ()[1]);
    reco_beam_spacePts_Z.push_back(sp->XYZ()[2]);
    
  }
}

void pduneana::PDSPAnalyzer::BeamTrackInfo(
    const art::Event & evt, const recob::Track * thisTrack,
    detinfo::DetectorClocksData const& clockData) {

  auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto const detProp
    = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(
        evt, clockData);
  auto allHits = evt.getValidHandle<std::vector<recob::Hit> >(fHitTag);

  //Ajib's michel identifier
  if (!fSkipMVA) {
    std::pair<double, int> vertex_michel_score =
        trackUtil.GetVertexMichelScore(*thisTrack, evt, fTrackerTag, fHitTag/*,
                                       0., -500., 500., 0., 500., 0.*/);
    reco_beam_vertex_nHits = vertex_michel_score.second;
    reco_beam_vertex_michel_score = vertex_michel_score.first;
    
    std::pair<double, double> vertex_michel_score_weight_by_charge =
        trackUtil.GetVertexMichelScore_weight_by_charge(*thisTrack, evt, fTrackerTag, fHitTag);
    reco_beam_vertex_michel_score_weight_by_charge = (vertex_michel_score_weight_by_charge.second != 0 ? vertex_michel_score_weight_by_charge.first/vertex_michel_score_weight_by_charge.second : -999.);
  }

  if (fVerbose) std::cout << "Beam particle is track-like " << thisTrack->ID() << std::endl;
  reco_beam_type = 13;

  //Flag to tell if it passes the beam cuts
  //-- should probably take out or make it more flexible for other runs
  reco_beam_passes_beam_cuts = beam_cuts.IsBeamlike( *thisTrack, evt, "1" );
  if (fVerbose) std::cout << "Beam Cuts " << reco_beam_passes_beam_cuts << std::endl;


  reco_beam_trackID = thisTrack->ID();

  //Using default pandora info -- not SCE-corrected
  reco_beam_startX = thisTrack->Trajectory().Start().X();
  reco_beam_startY = thisTrack->Trajectory().Start().Y();
  reco_beam_startZ = thisTrack->Trajectory().Start().Z();
  reco_beam_endX = thisTrack->Trajectory().End().X();
  reco_beam_endY = thisTrack->Trajectory().End().Y();
  reco_beam_endZ = thisTrack->Trajectory().End().Z();

  auto startDir = thisTrack->StartDirection();
  auto endDir   = thisTrack->EndDirection();

  //try flipping
  if( reco_beam_startZ > reco_beam_endZ ){
    reco_beam_flipped = true;
    reco_beam_endX = thisTrack->Trajectory().Start().X();
    reco_beam_endY = thisTrack->Trajectory().Start().Y();
    reco_beam_endZ = thisTrack->Trajectory().Start().Z();
    reco_beam_startX = thisTrack->Trajectory().End().X();
    reco_beam_startY = thisTrack->Trajectory().End().Y();
    reco_beam_startZ = thisTrack->Trajectory().End().Z();

    reco_beam_trackDirX =  -1. * endDir.X();
    reco_beam_trackDirY =  -1. * endDir.Y();
    reco_beam_trackDirZ =  -1. * endDir.Z();

    reco_beam_trackEndDirX =  -1. * startDir.X();
    reco_beam_trackEndDirY =  -1. * startDir.Y();
    reco_beam_trackEndDirZ =  -1. * startDir.Z();
  }
  else{
    reco_beam_flipped = false;
    reco_beam_trackDirX    =  startDir.X();
    reco_beam_trackDirY    =  startDir.Y();
    reco_beam_trackDirZ    =  startDir.Z();
    reco_beam_trackEndDirX =  endDir.X();
    reco_beam_trackEndDirY =  endDir.Y();
    reco_beam_trackEndDirZ =  endDir.Z();
  }

  reco_beam_len  = thisTrack->Length();
  trkf::TrackMomentumCalculator track_p_calc;
  
  //Calculates momentum according to CSDA range
  reco_beam_momByRange_proton = track_p_calc.GetTrackMomentum(
      thisTrack->Length(), 2212);
  reco_beam_momByRange_muon = track_p_calc.GetTrackMomentum(
      thisTrack->Length(), 13);
  ////////////////////////////////////////////////////////////////


  //To do: remove
  std::map< const recob::Hit *, int > hitsToSlices;
  std::map< int, std::vector< const recob::Hit * > > slicesToHits;

  //Looking at the hits in the beam track
  std::map< size_t, const recob::Hit * > trajPtsToHits = trackUtil.GetRecoHitsFromTrajPoints( *thisTrack, evt, fTrackerTag );
  //double max_X = 0.;
  //double max_Y = 0.;
  //double max_Z = 0.;

  for( auto it = trajPtsToHits.begin(); it != trajPtsToHits.end(); ++it ){

    const recob::Hit * theHit = it->second;
    size_t i = it->first;

    //double x = thisTrack->Trajectory().LocationAtPoint(i).X();
    //double y = thisTrack->Trajectory().LocationAtPoint(i).Y();
    double z = thisTrack->Trajectory().LocationAtPoint(i).Z();

    //if( fSaveHits ){
    //  //saving all hit coordinates for beamtrack
    //  reco_beam_spacePts_X.push_back(x);
    //  reco_beam_spacePts_Y.push_back(y);
    //  reco_beam_spacePts_Z.push_back(z);
    //}

    //This creates the slices for the thin slice method.
    //To do: remove this
    int slice = std::floor((z - fZ0) / fPitch);
    hitsToSlices[theHit] = slice;
    slicesToHits[slice].push_back(theHit);
  }

  //Last point is the vertex slice
  //To do: remove
  reco_beam_vertex_slice = slicesToHits.rbegin()->first;

  //Primary Track Calorimetry
  //SCE-corrected
  auto calo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag,
                                                fCalorimetryTagSCE);
  //Get the correct plane (collection) because it's not always the same
  bool found_calo = false;
  size_t index = 0;
  for ( index = 0; index < calo.size(); ++index) {
    if (calo[index].PlaneID().Plane == 2) {
      found_calo = true;
      break; 
    }
  }

  if (found_calo) {
    //Using SCE-corrected (from calorimetry) to calculate momentum by range
    reco_beam_momByRange_alt_proton = track_p_calc.GetTrackMomentum(
        calo[index].Range(), 2212);
    reco_beam_momByRange_alt_muon = track_p_calc.GetTrackMomentum(
        calo[index].Range(), 13);
    reco_beam_alt_len = calo[index].Range();

    auto calo_dQdX = calo[index].dQdx();
    auto calo_dEdX = calo[index].dEdx();
    auto calo_range = calo[index].ResidualRange();
    auto TpIndices = calo[index].TpIndices();

    //For Prod3 only
    if (fCalorimetryTagSCE == "pandoracali") {
      auto pandoracalo = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackerTag, "pandoracalo");
      size_t this_index = 0;
      for ( this_index = 0; this_index < pandoracalo.size(); ++this_index) {
        if (pandoracalo[this_index].PlaneID().Plane == 2) {
          break; 
        }
      }
      TpIndices = pandoracalo[this_index].TpIndices();
    }
    std::cout << calo_dQdX.size() << std::endl;
    std::cout << calo[index].PlaneID().Plane << std::endl;

    auto theXYZPoints = calo[index].XYZ();
    std::vector< size_t > calo_hit_indices;
    std::vector<std::vector<int>> reco_hit_IDE_IDs;
    std::vector<std::vector<int>> reco_hit_IDE_origins;
    std::vector<std::vector<double>> reco_hit_IDE_electrons;
    std::vector<std::vector<double>> reco_hit_IDE_energies;
    for( size_t i = 0; i < calo_dQdX.size(); ++i ){
      if (fVerbose) std::cout << i << std::endl;
      reco_beam_dQdX_SCE.push_back( calo_dQdX[i] );
      reco_beam_dEdX_SCE.push_back( calo_dEdX[i] );
      reco_beam_resRange_SCE.push_back( calo_range[i] );
      reco_beam_TrkPitch_SCE.push_back( calo[index].TrkPitchVec()[i] );

      const recob::Hit & theHit = (*allHits)[ TpIndices[i] ];
      reco_beam_dQ.push_back(theHit.Integral());
      reco_beam_calo_TPC.push_back(theHit.WireID().TPC);
      if (theHit.WireID().TPC == 1) {
        reco_beam_calo_wire.push_back( theHit.WireID().Wire );
      }
      else if (theHit.WireID().TPC == 5) {
        reco_beam_calo_wire.push_back( theHit.WireID().Wire + 479);
      }
      //Need other TPCs?
      else {
        reco_beam_calo_wire.push_back(theHit.WireID().Wire );
      }
      reco_beam_calo_tick.push_back( theHit.PeakTime() );
      calo_hit_indices.push_back( TpIndices[i] );

      reco_beam_calo_wire_z.push_back(
          wireReadout.Wire(theHit.WireID()).GetCenter().Z());
      reco_beam_calo_X.push_back(theXYZPoints[i].X());
      reco_beam_calo_Y.push_back(theXYZPoints[i].Y());
      reco_beam_calo_Z.push_back(theXYZPoints[i].Z());
      if (fVerbose)
        std::cout << theXYZPoints[i].X() << " " << theXYZPoints[i].Y() << " " <<
                     theXYZPoints[i].Z() << " " << theHit.WireID().Wire << " " <<
                     wireReadout.Wire(theHit.WireID()).GetCenter().Z() << " " <<
                     theHit.WireID().TPC << " " << std::endl;

      //truth infos
      if (!evt.isRealData()) {
        reco_hit_IDE_IDs.push_back(std::vector<int>());
        reco_hit_IDE_origins.push_back(std::vector<int>());
        reco_hit_IDE_electrons.push_back(std::vector<double>());
        reco_hit_IDE_energies.push_back(std::vector<double>());
        for (const auto & ide : bt_serv->HitToEveTrackIDEs(clockData, theHit)) {
          int track_ID = ide.trackID;
          int origin = pi_serv->TrackIdToMCTruth_P(track_ID)->Origin();
          float n_electrons = ide.numElectrons;

          reco_hit_IDE_IDs.back().push_back(track_ID);
          reco_hit_IDE_origins.back().push_back(origin);
          reco_hit_IDE_electrons.back().push_back(n_electrons);
          reco_hit_IDE_energies.back().push_back(ide.energy);
        }
      }
    }

    //Getting the SCE corrected start/end positions & directions
    std::sort(theXYZPoints.begin(), theXYZPoints.end(), [](auto a, auto b)
        {return (a.Z() < b.Z());});

    //Getting position + direction from SCE-corrected calorimetry
    //Attempt to use a few different ways to do this
    //
    //i.e. just gitting the line from start to end
    //or the line between first two points
    //or fitting a line to the first 3-4 points
    if (theXYZPoints.size()) {
      reco_beam_calo_startX = theXYZPoints[0].X();
      reco_beam_calo_startY = theXYZPoints[0].Y();
      reco_beam_calo_startZ = theXYZPoints[0].Z();
      reco_beam_calo_endX = theXYZPoints.back().X();
      reco_beam_calo_endY = theXYZPoints.back().Y();
      reco_beam_calo_endZ = theXYZPoints.back().Z();

      TVector3 dir((theXYZPoints.back().X() - theXYZPoints[0].X()),
                   (theXYZPoints.back().Y() - theXYZPoints[0].Y()),
                   (theXYZPoints.back().Z() - theXYZPoints[0].Z()));
      reco_beam_calo_startDirX.push_back(dir.Unit().X());
      reco_beam_calo_endDirX.push_back(dir.Unit().X());
      reco_beam_calo_startDirY.push_back(dir.Unit().Y());
      reco_beam_calo_endDirY.push_back(dir.Unit().Y());
      reco_beam_calo_startDirZ.push_back(dir.Unit().Z());
      reco_beam_calo_endDirZ.push_back(dir.Unit().Z());
    }
    else {
      reco_beam_calo_startDirX.push_back(-999.);
      reco_beam_calo_endDirX.push_back(-999.);
      reco_beam_calo_startDirY.push_back(-999.);
      reco_beam_calo_endDirY.push_back(-999.);
      reco_beam_calo_startDirZ.push_back(-999.);
      reco_beam_calo_endDirZ.push_back(-999.);
    }

    if (theXYZPoints.size() > 1) {
      TVector3 start_p1(theXYZPoints[0].X(),
          theXYZPoints[0].Y(), theXYZPoints[0].Z());
      TVector3 start_p2(theXYZPoints[1].X(),
          theXYZPoints[1].Y(), theXYZPoints[1].Z());
      TVector3 start_diff = start_p2 - start_p1;

      reco_beam_calo_startDirX.push_back(start_diff.Unit().X());
      reco_beam_calo_startDirY.push_back(start_diff.Unit().Y());
      reco_beam_calo_startDirZ.push_back(start_diff.Unit().Z());

      size_t nPoints = theXYZPoints.size();
      TVector3 end_p1(theXYZPoints[nPoints - 2].X(),
          theXYZPoints[nPoints - 2].Y(), theXYZPoints[nPoints - 2].Z());
      TVector3 end_p2(theXYZPoints[nPoints - 1].X(),
          theXYZPoints[nPoints - 1].Y(), theXYZPoints[nPoints - 1].Z());
      TVector3 end_diff = end_p2 - end_p1;

      reco_beam_calo_endDirX.push_back(end_diff.Unit().X());
      reco_beam_calo_endDirY.push_back(end_diff.Unit().Y());
      reco_beam_calo_endDirZ.push_back(end_diff.Unit().Z());
    }
    else {
      reco_beam_calo_startDirX.push_back(-999.);
      reco_beam_calo_endDirX.push_back(-999.);
      reco_beam_calo_startDirY.push_back(-999.);
      reco_beam_calo_endDirY.push_back(-999.);
      reco_beam_calo_startDirZ.push_back(-999.);
      reco_beam_calo_endDirZ.push_back(-999.);
    }

    if (theXYZPoints.size() > 2) {
      std::vector<TVector3> input;
      for (size_t iP = 0; iP < 3; ++iP) {
        input.push_back(TVector3(theXYZPoints[iP].X(),
                                 theXYZPoints[iP].Y(),
                                 theXYZPoints[iP].Z()));
      }

      TVector3 startDiff = FitLine(input);
      reco_beam_calo_startDirX.push_back(startDiff.Unit().X());
      reco_beam_calo_startDirY.push_back(startDiff.Unit().Y());
      reco_beam_calo_startDirZ.push_back(startDiff.Unit().Z());

      std::vector<TVector3> end_input;
      size_t nPoints = theXYZPoints.size();
      for (size_t iP = nPoints - 3; iP < nPoints; ++iP) {
        end_input.push_back(TVector3(theXYZPoints[iP].X(),
                                     theXYZPoints[iP].Y(),
                                     theXYZPoints[iP].Z()));
      }

      TVector3 endDiff = FitLine(end_input);
      reco_beam_calo_endDirX.push_back(endDiff.Unit().X());
      reco_beam_calo_endDirY.push_back(endDiff.Unit().Y());
      reco_beam_calo_endDirZ.push_back(endDiff.Unit().Z());
    }
    else {
      reco_beam_calo_startDirX.push_back(-999.);
      reco_beam_calo_endDirX.push_back(-999.);
      reco_beam_calo_startDirY.push_back(-999.);
      reco_beam_calo_endDirY.push_back(-999.);
      reco_beam_calo_startDirZ.push_back(-999.);
      reco_beam_calo_endDirZ.push_back(-999.);
    }

    if (theXYZPoints.size() > 3) {
      std::vector<TVector3> input;
      for (size_t iP = 0; iP < 4; ++iP) {
        input.push_back(TVector3(theXYZPoints[iP].X(),
                                 theXYZPoints[iP].Y(),
                                 theXYZPoints[iP].Z()));
      }

      TVector3 startDiff = FitLine(input);
      reco_beam_calo_startDirX.push_back(startDiff.Unit().X());
      reco_beam_calo_startDirY.push_back(startDiff.Unit().Y());
      reco_beam_calo_startDirZ.push_back(startDiff.Unit().Z());

      std::vector<TVector3> end_input;
      size_t nPoints = theXYZPoints.size();
      for (size_t iP = nPoints - 4; iP < nPoints; ++iP) {
        end_input.push_back(TVector3(theXYZPoints[iP].X(),
                                     theXYZPoints[iP].Y(),
                                     theXYZPoints[iP].Z()));
      }

      TVector3 endDiff = FitLine(end_input);
      reco_beam_calo_endDirX.push_back(endDiff.Unit().X());
      reco_beam_calo_endDirY.push_back(endDiff.Unit().Y());
      reco_beam_calo_endDirZ.push_back(endDiff.Unit().Z());

    }
    else {
      reco_beam_calo_startDirX.push_back(-999.);
      reco_beam_calo_endDirX.push_back(-999.);
      reco_beam_calo_startDirY.push_back(-999.);
      reco_beam_calo_endDirY.push_back(-999.);
      reco_beam_calo_startDirZ.push_back(-999.);
      reco_beam_calo_endDirZ.push_back(-999.);
    }
    ////////////////////////////////////////////

    //New Calibration
    std::cout << "Getting reco beam calo" << std::endl;
    if (fRecalibrate){
      std::vector< float > new_dEdX = calibration_SCE.GetCalibratedCalorimetry(*thisTrack, evt, fTrackerTag, fCalorimetryTagSCE, 2, -10.);
      std::cout << new_dEdX.size() << " " << reco_beam_resRange_SCE.size() << std::endl;
      for( size_t i = 0; i < new_dEdX.size(); ++i ){ reco_beam_calibrated_dEdX_SCE.push_back( new_dEdX[i] ); }
      //std::cout << "got calibrated dedx" << std::endl;

      std::vector<double> new_dQdX = calibration_SCE.CalibratedQdX(
        *thisTrack, evt, fTrackerTag,
        fCalorimetryTagSCE, 2, -10.);
      for (auto dqdx : new_dQdX) {
        reco_beam_calibrated_dQdX_SCE.push_back(dqdx);
      }

      std::vector<double> efield = calibration_SCE.GetEFieldVector(
        *thisTrack, evt, fTrackerTag, fCalorimetryTagSCE, 2, -10.);
      for (auto ef : efield) {
        reco_beam_EField_SCE.push_back(ef);
      }
    }
    else{
      for (size_t i = 0; i < calo_dQdX.size(); ++i){
        reco_beam_calibrated_dEdX_SCE.push_back(calo_dEdX[i]);
        reco_beam_calibrated_dQdX_SCE.push_back(calo_dQdX[i]);

        // Electric Field in the drift region in KV/cm
        double E_field_nominal = detProp.Efield();

        geo::Vector_t E_field_offsets
            = (sce->EnableCalEfieldSCE() && fSCE ?
               sce->GetCalEfieldOffsets(
                   geo::Point_t{theXYZPoints[i].X(),
                                theXYZPoints[i].Y(),
                                theXYZPoints[i].Z()},
                   calo[index].PlaneID().TPC) :
               geo::Vector_t{0., 0., 0.});
        TVector3 E_field_vector = {
            E_field_nominal*(1 + E_field_offsets.X()),
            E_field_nominal*E_field_offsets.Y(),
            E_field_nominal*E_field_offsets.Z()};

        double E_field = E_field_vector.Mag();
        reco_beam_EField_SCE.push_back(E_field);
      }
    }
    ////////////////////////////////////////////

    //Get the chi2-based PID for the SCE-corrected beam track
    std::pair<double, int> pid_chi2_ndof = trackUtil.Chi2PID(
        reco_beam_calibrated_dEdX_SCE, reco_beam_resRange_SCE, templates[2212]);
    reco_beam_Chi2_proton = pid_chi2_ndof.first;
    reco_beam_Chi2_ndof = pid_chi2_ndof.second;

    pid_chi2_ndof = trackUtil.Chi2PID(
        reco_beam_calibrated_dEdX_SCE, reco_beam_resRange_SCE, templates[13]);
    reco_beam_Chi2_muon = pid_chi2_ndof.first;

    if (fVerbose)
      std::cout << "Calo check: " << reco_beam_calibrated_dEdX_SCE.size() << " " <<
                   reco_beam_TrkPitch_SCE.size() << std::endl;
    std::vector< calo_point > reco_beam_calo_points;
    //Doing thin slice method
    if (reco_beam_calibrated_dEdX_SCE.size() &&
        reco_beam_calibrated_dEdX_SCE.size() == reco_beam_TrkPitch_SCE.size() &&
        reco_beam_calibrated_dEdX_SCE.size() == reco_beam_calo_wire.size()) {

      for( size_t i = 0; i < reco_beam_calibrated_dEdX_SCE.size(); ++i ){
        reco_beam_calo_points.push_back(
          calo_point(reco_beam_calo_wire[i], reco_beam_calo_tick[i],
                     reco_beam_TrkPitch_SCE[i],
                     reco_beam_dQdX_SCE[i], reco_beam_dEdX_SCE[i],
                     reco_beam_dQ[i],
                     reco_beam_calibrated_dQdX_SCE[i],
                     reco_beam_calibrated_dEdX_SCE[i],
                     reco_beam_resRange_SCE[i], calo_hit_indices[i],
                     reco_beam_calo_wire_z[i], reco_beam_calo_TPC[i],
                     reco_beam_EField_SCE[i], reco_beam_calo_X[i],
                     reco_beam_calo_Y[i], reco_beam_calo_Z[i]));
        if (!evt.isRealData()) {
          reco_beam_calo_points.back().SetIDE_IDs(reco_hit_IDE_IDs[i]);
          reco_beam_calo_points.back().SetIDE_origins(reco_hit_IDE_origins[i]);
          reco_beam_calo_points.back().SetIDE_electrons(reco_hit_IDE_electrons[i]);
          reco_beam_calo_points.back().SetIDE_energies(reco_hit_IDE_energies[i]);
        }
      }

      //std::cout << "N Calo points: " << reco_beam_calo_points.size() << std::endl;
      //Sort
      std::sort( reco_beam_calo_points.begin(), reco_beam_calo_points.end(), [](calo_point a, calo_point b) {return ( a.z < b.z );} );

      //And also put these in the right order
      for( size_t i = 0; i < reco_beam_calo_points.size(); ++i ){
        calo_point thePoint = reco_beam_calo_points[i];
        reco_beam_calo_wire[i] = thePoint.wire;
        reco_beam_calo_tick[i] = thePoint.tick;
        reco_beam_calibrated_dQdX_SCE[i] = thePoint.calibrated_dQdX;
        reco_beam_calibrated_dEdX_SCE[i] = thePoint.calibrated_dEdX;
        reco_beam_TrkPitch_SCE[i] = thePoint.pitch;
        calo_hit_indices[i] = thePoint.hit_index;
        reco_beam_calo_wire_z[i] = thePoint.wire_z;
        reco_beam_calo_TPC[i] = thePoint.tpc;
        reco_beam_resRange_SCE[i] = thePoint.res_range;
        reco_beam_dQdX_SCE[i] = thePoint.dQdX;
        reco_beam_dQ[i] = thePoint.dQ;
        reco_beam_dEdX_SCE[i] = thePoint.dEdX;
        reco_beam_EField_SCE[i] = thePoint.EField;
        reco_beam_calo_X[i] = thePoint.x;
        reco_beam_calo_Y[i] = thePoint.y;
        reco_beam_calo_Z[i] = thePoint.z;
        //reco_beam_calo_x[i]
        if (!evt.isRealData() && fSaveHitIDEInfo) {
          reco_beam_hit_IDE_IDs.push_back(std::vector<int>());
          reco_beam_hit_IDE_origins.push_back(std::vector<int>());
          reco_beam_hit_IDE_electrons.push_back(std::vector<double>());
          reco_beam_hit_IDE_energies.push_back(std::vector<double>());
          double total_cosmic = 0.;
          double total_beam = 0.;
          double total_cosmic_energy = 0.;
          double total_beam_energy = 0.;
          for (size_t j = 0; j < thePoint.IDE_IDs.size(); ++j) {
            reco_beam_hit_IDE_IDs.back().push_back(thePoint.IDE_IDs[j]);
            reco_beam_hit_IDE_electrons.back().push_back(thePoint.IDE_electrons[j]);
            reco_beam_hit_IDE_energies.back().push_back(thePoint.IDE_energies[j]);
            reco_beam_hit_IDE_origins.back().push_back(thePoint.IDE_origins[j]);
            if (thePoint.IDE_origins[j] == 2) {
              total_cosmic += thePoint.IDE_electrons[j];
              total_cosmic_energy += thePoint.IDE_energies[j];
            }
            else if (thePoint.IDE_origins[j] == 4) {
              total_beam += thePoint.IDE_electrons[j];
              total_beam_energy += thePoint.IDE_energies[j];
            }
          }
          reco_beam_hit_IDE_cosmic_electrons.push_back(total_cosmic);
          reco_beam_hit_IDE_beam_electrons.push_back(total_beam);
          reco_beam_hit_IDE_cosmic_energies.push_back(total_cosmic_energy);
          reco_beam_hit_IDE_beam_energies.push_back(total_beam_energy);
        }
      }

      //Get the initial Energy KE
      //Assumes it's a pion -- maybe make configurable To do
      //double mass = 0.;
      double init_KE = 0.;
      //std::cout << "Has BI? " << fMCHasBI << " " << evt.isRealData() << std::endl;
      if (evt.isRealData() || fMCHasBI) {
        double mass = 139.57;

        init_KE =  sqrt( 1.e6*beam_inst_P*beam_inst_P + mass*mass ) - mass;
       // std::cout << "MC has BI: " << init_KE << std::endl;
      }
      else{
        init_KE = sqrt(1.e6*true_beam_startP*true_beam_startP +
                       true_beam_mass*true_beam_mass) - true_beam_mass;
      }

      reco_beam_incidentEnergies.push_back( init_KE );
      for( size_t i = 0; i < reco_beam_calo_points.size() - 1; ++i ){ //-1 to not count the last slice
        //use dedx * pitch or new hit calculation?
        if (reco_beam_calo_points[i].calibrated_dEdX < 0.) continue;
        double this_energy = reco_beam_incidentEnergies.back() - ( reco_beam_calo_points[i].calibrated_dEdX * reco_beam_calo_points[i].pitch );
        reco_beam_incidentEnergies.push_back( this_energy );
      }
      if( reco_beam_incidentEnergies.size() ) reco_beam_interactingEnergy = reco_beam_incidentEnergies.back();
    }
  }

  //Uncorrected
  auto calo_NoSCE = trackUtil.GetRecoTrackCalorimetry(*thisTrack, evt,
                                                      fTrackerTag,
                                                      fCalorimetryTagNoSCE);
  found_calo = false;
  index = 0;
  for (index = 0; index < calo_NoSCE.size(); ++index) {
    if (calo_NoSCE[index].PlaneID().Plane == 2) {
      found_calo = true;
      break; 
    }
  }

  if (found_calo) {
    auto calo_dQdX = calo_NoSCE[index].dQdx();
    auto calo_dEdX = calo_NoSCE[index].dEdx();
    auto calo_range = calo_NoSCE[index].ResidualRange();
    auto TpIndices = calo_NoSCE[index].TpIndices();

    std::vector< size_t > calo_hit_indices;
    for( size_t i = 0; i < calo_dQdX.size(); ++i ){
      if (fVerbose) std::cout << i << std::endl;
      reco_beam_dQdX_NoSCE.push_back( calo_dQdX[i] );
      reco_beam_dEdX_NoSCE.push_back( calo_dEdX[i] );
      reco_beam_resRange_NoSCE.push_back( calo_range[i] );
      reco_beam_TrkPitch_NoSCE.push_back( calo_NoSCE[index].TrkPitchVec()[i] );
      calo_hit_indices.push_back( TpIndices[i] );
      const recob::Hit & theHit = (*allHits)[ TpIndices[i] ];
      reco_beam_dQ_NoSCE.push_back(theHit.Integral());
      reco_beam_calo_TPC_NoSCE.push_back(theHit.WireID().TPC);
      if (theHit.WireID().TPC == 1) {
        reco_beam_calo_wire_NoSCE.push_back( theHit.WireID().Wire );
      }
      else if (theHit.WireID().TPC == 5) {
        reco_beam_calo_wire_NoSCE.push_back( theHit.WireID().Wire + 479);
      }
      //Need other TPCs?
      else {
        reco_beam_calo_wire_NoSCE.push_back(theHit.WireID().Wire );
      }
      reco_beam_calo_wire_z_NoSCE.push_back(
          wireReadout.Wire(theHit.WireID()).GetCenter().Z());

    }

    if (fRecalibrate) {
      std::vector< float > new_dEdX = calibration_NoSCE.GetCalibratedCalorimetry(  *thisTrack, evt, fTrackerTag, fCalorimetryTagNoSCE, 2, -1.);
      for( size_t i = 0; i < new_dEdX.size(); ++i ){ reco_beam_calibrated_dEdX_NoSCE.push_back( new_dEdX[i] ); }
    }
    else {
      for (auto dedx : reco_beam_dEdX_NoSCE) {
        reco_beam_calibrated_dEdX_NoSCE.push_back(dedx);
      }
    }
    ////////////////////////////////////////////

    std::vector< calo_point > reco_beam_calo_points;
    if (reco_beam_calibrated_dEdX_NoSCE.size() &&
        reco_beam_calibrated_dEdX_NoSCE.size() == reco_beam_TrkPitch_NoSCE.size() &&
        reco_beam_calibrated_dEdX_NoSCE.size() == reco_beam_calo_wire.size()) {

      for( size_t i = 0; i < reco_beam_calibrated_dEdX_NoSCE.size(); ++i ){
        reco_beam_calo_points.push_back(
          calo_point(reco_beam_calo_wire_NoSCE[i], 0.,
                     reco_beam_TrkPitch_NoSCE[i],
                     reco_beam_dQdX_NoSCE[i], reco_beam_dEdX_NoSCE[i],
                     reco_beam_dQ_NoSCE[i], 0.,
                     reco_beam_calibrated_dEdX_NoSCE[i],
                     reco_beam_resRange_NoSCE[i], calo_hit_indices[i],
                     reco_beam_calo_wire_z_NoSCE[i], reco_beam_calo_TPC_NoSCE[i],
                     0., 0., 0., 0.));
      }

      //std::cout << "N Calo points: " << reco_beam_calo_points.size() << std::endl;
      //Sort
      std::sort( reco_beam_calo_points.begin(), reco_beam_calo_points.end(), [](calo_point a, calo_point b) {return ( a.z < b.z );} );

      //And also put these in the right order
      for( size_t i = 0; i < reco_beam_calo_points.size(); ++i ){
        calo_point thePoint = reco_beam_calo_points[i];
        reco_beam_calo_wire_NoSCE[i] = thePoint.wire;
        reco_beam_calibrated_dEdX_NoSCE[i] = thePoint.calibrated_dEdX;
        reco_beam_TrkPitch_NoSCE[i] = thePoint.pitch;
        calo_hit_indices[i] = thePoint.hit_index;
        reco_beam_calo_wire_z_NoSCE[i] = thePoint.z;
        reco_beam_calo_TPC_NoSCE[i] = thePoint.tpc;
        reco_beam_resRange_NoSCE[i] = thePoint.res_range;
        reco_beam_dQdX_NoSCE[i] = thePoint.dQdX;
        reco_beam_dEdX_NoSCE[i] = thePoint.dEdX;
      }
    }
  }
}

//info from reconstructed beam shower info
void pduneana::PDSPAnalyzer::BeamShowerInfo(const art::Event & evt, const recob::Shower* thisShower) {
  reco_beam_type = 11;
  reco_beam_trackID = thisShower->ID();
  reco_beam_startX = thisShower->ShowerStart().X();
  reco_beam_startY = thisShower->ShowerStart().Y();
  reco_beam_startZ = thisShower->ShowerStart().Z();
  reco_beam_trackDirX = thisShower->Direction().X();
  reco_beam_trackDirY = thisShower->Direction().Y();
  reco_beam_trackDirZ = thisShower->Direction().Z();

  auto allHits = evt.getValidHandle<std::vector<recob::Hit> >(fHitTag);
  auto calo = showerUtil.GetRecoShowerCalorimetry(*thisShower, evt, fShowerTag, "pandoraShowercalo");
  bool found_calo = false;
  size_t index = 0;
  for ( index = 0; index < calo.size(); ++index) {
    if (calo[index].PlaneID().Plane == 2) {
      found_calo = true;
      break; 
    }
  }
  if (found_calo) {
    auto TpIndices = calo[index].TpIndices();
    for( size_t i = 0; i < TpIndices.size(); ++i ){
      const recob::Hit & theHit = (*allHits)[ TpIndices[i] ];
      reco_beam_calo_TPC.push_back(theHit.WireID().TPC);
      if (theHit.WireID().TPC == 1) {
        reco_beam_calo_wire.push_back( theHit.WireID().Wire );
      }
      else if (theHit.WireID().TPC == 5) {
        reco_beam_calo_wire.push_back( theHit.WireID().Wire + 479);
      }
      //Need other TPCs?
      else {
        reco_beam_calo_wire.push_back(theHit.WireID().Wire );
      }

      reco_beam_calo_tick.push_back( theHit.PeakTime() );
    }
  }

  if (fVerbose) {
    std::cout << "Beam particle is shower-like" << std::endl;
    std::cout << thisShower->ShowerStart().X() << " " << thisShower->ShowerStart().Y() << " " << thisShower->ShowerStart().Z() << std::endl;
    std::cout << thisShower->Direction().X() << " " << thisShower->Direction().Y() << " " << thisShower->Direction().Z() << std::endl;
  }
}

void pduneana::PDSPAnalyzer::CheckEff(const art::Event & evt, int rightTrackID) {
  auto hitHandler = evt.getValidHandle<std::vector<recob::Hit> >("hitpdune");
  auto spHandler = evt.getValidHandle< std::vector<recob::SpacePoint> >("hitpdune");
  art::FindManyP<recob::Hit> hitFromSP(spHandler, evt, "hitpdune");
  art::FindManyP<recob::Track> trackFromHit(hitHandler, evt, "pandoraTrack");
  if (hitFromSP.size() > 0) {
    for_truncation_method = 0;
    auto const & hits = hitFromSP.at(0);
    std::vector<int> keys;
    for (auto const & hit: hits){
      auto const & tracks = trackFromHit.at(hit.key()); // const std::vector<art::Ptr<recob::Track>, std::allocator<art::Ptr<recob::Track> > >
      if (!tracks.empty()){
        keys.push_back(tracks[0].key());
      }
    }
    
    // find the most common element (mode) in keys
    int repetition = 0;
    int mode = -2;
    std::map<int,int> mmap;
    //for (std::vector<int>::iterator vi = keys.begin(); vi != keys.end(); vi++) {
    for (auto vi: keys) {
      mmap[vi]++;
      if (mmap[vi] > repetition) {
        repetition = mmap[vi];
        mode = vi;
      }
    }
    
    if (repetition > hits.size()/2.) { // if repetition of mode in keys > 1/2 size of tagged hits
      for_truncation_method = 1; // reconstructed
      if (mode == rightTrackID) {
        for_truncation_method = 2; // identified
      }
    }
  }
}

//Info from the true beam from event only in MC
void pduneana::PDSPAnalyzer::TrueBeamInfo(
    const art::Event & evt, const simb::MCParticle* true_beam_particle,
    detinfo::DetectorClocksData const& clockData,
    const sim::ParticleList & plist,
    std::map<int, std::vector<int>> & trueToPFPs,
    anab::MVAReader<recob::Hit,4> * hitResults) {
  auto pfpVec
      = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
  true_beam_endProcess = true_beam_particle->EndProcess();
  true_beam_PDG         = true_beam_particle->PdgCode();
  true_beam_mass        = true_beam_particle->Mass()*1.e3;
  true_beam_ID          = true_beam_particle->TrackId();
  std::vector<int> primary_hierarchy = PrimaryHierarchy(true_beam_ID, plist);
  true_beam_hierarchy.insert(true_beam_hierarchy.end(),
                                primary_hierarchy.begin(),
                                primary_hierarchy.end());
  true_beam_hierarchy_size = true_beam_hierarchy.size();
  reco_beam_true_byHits_in_primary_hierarchy
      = (std::find(true_beam_hierarchy.begin(), true_beam_hierarchy.end(),
                   reco_beam_true_byHits_ID)
         != true_beam_hierarchy.end());
  //std::cout << "Primary Hierarchy has " << primary_hierarchy.size() << std::endl;
  //for (auto & i : primary_hierarchy) {
  //  std::cout << i << std::endl;
  //}
  true_beam_endX = true_beam_particle->EndX();
  true_beam_endY = true_beam_particle->EndY();
  true_beam_endZ = true_beam_particle->EndZ();
  true_beam_startX     = true_beam_particle->Position(0).X();
  true_beam_startY     = true_beam_particle->Position(0).Y();
  true_beam_startZ     = true_beam_particle->Position(0).Z();

  auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto offset = sce->GetPosOffsets(
      {true_beam_endX, true_beam_endY, true_beam_endZ});
  true_beam_endX_SCE = true_beam_endX - offset.X();
  true_beam_endY_SCE = true_beam_endY + offset.Y();
  true_beam_endZ_SCE = true_beam_endZ + offset.Z();


  true_beam_startPx    = true_beam_particle->Px();
  true_beam_startPy    = true_beam_particle->Py();
  true_beam_startPz    = true_beam_particle->Pz();
  true_beam_startP     = true_beam_particle->P();

  size_t true_np = true_beam_particle->NumberTrajectoryPoints();

  true_beam_endPx    = true_beam_particle->Px(true_np-2);
  true_beam_endPy    = true_beam_particle->Py(true_np-2);
  true_beam_endPz    = true_beam_particle->Pz(true_np-2);
  true_beam_endP     = true_beam_particle->P(true_np-2);
  true_beam_endP2     = true_beam_particle->P(true_np-1);

  true_beam_last_len
      = sqrt(std::pow((true_beam_particle->Position(true_np-2).X() -
                       true_beam_particle->Position(true_np-1).X()), 2) +
             std::pow((true_beam_particle->Position(true_np-2).Y() -
                       true_beam_particle->Position(true_np-1).Y()), 2) +
             std::pow((true_beam_particle->Position(true_np-2).Z() -
                       true_beam_particle->Position(true_np-1).Z()), 2));

  true_beam_startDirX  = true_beam_startPx / true_beam_startP;
  true_beam_startDirY  = true_beam_startPy / true_beam_startP;
  true_beam_startDirZ  = true_beam_startPz / true_beam_startP;

  true_beam_nHits = truthUtil.GetMCParticleHits( clockData, *true_beam_particle, evt, fHitTag ).size();

  //Checking all reconstructed objects to see what the true beam particle
  //contributed to
  true_beam_reco_byHits_PFP_ID.push_back( std::vector< int >() );
  true_beam_reco_byHits_PFP_nHits.push_back( std::vector< int >() );
  true_beam_reco_byHits_allTrack_ID.push_back( std::vector< int >() );
  if( fTrueToReco ){
    for( size_t i = 0; i < trueToPFPs[ true_beam_ID ].size(); ++i ){
      true_beam_reco_byHits_PFP_ID.back().push_back( trueToPFPs[ true_beam_ID ][i] );

      const recob::PFParticle * thePFP = &(pfpVec->at( trueToPFPs[ true_beam_ID ][i] ));
      true_beam_reco_byHits_PFP_nHits.back().push_back(
        pfpUtil.GetPFParticleHits_Ptrs( *thePFP, evt, fPFParticleTag ).size()
      );

      const recob::Track* pandora2Track = 0x0;

      try{
        pandora2Track = pfpUtil.GetPFParticleTrack( *thePFP, evt, fPFParticleTag, "pandora2Track" );
      }
      catch( const cet::exception &e ){
        MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Track object not found, moving on" << std::endl;
      }

      if( pandora2Track ){
        true_beam_reco_byHits_allTrack_ID.back().push_back( pandora2Track->ID() );
      }
      else{
        true_beam_reco_byHits_allTrack_ID.back().push_back( -1 );
      }
    }
  }

  //Truth thin slice info
  //Go through the true processes within the MCTrajectory
  const simb::MCTrajectory & true_beam_trajectory = true_beam_particle->Trajectory();
  auto true_beam_proc_map = true_beam_trajectory.TrajectoryProcesses();
  if (fVerbose) std::cout << "Processes: " << std::endl;

  for( auto itProc = true_beam_proc_map.begin(); itProc != true_beam_proc_map.end(); ++itProc ){
    int index = itProc->first;
    std::string process = true_beam_trajectory.KeyToProcess(itProc->second);
    if (fVerbose) std::cout << index << " " << process << std::endl;

    true_beam_processes.push_back( process );

    if( process == "hadElastic" ){

      ++true_beam_nElasticScatters;

      double process_X = true_beam_trajectory.X( index );
      double process_Y = true_beam_trajectory.Y( index );
      double process_Z = true_beam_trajectory.Z( index );

      double PX      = true_beam_trajectory.Px( index );
      double next_PX = true_beam_trajectory.Px( index + 1 );
      double PY      = true_beam_trajectory.Py( index );
      double next_PY = true_beam_trajectory.Py( index + 1 );
      double PZ      = true_beam_trajectory.Pz( index );
      double next_PZ = true_beam_trajectory.Pz( index + 1 );

      double total_P = sqrt( PX*PX + PY*PY + PZ*PZ );
      double total_next_P = sqrt( next_PX*next_PX + next_PY*next_PY + next_PZ*next_PZ );

      //Get the angle between the direction of this step and the next
      true_beam_elastic_costheta.push_back(
        ( ( PX * next_PX ) + ( PY * next_PY ) + ( PZ * next_PZ ) ) / ( total_P * total_next_P )
      );

      double total_E = sqrt(total_P*total_P*1.e6 +
                            true_beam_mass*true_beam_mass);
      double total_next_E = sqrt(total_next_P*total_next_P*1.e6 +
                                 true_beam_mass*true_beam_mass);

      true_beam_elastic_X.push_back( process_X );
      true_beam_elastic_Y.push_back( process_Y );
      true_beam_elastic_Z.push_back( process_Z );

      true_beam_elastic_deltaE.push_back(total_E - total_next_E);

      std::vector<const sim::IDE *> ides_between_points =
          truthUtil.GetSimIDEsBetweenPoints(
              *true_beam_particle, true_beam_trajectory.Position(index),
              true_beam_trajectory.Position(index +1));

      double total_edep = 0.;
      for (size_t i = 0; i < ides_between_points.size(); ++i) {
        total_edep += ides_between_points[i]->energy;
      }
      true_beam_elastic_IDE_edep.push_back(total_edep);

    }
  }
  if( true_beam_endProcess.find( "Inelastic" ) == std::string::npos ){
    true_beam_processes.push_back( true_beam_endProcess );
  }


  //Looking at info from IDEs -- simulated energy depositions
  //To do: might remove because the IDEs are 'eaten up' due to various
  //changes in the simulation. So their original intent
  //(determining slices for thin slice method) was made useless
  if (fVerbose) std::cout << "Looking at IDEs" << std::endl;
  auto view2_IDEs = bt_serv->TrackIdToSimIDEs_Ps( true_beam_ID, geo::View_t(2) );
  if (fVerbose) std::cout << "N view2 IDEs: " << view2_IDEs.size() << std::endl;

  true_beam_IDE_totalDep = 0.;
  for (size_t i = 0; i < view2_IDEs.size(); ++i) {
    const sim::IDE * ide = view2_IDEs[i];

    true_beam_IDE_totalDep += ide->energy;
  }
  //Sort based on ide z-position
  std::sort( view2_IDEs.begin(), view2_IDEs.end(), sort_IDEs );

  //This is an attempt to remove IDEs from things like delta-rays
  //that have a large gap in z to the previous
  size_t remove_index = 0;
  bool   do_remove = false;
  if( view2_IDEs.size() ){
    for( size_t i = 1; i < view2_IDEs.size()-1; ++i ){
      const sim::IDE * prev_IDE = view2_IDEs[i-1];
      const sim::IDE * this_IDE = view2_IDEs[i];

      if (this_IDE->trackID < 0.) {
        em_energy += this_IDE->energy;
      }


      if( this_IDE->trackID < 0 && ( this_IDE->z - prev_IDE->z ) > 5 ){
        remove_index = i;
        do_remove = true;
        break;   
      }
    }
  }

  if( do_remove ){
    view2_IDEs.erase( view2_IDEs.begin() + remove_index, view2_IDEs.end() );
  }

  double init_KE = sqrt(1.e6 * true_beam_startP*true_beam_startP +
                        true_beam_mass*true_beam_mass) - true_beam_mass;

  //slice up the view2_IDEs up by the wire pitch
  auto sliced_ides = slice_IDEs( view2_IDEs, fZ0, fPitch, true_beam_endZ);
  //Get the momentum at the start of the slices.
  //Get the first slice
  if (fVerbose) std::cout << "size: " << sliced_ides.size() << std::endl;
  if (sliced_ides.size()) {
    auto first_slice = sliced_ides.begin();
    if (fVerbose) std::cout << "Got first slice" << std::endl;

    //Check it has any IDEs
    auto theIDEs = first_slice->second;
    if (fVerbose) std::cout << "Got ides" << std::endl;
    if (theIDEs.size()) {
      //Get the first ide z position
      double ide_z = theIDEs[0]->z;

      //Go through the trajectory position
      //and check for the position that comes immediately before the
      //first ide
      for (size_t i = 1; i < true_beam_trajectory.size(); ++i) {
        double z0 = true_beam_trajectory.Z(i-1);
        double z1 = true_beam_trajectory.Z(i);

        double x0 = true_beam_trajectory.X(i-1);
        double y0 = true_beam_trajectory.Y(i-1);
        double x1 = true_beam_trajectory.X(i);
        double y1 = true_beam_trajectory.Y(i);
        geo::Point_t test_point{x0, y0, z0};
        geo::Point_t test_point_1{x1, y1, z1};
        //const TGeoMaterial * mat = geom->Material(test_point);
        if (z0 < ide_z && z1 > ide_z) {
          //const TGeoMaterial * mat1 = geom->Material(test_point_1);
          init_KE = 1.e3 * true_beam_trajectory.E(i-1) - true_beam_mass;
          if (fVerbose) {
            std::cout << "Found matching position" << z0 << " " << ide_z <<
                         " " << z1 << std::endl;
            std::cout << "init KE: " << init_KE << " " <<
                         (1.e3*true_beam_trajectory.E(i) - true_beam_mass) << std::endl;
          }
          break;
        }
      }
    }
  }

  //Go through the sliced up IDEs to create the thin targets 
  true_beam_incidentEnergies.push_back( init_KE );
  for( auto it = sliced_ides.begin(); it != sliced_ides.end(); ++it ){

    auto theIDEs = it->second;

    true_beam_slices.push_back( it->first );

    double deltaE = 0.;
    for( size_t i = 0; i < theIDEs.size(); ++i ){
      deltaE += theIDEs[i]->energy;
    }

    true_beam_slices_deltaE.push_back( deltaE );
    true_beam_incidentEnergies.push_back( true_beam_incidentEnergies.back() - deltaE );
  }
  true_beam_incidentEnergies.pop_back();
  if( true_beam_incidentEnergies.size() ) true_beam_interactingEnergy = true_beam_incidentEnergies.back();

  //Save the trajectory points
  art::ServiceHandle <geo::Geometry> geo_serv;
  for (size_t i = 0; i < true_beam_trajectory.size(); ++i) {
    true_beam_traj_X.push_back(true_beam_trajectory.X(i));
    true_beam_traj_Y.push_back(true_beam_trajectory.Y(i));
    true_beam_traj_Z.push_back(true_beam_trajectory.Z(i));
    true_beam_traj_Px.push_back(true_beam_trajectory.Px(i));
    true_beam_traj_Py.push_back(true_beam_trajectory.Py(i));
    true_beam_traj_Pz.push_back(true_beam_trajectory.Pz(i));
    
    true_beam_traj_KE.push_back(true_beam_trajectory.E(i)*1.e3 - true_beam_mass);

    auto offset = sce->GetPosOffsets(
        {true_beam_trajectory.X(i), true_beam_trajectory.Y(i),
         true_beam_trajectory.Z(i)});
    true_beam_traj_X_SCE.push_back(true_beam_trajectory.X(i) - offset.X());
    true_beam_traj_Y_SCE.push_back(true_beam_trajectory.Y(i) + offset.Y());
    true_beam_traj_Z_SCE.push_back(true_beam_trajectory.Z(i) + offset.Z());


    geo::Point_t test_point{true_beam_trajectory.X(i),
                            true_beam_trajectory.Y(i),
                            true_beam_trajectory.Z(i)};
    const TGeoMaterial * test_material = geo_serv->Material(test_point);
    if (fVerbose) std::cout << test_material->GetName() << std::endl;
    if (!strcmp(test_material->GetName(), "ALUMINUM_Al")) true_beam_is_scraper = true;

  }
  true_beam_len = 0;
  for (size_t i = 1; i < true_beam_trajectory.size(); ++i) {
    true_beam_len += sqrt( pow( true_beam_traj_X.at(i)-true_beam_traj_X.at(i-1), 2)
                       + pow( true_beam_traj_Y.at(i)-true_beam_traj_Y.at(i-1), 2)
                       + pow( true_beam_traj_Z.at(i)-true_beam_traj_Z.at(i-1), 2));
  }

  //Look through the daughters
  for( int i = 0; i < true_beam_particle->NumberDaughters(); ++i ){
    int daughterID = true_beam_particle->Daughter(i);

    if (fVerbose) std::cout << "Daughter " << i << " ID: " << daughterID << std::endl;
    auto part = plist[ daughterID ];
    int pid = part->PdgCode();

    std::string process = part->Process();

    if (process == "muIoni" || process == "hIoni")
      continue;

    true_beam_daughter_PDG.push_back(pid);
    true_beam_daughter_ID.push_back( part->TrackId() );

    true_beam_daughter_len.push_back( part->Trajectory().TotalLength() );

    true_beam_daughter_startX.push_back( part->Position(0).X() );
    true_beam_daughter_startY.push_back( part->Position(0).Y() );
    true_beam_daughter_startZ.push_back( part->Position(0).Z() );

    true_beam_daughter_endX.push_back( part->EndX() );
    true_beam_daughter_endY.push_back( part->EndY() );
    true_beam_daughter_endZ.push_back( part->EndZ() );

    true_beam_daughter_startPx.push_back( part->Px() );
    true_beam_daughter_startPy.push_back( part->Py() );
    true_beam_daughter_startPz.push_back( part->Pz() );
    true_beam_daughter_startP.push_back( part->P() );

    true_beam_daughter_Process.push_back( part->Process() );
    true_beam_daughter_endProcess.push_back( part->EndProcess() );

    if (fVerbose) {
      std::cout << "Process: " << part->Process() << std::endl;
      std::cout << "PID: " << pid << std::endl;
      std::cout << "Start: " << part->Position(0).X() << " " << part->Position(0).Y() << " " << part->Position(0).Z() << std::endl;
      std::cout << "End: " << part->EndPosition().X() << " " << part->EndPosition().Y() << " " << part->EndPosition().Z() << std::endl;
      std::cout << "Len: " << part->Trajectory().TotalLength() << std::endl;
    }

    if( part->Process().find( "Inelastic" ) != std::string::npos ){
      if( pid == 211  ) ++true_daughter_nPiPlus;
      if( pid == -211 ) ++true_daughter_nPiMinus;
      if( pid == 111  ) ++true_daughter_nPi0;
      if( pid == 2212 ) ++true_daughter_nProton;
      if( pid == 2112 ) ++true_daughter_nNeutron;
      if( pid > 2212  ) ++true_daughter_nNucleus;
    }

    //Look for the gammas coming out of the pi0s
    if( pid == 111 ){
      //std::cout << "Found pi0. Looking at true daughters" << std::endl;
      for( int j = 0; j < part->NumberDaughters(); ++j ){
        int pi0_decay_daughter_ID = part->Daughter(j);
        auto pi0_decay_part = plist[ pi0_decay_daughter_ID ];
        true_beam_Pi0_decay_PDG.push_back( pi0_decay_part->PdgCode() );
        true_beam_Pi0_decay_ID.push_back( pi0_decay_part->TrackId() );
        true_beam_Pi0_decay_startP.push_back( pi0_decay_part->P() );
        true_beam_Pi0_decay_startPx.push_back( pi0_decay_part->Px() );
        true_beam_Pi0_decay_startPy.push_back( pi0_decay_part->Py() );
        true_beam_Pi0_decay_startPz.push_back( pi0_decay_part->Pz() );
        true_beam_Pi0_decay_startX.push_back( pi0_decay_part->Position(0).X() );
        true_beam_Pi0_decay_startY.push_back( pi0_decay_part->Position(0).Y() );
        true_beam_Pi0_decay_startZ.push_back( pi0_decay_part->Position(0).Z() );
        true_beam_Pi0_decay_parID.push_back( pi0_decay_part->Mother() );

        true_beam_Pi0_decay_len.push_back( pi0_decay_part->Trajectory().TotalLength() );
        true_beam_Pi0_decay_nHits.push_back( truthUtil.GetMCParticleHits( clockData, *pi0_decay_part, evt, fHitTag ).size() );

        true_beam_Pi0_decay_reco_byHits_PFP_ID.push_back( std::vector<int>() );
        true_beam_Pi0_decay_reco_byHits_PFP_nHits.push_back( std::vector<int>() );
        true_beam_Pi0_decay_reco_byHits_PFP_trackScore.push_back( std::vector<double>() );

        true_beam_Pi0_decay_reco_byHits_allTrack_ID.push_back( std::vector<int>() );
        true_beam_Pi0_decay_reco_byHits_allTrack_startX.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allTrack_startY.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allTrack_startZ.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allTrack_len.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allTrack_endX.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allTrack_endY.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allTrack_endZ.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allShower_ID.push_back( std::vector<int>() );
        true_beam_Pi0_decay_reco_byHits_allShower_startX.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allShower_startY.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allShower_startZ.push_back( std::vector<double>() );
        true_beam_Pi0_decay_reco_byHits_allShower_len.push_back( std::vector<double>() );
        if( fTrueToReco ){
          for( size_t k = 0; k < trueToPFPs[ pi0_decay_part->TrackId() ].size(); ++k ){
            true_beam_Pi0_decay_reco_byHits_PFP_ID.back().push_back( trueToPFPs[ pi0_decay_part->TrackId() ][k] );

            const recob::PFParticle * thePFP = &(pfpVec->at( trueToPFPs[ pi0_decay_part->TrackId() ][k] ));
            true_beam_Pi0_decay_reco_byHits_PFP_nHits.back().push_back(
              pfpUtil.GetPFParticleHits_Ptrs( *thePFP, evt, fPFParticleTag ).size()
            );

            if (!fSkipMVA) {
              cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle( *thePFP, evt, *hitResults, pfpUtil, fPFParticleTag );
              double track_score = ((theCNNResults.nHits > 0) ? 
                                    (theCNNResults.track / theCNNResults.nHits) :
                                    -999.);
              true_beam_Pi0_decay_reco_byHits_PFP_trackScore.back().push_back(track_score);
            }
            else {
              true_beam_Pi0_decay_reco_byHits_PFP_trackScore.back().push_back(-999.);
            }

            const recob::Track* pandora2Track = 0x0;
            try{
              pandora2Track = pfpUtil.GetPFParticleTrack( *thePFP, evt, fPFParticleTag, "pandora2Track" );
            }
            catch( const cet::exception &e ){
              MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Track object not found, moving on" << std::endl;
            }

            if( pandora2Track ){
              true_beam_Pi0_decay_reco_byHits_allTrack_ID.back().push_back( pandora2Track->ID() );
              true_beam_Pi0_decay_reco_byHits_allTrack_startX.back().push_back( pandora2Track->Trajectory().Start().X() );
              true_beam_Pi0_decay_reco_byHits_allTrack_startY.back().push_back( pandora2Track->Trajectory().Start().Y() );
              true_beam_Pi0_decay_reco_byHits_allTrack_startZ.back().push_back( pandora2Track->Trajectory().Start().Z() );
              true_beam_Pi0_decay_reco_byHits_allTrack_endX.back().push_back( pandora2Track->Trajectory().End().X() );
              true_beam_Pi0_decay_reco_byHits_allTrack_endY.back().push_back( pandora2Track->Trajectory().End().Y() );
              true_beam_Pi0_decay_reco_byHits_allTrack_endZ.back().push_back( pandora2Track->Trajectory().End().Z() );
              true_beam_Pi0_decay_reco_byHits_allTrack_len.back().push_back( pandora2Track->Length() );
     
            }
            else{
              true_beam_Pi0_decay_reco_byHits_allTrack_ID.back().push_back( -999 );
              true_beam_Pi0_decay_reco_byHits_allTrack_startX.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_startY.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_startZ.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_endX.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_endY.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_endZ.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allTrack_len.back().push_back( -999. );
            }

            const recob::Shower* pandora2Shower = 0x0;
            try{
              pandora2Shower = pfpUtil.GetPFParticleShower( *thePFP, evt, fPFParticleTag, "pandora2Shower" );
            }
            catch( const cet::exception &e ){
              MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Shower object not found, moving on" << std::endl;
            }

            if( pandora2Shower ){
              true_beam_Pi0_decay_reco_byHits_allShower_ID.back().push_back( pandora2Shower->ID() );
              true_beam_Pi0_decay_reco_byHits_allShower_startX.back().push_back( pandora2Shower->ShowerStart().X() );
              true_beam_Pi0_decay_reco_byHits_allShower_startY.back().push_back( pandora2Shower->ShowerStart().Y() );
              true_beam_Pi0_decay_reco_byHits_allShower_startZ.back().push_back( pandora2Shower->ShowerStart().Z() );
              true_beam_Pi0_decay_reco_byHits_allShower_len.back().push_back( pandora2Shower->Length() );
            }
            else{
              true_beam_Pi0_decay_reco_byHits_allShower_ID.back().push_back( -999 );
              true_beam_Pi0_decay_reco_byHits_allShower_startX.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allShower_startY.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allShower_startZ.back().push_back( -999. );
              true_beam_Pi0_decay_reco_byHits_allShower_len.back().push_back( -999. );
            }

          }
        }
      }
    }

    //Look another level down
    for( int j = 0; j < part->NumberDaughters(); ++j ){
      int grand_daughter_ID = part->Daughter(j);
      auto grand_daughter_part = plist[ grand_daughter_ID ];
      true_beam_grand_daughter_PDG.push_back( grand_daughter_part->PdgCode() );
      true_beam_grand_daughter_ID.push_back(  grand_daughter_part->TrackId() );
      true_beam_grand_daughter_parID.push_back(  part->TrackId() );
      true_beam_grand_daughter_nHits.push_back( truthUtil.GetMCParticleHits( clockData, *grand_daughter_part, evt, fHitTag ).size() );
      true_beam_grand_daughter_Process.push_back( grand_daughter_part->Process() );
      true_beam_grand_daughter_endProcess.push_back( grand_daughter_part->EndProcess() );
    }



    //Same dumb cross checking of true particles to reconstructed objects
    //To do: might remove
    true_beam_daughter_reco_byHits_PFP_ID.push_back( std::vector<int>() );
    true_beam_daughter_reco_byHits_PFP_nHits.push_back( std::vector<int>() );
    true_beam_daughter_reco_byHits_PFP_trackScore.push_back( std::vector<double>() );

    true_beam_daughter_reco_byHits_allTrack_ID.push_back( std::vector<int>() );
    true_beam_daughter_reco_byHits_allTrack_startX.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allTrack_startY.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allTrack_startZ.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allTrack_len.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allTrack_endX.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allTrack_endY.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allTrack_endZ.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allShower_ID.push_back( std::vector<int>() );
    true_beam_daughter_reco_byHits_allShower_startX.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allShower_startY.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allShower_startZ.push_back( std::vector<double>() );
    true_beam_daughter_reco_byHits_allShower_len.push_back( std::vector<double>() );
    if( fTrueToReco ){
      for( size_t j = 0; j < trueToPFPs[ part->TrackId() ].size(); ++j ){
        true_beam_daughter_reco_byHits_PFP_ID.back().push_back( trueToPFPs[ part->TrackId() ][j] );

        const recob::PFParticle * thePFP = &(pfpVec->at( trueToPFPs[ part->TrackId() ][j] ));
        true_beam_daughter_reco_byHits_PFP_nHits.back().push_back(
          pfpUtil.GetPFParticleHits_Ptrs( *thePFP, evt, fPFParticleTag ).size()
        );

        if (!fSkipMVA) {
          cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle( *thePFP, evt, *hitResults, pfpUtil, fPFParticleTag );
          double track_score = ((theCNNResults.nHits > 0) ? 
                                (theCNNResults.track / theCNNResults.nHits) :
                                -999.);
          true_beam_daughter_reco_byHits_PFP_trackScore.back().push_back(track_score);
        }
        else {
          true_beam_daughter_reco_byHits_PFP_trackScore.back().push_back(-999.);
        }

        //cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle( *thePFP, evt, *hitResults, pfpUtil, fPFParticleTag );
        //true_beam_daughter_reco_byHits_PFP_trackScore.back().push_back( ( ( theCNNResults.nHits > 0 ) ? ( theCNNResults.track / theCNNResults.nHits ) : -999. ) );

        const recob::Track* pandora2Track = 0x0;
        try{
           pandora2Track = pfpUtil.GetPFParticleTrack( *thePFP, evt, fPFParticleTag, "pandora2Track" );
        }
        catch( const cet::exception &e ){
          MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Track object not found, moving on" << std::endl;
        }

        if( pandora2Track ){
          true_beam_daughter_reco_byHits_allTrack_ID.back().push_back( pandora2Track->ID() );
          true_beam_daughter_reco_byHits_allTrack_startX.back().push_back( pandora2Track->Trajectory().Start().X() );
          true_beam_daughter_reco_byHits_allTrack_startY.back().push_back( pandora2Track->Trajectory().Start().Y() );
          true_beam_daughter_reco_byHits_allTrack_startZ.back().push_back( pandora2Track->Trajectory().Start().Z() );
          true_beam_daughter_reco_byHits_allTrack_endX.back().push_back( pandora2Track->Trajectory().End().X() );
          true_beam_daughter_reco_byHits_allTrack_endY.back().push_back( pandora2Track->Trajectory().End().Y() );
          true_beam_daughter_reco_byHits_allTrack_endZ.back().push_back( pandora2Track->Trajectory().End().Z() );
          true_beam_daughter_reco_byHits_allTrack_len.back().push_back( pandora2Track->Length() );
 
        }
        else{
          true_beam_daughter_reco_byHits_allTrack_ID.back().push_back( -999 );
          true_beam_daughter_reco_byHits_allTrack_startX.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_startY.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_startZ.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_endX.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_endY.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_endZ.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allTrack_len.back().push_back( -999. );
        }

        const recob::Shower* pandora2Shower = 0x0;
        try{
          pandora2Shower = pfpUtil.GetPFParticleShower( *thePFP, evt, fPFParticleTag, "pandora2Shower" );
        }
        catch( const cet::exception &e ){
          MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Shower object not found, moving on" << std::endl;
        }
 
        if( pandora2Shower ){
          true_beam_daughter_reco_byHits_allShower_ID.back().push_back( pandora2Shower->ID() );
          true_beam_daughter_reco_byHits_allShower_startX.back().push_back( pandora2Shower->ShowerStart().X() );
          true_beam_daughter_reco_byHits_allShower_startY.back().push_back( pandora2Shower->ShowerStart().Y() );
          true_beam_daughter_reco_byHits_allShower_startZ.back().push_back( pandora2Shower->ShowerStart().Z() );
          true_beam_daughter_reco_byHits_allShower_len.back().push_back( pandora2Shower->Length() );
 
        }
        else{
          true_beam_daughter_reco_byHits_allShower_ID.back().push_back( -999 );
          true_beam_daughter_reco_byHits_allShower_startX.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allShower_startY.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allShower_startZ.back().push_back( -999. );
          true_beam_daughter_reco_byHits_allShower_len.back().push_back( -999. );
        }

      }
    }
    true_beam_daughter_nHits.push_back( truthUtil.GetMCParticleHits( clockData, *part, evt, fHitTag ).size() );

  }
}

void pduneana::PDSPAnalyzer::BeamInstInfo(const art::Event & evt) {
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;

  //Get beam instrumentation info
  //Works for both MC and data, can probably merge these (To do)
  if (evt.isRealData()) {
    auto beamHandle
        = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(
            fBeamModuleLabel);

    if (beamHandle.isValid()) {
      art::fill_ptr_vector(beamVec, beamHandle);
    }
  }
  else {
    try {
      auto beamHandle
          = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(
              "generator");

      if( beamHandle.isValid()){
        art::fill_ptr_vector(beamVec, beamHandle);
      }
      fMCHasBI = true;
    }
    catch (const cet::exception &e) {
      MF_LOG_WARNING("PDSPAnalyzer") << "BeamEvent generator object not " <<
                                        "found, moving on" << std::endl;
      fMCHasBI = false;
      return;
    }
  }

  //High level info about the beam trigger
  const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0)); 
  beam_inst_trigger = beamEvent.GetTimingTrigger();
  if (evt.isRealData() && !fBeamlineUtils.IsGoodBeamlineTrigger(evt)) {
    std::cout << "Matched? " << beamEvent.CheckIsMatched() << std::endl;
    std::vector< double > momenta = beamEvent.GetRecoBeamMomenta();
    std::cout << momenta.size() << std::endl;
    MF_LOG_WARNING("PDSPAnalyzer") << "Failed beam quality check" << std::endl;
    beam_inst_valid = false;
    return;
  }

  //Number of tracks reconstructed in the final four fiber monitors
  int nTracks = beamEvent.GetBeamTracks().size();

  //Momentum reconstructed from spectrometer
  std::vector< double > momenta = beamEvent.GetRecoBeamMomenta();
  int nMomenta = momenta.size();

  if (fVerbose) {
    std::cout << "Got beam event" << std::endl;
    std::cout << "Got " << nTracks << " Tracks" << std::endl;
    std::cout << "Got " << nMomenta << " Momenta" << std::endl;
  }

  if( nMomenta > 0 ){
    beam_inst_P = momenta[0];
    if (fVerbose) std::cout << "reco P " << beam_inst_P << std::endl;

    if (!evt.isRealData()) {
      beam_inst_P *= fBeamInstPFix;
    }
  }

  //Time of flight + channel combinations
  const std::vector<double> the_tofs = beamEvent.GetTOFs();
  const std::vector<int> the_chans = beamEvent.GetTOFChans();
  for (size_t iTOF = 0; iTOF < the_tofs.size(); ++iTOF) {
    beam_inst_TOF.push_back(the_tofs[iTOF]);
    beam_inst_TOF_Chan.push_back(the_chans[iTOF]);
  }

  beam_inst_C0 = beamEvent.GetCKov0Status();
  beam_inst_C1 = beamEvent.GetCKov1Status();
  beam_inst_C0_pressure = beamEvent.GetCKov0Pressure();
  beam_inst_C1_pressure = beamEvent.GetCKov1Pressure();

  //position from the projected track into the TPC -- just using the first
  //index
  if( nTracks > 0 ){
    beam_inst_X = beamEvent.GetBeamTracks()[0].Trajectory().End().X();
    beam_inst_Y = beamEvent.GetBeamTracks()[0].Trajectory().End().Y();
    beam_inst_Z = beamEvent.GetBeamTracks()[0].Trajectory().End().Z();

    beam_inst_dirX = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().X();
    beam_inst_dirY = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().Y();
    beam_inst_dirZ = beamEvent.GetBeamTracks()[0].Trajectory().EndDirection().Z();
  }

  beam_inst_nTracks = nTracks;
  beam_inst_nMomenta = nMomenta;

  //beamline PID 
  if (evt.isRealData()) {
    std::vector< int > pdg_cands = fBeamlineUtils.GetPID( beamEvent, fBeamPIDMomentum );
    beam_inst_PDG_candidates.insert( beam_inst_PDG_candidates.end(), pdg_cands.begin(), pdg_cands.end() );
  }

  //number of fibers struck in the spectrometer monitors
  beam_inst_nFibersP1 = beamEvent.GetActiveFibers( "XBPF022697" ).size();
  beam_inst_nFibersP2 = beamEvent.GetActiveFibers( "XBPF022701" ).size();
  beam_inst_nFibersP3 = beamEvent.GetActiveFibers( "XBPF022702" ).size();
}

//Reconstructed info from all of the daughter PFP objects associated to the beam particle
void pduneana::PDSPAnalyzer::DaughterPFPInfo(
    const art::Event & evt, const recob::PFParticle* particle,
    detinfo::DetectorClocksData const& clockData,
    anab::MVAReader<recob::Hit,4> * hitResults) {

  //Get all PFParticles in the event
  auto pfpVec
      = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);

  auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();
  auto const detProp
      = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(
          evt, clockData);

  //Util for momentum by range calculation
  trkf::TrackMomentumCalculator track_p_calc;

  //Iterate over all daughters
  for (size_t daughterID : particle->Daughters()) {
    const recob::PFParticle * daughterPFP = &(pfpVec->at( daughterID ));
    reco_daughter_PFP_ID.push_back( daughterID );

    //Get hits info from daughter pfp
    const std::vector< art::Ptr< recob::Hit > > daughterPFP_hits = pfpUtil.GetPFParticleHits_Ptrs( *daughterPFP, evt, fPFParticleTag );
    if (fVerbose) std::cout << "Got " << daughterPFP_hits.size() << " hits from daughter " << daughterID << std::endl;

    reco_daughter_PFP_nHits.push_back( daughterPFP_hits.size() );
    size_t nHits_coll = 0;
    for (size_t i = 0; i < daughterPFP_hits.size(); ++i) {
      if (daughterPFP_hits[i]->View() == 2) {
        ++nHits_coll;
      }
    }
    reco_daughter_PFP_nHits_collection.push_back(nHits_coll);

    //Get CNN output for the daughter particle
    if (!fSkipMVA) {
      cnnOutput2D theCNNResults = GetCNNOutputFromPFParticle(
          *daughterPFP, evt, *hitResults, pfpUtil, fPFParticleTag);
      double track_score = (theCNNResults.nHits > 0 ?
                            (theCNNResults.track / theCNNResults.nHits) :
                            -999.);
      double em_score = (theCNNResults.nHits > 0 ?
                            (theCNNResults.em / theCNNResults.nHits) :
                            -999.);
      double michel_score = (theCNNResults.nHits > 0 ?
                            (theCNNResults.michel / theCNNResults.nHits) :
                            -999.);
      reco_daughter_PFP_trackScore.push_back(track_score);
      reco_daughter_PFP_emScore.push_back(em_score);
      reco_daughter_PFP_michelScore.push_back(michel_score);

      cnnOutput2D cnn_collection = GetCNNOutputFromPFParticleFromPlane(
          *daughterPFP, evt, *hitResults, pfpUtil, fPFParticleTag, 2);
      double track_score_collection = (cnn_collection.nHits > 0 ?
                            (cnn_collection.track / cnn_collection.nHits) :
                            -999.);
      double em_score_collection = (cnn_collection.nHits > 0 ?
                            (cnn_collection.em / cnn_collection.nHits) :
                            -999.);
      double michel_score_collection = (cnn_collection.nHits > 0 ?
                            (cnn_collection.michel / cnn_collection.nHits) :
                            -999.);
      reco_daughter_PFP_trackScore_collection.push_back(track_score_collection);
      reco_daughter_PFP_emScore_collection.push_back(em_score_collection);
      reco_daughter_PFP_michelScore_collection.push_back(michel_score_collection);

      reco_daughter_PFP_trackScore_weight_by_charge.push_back(theCNNResults.track_weight_by_charge);
      reco_daughter_PFP_emScore_weight_by_charge.push_back(theCNNResults.em_weight_by_charge);
      reco_daughter_PFP_michelScore_weight_by_charge.push_back(theCNNResults.michel_weight_by_charge);

      reco_daughter_PFP_trackScore_collection_weight_by_charge.push_back(cnn_collection.track_weight_by_charge);
      reco_daughter_PFP_emScore_collection_weight_by_charge.push_back(cnn_collection.em_weight_by_charge);
      reco_daughter_PFP_michelScore_collection_weight_by_charge.push_back(cnn_collection.michel_weight_by_charge);

    }
    else{
      reco_daughter_PFP_trackScore.push_back( -999. );
      reco_daughter_PFP_emScore.push_back( -999. );
      reco_daughter_PFP_michelScore.push_back( -999. );
      reco_daughter_PFP_trackScore_collection.push_back( -999. );
      reco_daughter_PFP_emScore_collection.push_back( -999. );
      reco_daughter_PFP_michelScore_collection.push_back( -999. );

      reco_daughter_PFP_trackScore_weight_by_charge.push_back( -999. );
      reco_daughter_PFP_emScore_weight_by_charge.push_back( -999. );
      reco_daughter_PFP_michelScore_weight_by_charge.push_back( -999. );
      reco_daughter_PFP_trackScore_collection_weight_by_charge.push_back( -999. );
      reco_daughter_PFP_emScore_collection_weight_by_charge.push_back( -999. );
      reco_daughter_PFP_michelScore_collection_weight_by_charge.push_back( -999. );
    }



    //Matching by hits to get true daughter info
    if( !evt.isRealData() ){
      DaughterMatchInfo(evt, daughterPFP, clockData);
    }

    if (fSaveHits) {
      reco_daughter_spacePts_X.push_back(std::vector<double>());
      reco_daughter_spacePts_Y.push_back(std::vector<double>());
      reco_daughter_spacePts_Z.push_back(std::vector<double>());

      // SparseNet features      
      sparsenet_features_charge.push_back(std::vector<double>());
      sparsenet_features_angle.push_back(std::vector<double>());
      sparsenet_features_dot_product.push_back(std::vector<double>());
      sparsenet_features_neighboring_nodes_3.push_back(std::vector<double>());
      sparsenet_features_neighboring_nodes_10.push_back(std::vector<double>());
      sparsenet_features_neighboring_nodes_30.push_back(std::vector<double>());
      sparsenet_features_charge_distance_3.push_back(std::vector<double>());
      sparsenet_features_charge_distance_10.push_back(std::vector<double>());
      sparsenet_features_charge_distance_30.push_back(std::vector<double>());


      const auto space_pts = pfpUtil.GetPFParticleSpacePoints(*daughterPFP, evt, fPFParticleTag);


      // We can calculate the number of neighbours for each space point with some radius 
      // Store the number of neighbours for each spacepoint ID, for each slice. If we 
      // only want the beam slice, or all of the slices together, this will have one 
      // element 
      std::map<unsigned int,std::vector<std::map<int,unsigned int>>> neighbourMap; 
      // Get number of neighbours for each spacepoint. Slice number 0 is for all spacepoints
      neighbourMap.insert(std::make_pair(0,GetNeighboursForRadii(evt,fNeighbourRadii, fChargeRadii, space_pts))); 

      // Get map of two neighrest spacepoints
      std::map<int,std::pair<int,int>> twoNearest = GetTwoNearestNeighbours(evt,space_pts);


      // Get charge map 
      std::map<unsigned int,float> chargeMap = GetSpacePointChargeMap(evt, "pandora");

      // Loop over all spacepoints
      for (const auto * sp : space_pts) {
        reco_daughter_spacePts_X.back().push_back(sp->XYZ()[0]);
        reco_daughter_spacePts_Y.back().push_back(sp->XYZ()[1]);
        reco_daughter_spacePts_Z.back().push_back(sp->XYZ()[2]);        
          
        
        // Add the number of neighbouring hits to the features
        for(unsigned int m = 0; m < fNeighbourRadii.size(); ++m){
          if (m==0)
            sparsenet_features_neighboring_nodes_3.back().push_back(neighbourMap.at(0)[m].at(sp->ID()));
          if (m==1)
            sparsenet_features_neighboring_nodes_10.back().push_back(neighbourMap.at(0)[m].at(sp->ID()));
          if (m==2)
            sparsenet_features_neighboring_nodes_30.back().push_back(neighbourMap.at(0)[m].at(sp->ID()));
        }
        
 
        // Add the charge map to the features 
        sparsenet_features_charge.back().push_back(chargeMap.at(sp->ID()));
        
        // Add angle and dot product between node and its two nearest neighbours
        double angle = -999.;
        double dotProduct = -999.;        
        int n1ID = twoNearest[sp->ID()].first;
        int n2ID = twoNearest[sp->ID()].second;

        // AAA: horrible hack to get the spacepoint for a given ID
        //const recob::SpacePoint n1 = *(space_pts.at(n1ID));
        //const recob::SpacePoint n2 = *(space_pts.at(n2ID));     
        recob::SpacePoint n1;
        recob::SpacePoint n2;
        for (const auto * sp : space_pts) {
          if (sp->ID() == n1ID) 
            n1 = *sp; 
          if (sp->ID() == n2ID)
            n2 = *sp;
        }        
        GetAngleAndDotProduct(*sp,n1,n2,dotProduct,angle);      
        sparsenet_features_dot_product.back().push_back(dotProduct);
        sparsenet_features_angle.back().push_back(angle);
 
        // Adding to the features vec the charge within radius
        for(unsigned int m = 0; m < fChargeRadii.size(); ++m){
          if (m==0)
            sparsenet_features_charge_distance_3.back().push_back(neighbourMap.at(0)[m+fNeighbourRadii.size()].at(sp->ID()));
          if (m==1)
            sparsenet_features_charge_distance_10.back().push_back(neighbourMap.at(0)[m+fNeighbourRadii.size()].at(sp->ID()));
          if (m==2)
            sparsenet_features_charge_distance_30.back().push_back(neighbourMap.at(0)[m+fNeighbourRadii.size()].at(sp->ID()));
 
        } 


      } // Loop over spacepoints



    } // SaveHits

    //Getting the default pandora reconstruction type (shower/track)
    const recob::Track * pandoraTrack = pfpUtil.GetPFParticleTrack(
        *daughterPFP, evt, fPFParticleTag, "pandoraTrack");
    const recob::Shower * pandoraShower = pfpUtil.GetPFParticleShower(
        *daughterPFP, evt, fPFParticleTag, "pandoraShower");
    if (pandoraTrack != 0x0) {
      reco_daughter_pandora_type.push_back(13);
    }
    else if (pandoraShower != 0x0) {
      reco_daughter_pandora_type.push_back(11);
    }
    else {
      reco_daughter_pandora_type.push_back(-1);
    }


    //If it exists (might not need this check anymore), get the forced tracking from pandora
    try{
      const recob::Track* pandora2Track = pfpUtil.GetPFParticleTrack( *daughterPFP, evt, fPFParticleTag, "pandora2Track" );
      if (fVerbose) std::cout << "pandora2 track: " << pandora2Track << std::endl;

      if( pandora2Track ){
        reco_daughter_allTrack_ID.push_back( pandora2Track->ID() );

        //Calorimetry info
        auto dummy_caloSCE =
            trackUtil.GetRecoTrackCalorimetry(
                *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE);
        bool found_calo = false;
        size_t index = 0;
        for ( index = 0; index < dummy_caloSCE.size(); ++index) {
          if (dummy_caloSCE[index].PlaneID().Plane == 2) {
            found_calo = true;
            break; 
          }
        }

        reco_daughter_allTrack_resRange_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_dEdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_dQdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_calibrated_dQdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_calibrated_dEdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_EField_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_calo_X.push_back( std::vector<double>() );
        reco_daughter_allTrack_calo_Y.push_back( std::vector<double>() );
        reco_daughter_allTrack_calo_Z.push_back( std::vector<double>() );

        if (found_calo) {
          auto dummy_dEdx_SCE = dummy_caloSCE[index].dEdx();
          auto dummy_dQdx_SCE = dummy_caloSCE[index].dQdx();
          auto dummy_Range_SCE = dummy_caloSCE[index].ResidualRange();
          auto theXYZPoints = dummy_caloSCE[index].XYZ();

          //Momentum by range calculation using SCE corrected info
          reco_daughter_allTrack_momByRange_alt_proton.push_back(
              track_p_calc.GetTrackMomentum(dummy_caloSCE[index].Range(),
              2212));
          reco_daughter_allTrack_momByRange_alt_muon.push_back(
              track_p_calc.GetTrackMomentum(dummy_caloSCE[index].Range(),
              13));
          reco_daughter_allTrack_alt_len.push_back(dummy_caloSCE[index].Range() );

          for( size_t j = 0; j < dummy_dEdx_SCE.size(); ++j ){
            reco_daughter_allTrack_resRange_SCE.back().push_back( dummy_Range_SCE[j] );
            reco_daughter_allTrack_dEdX_SCE.back().push_back( dummy_dEdx_SCE[j] );
            reco_daughter_allTrack_dQdX_SCE.back().push_back( dummy_dQdx_SCE[j] );
            reco_daughter_allTrack_calo_X.back().push_back(theXYZPoints[j].X());
            reco_daughter_allTrack_calo_Y.back().push_back(theXYZPoints[j].Y());
            reco_daughter_allTrack_calo_Z.back().push_back(theXYZPoints[j].Z());
          }

          if (fRecalibrate) {
            std::vector<float> cali_dEdX_SCE
                = calibration_SCE.GetCalibratedCalorimetry(
                    *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 2);

            for( size_t j = 0; j < cali_dEdX_SCE.size(); ++j ){
              reco_daughter_allTrack_calibrated_dEdX_SCE.back().push_back(
                  cali_dEdX_SCE[j]);
            }
            std::vector<double> new_dQdX = calibration_SCE.CalibratedQdX(
                *pandora2Track, evt, "pandora2Track",
                fPandora2CaloSCE, 2, -10.);
            for (auto dqdx : new_dQdX) {
              reco_daughter_allTrack_calibrated_dQdX_SCE.back().push_back(dqdx);
            }
            std::vector<double> efield = calibration_SCE.GetEFieldVector(
                *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 2, -10.);
            for (auto ef : efield) {
              reco_daughter_allTrack_EField_SCE.back().push_back(ef);
            }

          }
          else {
            for (size_t j = 0; j < dummy_dQdx_SCE.size(); ++j) {
              reco_daughter_allTrack_calibrated_dEdX_SCE.back().push_back(
                  dummy_dEdx_SCE[j]);
              reco_daughter_allTrack_calibrated_dQdX_SCE.back().push_back(
                  dummy_dQdx_SCE[j]);
              double E_field_nominal = detProp.Efield();   // Electric Field in the drift region in KV/cm
              geo::Vector_t E_field_offsets
                  = (sce->EnableCalEfieldSCE() && fSCE ?
                     sce->GetCalEfieldOffsets(
                         geo::Point_t{theXYZPoints[j].X(), theXYZPoints[j].Y(),
                                      theXYZPoints[j].Z()},
                         dummy_caloSCE[index].PlaneID().TPC) :
                     geo::Vector_t{0., 0., 0.});
              TVector3 E_field_vector
                  = {E_field_nominal*(1 + E_field_offsets.X()),
                     E_field_nominal*E_field_offsets.Y(),
                     E_field_nominal*E_field_offsets.Z()};

              double E_field = E_field_vector.Mag();
              reco_daughter_allTrack_EField_SCE.back().push_back(E_field);
            }
          }

          //Chi2 based PID
          std::pair<double, int> this_chi2_ndof = trackUtil.Chi2PID(
              reco_daughter_allTrack_calibrated_dEdX_SCE.back(),
              reco_daughter_allTrack_resRange_SCE.back(), templates[2212]);

          reco_daughter_allTrack_Chi2_proton.push_back(this_chi2_ndof.first);
          reco_daughter_allTrack_Chi2_ndof.push_back(this_chi2_ndof.second);

          this_chi2_ndof = trackUtil.Chi2PID(
              reco_daughter_allTrack_calibrated_dEdX_SCE.back(),
              reco_daughter_allTrack_resRange_SCE.back(), templates[211]);

          reco_daughter_allTrack_Chi2_pion.push_back(this_chi2_ndof.first);
          reco_daughter_allTrack_Chi2_ndof_pion.push_back(this_chi2_ndof.second);

          this_chi2_ndof = trackUtil.Chi2PID(
              reco_daughter_allTrack_calibrated_dEdX_SCE.back(),
              reco_daughter_allTrack_resRange_SCE.back(), templates[13]);

          reco_daughter_allTrack_Chi2_muon.push_back(this_chi2_ndof.first);
          reco_daughter_allTrack_Chi2_ndof_muon.push_back(this_chi2_ndof.second);
        }
        else {
          reco_daughter_allTrack_momByRange_alt_proton.push_back(-999.);
          reco_daughter_allTrack_momByRange_alt_muon.push_back(-999.);
          reco_daughter_allTrack_alt_len.push_back(-999.);
          reco_daughter_allTrack_Chi2_proton.push_back(-999.);
          reco_daughter_allTrack_Chi2_ndof.push_back(-999);

          reco_daughter_allTrack_Chi2_pion.push_back(-999.);
          reco_daughter_allTrack_Chi2_ndof_pion.push_back(-999);

          reco_daughter_allTrack_Chi2_muon.push_back(-999.);
          reco_daughter_allTrack_Chi2_ndof_muon.push_back(-999);
        }

        //Calorimetry + chi2 for planes 0 and 1
        size_t plane0_index = 0;
        bool found_plane0 = false;
        for ( plane0_index = 0; plane0_index < dummy_caloSCE.size(); ++plane0_index) {
          if (dummy_caloSCE[plane0_index].PlaneID().Plane == 0) {
            found_plane0 = true;
            break; 
          }
        }
        size_t plane1_index = 0;
        bool found_plane1 = false;
        for ( plane1_index = 0; plane1_index < dummy_caloSCE.size(); ++plane1_index) {
          if (dummy_caloSCE[plane1_index].PlaneID().Plane == 1) {
            found_plane1 = true;
            break; 
          }
        }
        if (fVerbose)
            std::cout << "Plane 0, 1 " << plane0_index << " " <<
                         plane1_index << std::endl;


        reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.push_back(
            std::vector<double>());
        reco_daughter_allTrack_resRange_plane0.push_back(
            std::vector<double>());

        if (found_plane0) {
          auto resRange_plane0 = dummy_caloSCE[plane0_index].ResidualRange();
          //auto dEdX_plane0 = dummy_caloSCE[plane0_index].dEdx();

          for (size_t j = 0; j < resRange_plane0.size(); ++j) {
            reco_daughter_allTrack_resRange_plane0.back().push_back(
                resRange_plane0[j]);
          }

 
          std::vector<float> dEdX_plane0
              = (fRecalibrate ?
                 calibration_SCE.GetCalibratedCalorimetry(
                    *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 0) :
                 dummy_caloSCE[plane0_index].dEdx());
          for (auto & dedx : dEdX_plane0) {
            reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.back().push_back(
                dedx);
          }
          /*
          if (fRecalibrate){
            std::vector<float> dEdX_plane0 = calibration_SCE.GetCalibratedCalorimetry(
              *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 0);
            for (size_t j = 0; j < dEdX_plane0.size(); ++j) {
              reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.back().push_back(dEdX_plane0[j]);
            }
          }
          else{
            for (size_t j = 0; j < dEdX_plane0.size(); ++j){
              reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.back().push_back(dEdX_plane0[j]);
            }
          }*/

          std::pair<double, int> plane0_chi2_ndof = trackUtil.Chi2PID(
              reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.back(),
              reco_daughter_allTrack_resRange_plane0.back(), templates[2212]);
          reco_daughter_allTrack_Chi2_proton_plane0.push_back(
              plane0_chi2_ndof.first);
          reco_daughter_allTrack_Chi2_ndof_plane0.push_back(
              plane0_chi2_ndof.second);
        }
        else {
          reco_daughter_allTrack_Chi2_proton_plane0.push_back(
              -999.);
          reco_daughter_allTrack_Chi2_ndof_plane0.push_back(
              -999);
        }


        reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.push_back(
            std::vector<double>());

        reco_daughter_allTrack_resRange_plane1.push_back(
            std::vector<double>());

        if (found_plane1) {
          auto resRange_plane1 = dummy_caloSCE[plane1_index].ResidualRange();
          for (size_t j = 0; j < resRange_plane1.size(); ++j) {
            reco_daughter_allTrack_resRange_plane1.back().push_back(
                resRange_plane1[j]);
          }


          std::vector<float> dEdX_plane1
              = (fRecalibrate ? 
                 calibration_SCE.GetCalibratedCalorimetry(
                    *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 1) :
                 dummy_caloSCE[plane1_index].dEdx());
          for (auto & dedx : dEdX_plane1) {
            reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.back().push_back(
                dedx);
          }

          /*
          if (fRecalibrate){
            std::vector<float> dEdX_plane1 = calibration_SCE.GetCalibratedCalorimetry(
              *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 1);

            for (size_t j = 0; j < dEdX_plane1.size(); ++j) {
              reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.back().push_back(dEdX_plane1[j]);
            }
          }
          else{
            auto dEdX_plane1 = dummy_caloSCE[plane1_index].dEdx();
            for (size_t j = 0; j < dEdX_plane1.size(); ++j){
              reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.back().push_back(dEdX_plane1[j]);
            }
          } */

          std::pair<double, int> plane1_chi2_ndof = trackUtil.Chi2PID(
              reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.back(),
              reco_daughter_allTrack_resRange_plane1.back(), templates[2212]);
          reco_daughter_allTrack_Chi2_proton_plane1.push_back(
              plane1_chi2_ndof.first);
          reco_daughter_allTrack_Chi2_ndof_plane1.push_back(
              plane1_chi2_ndof.second);
        }
        else {
          reco_daughter_allTrack_Chi2_proton_plane1.push_back(
              -999.);
          reco_daughter_allTrack_Chi2_ndof_plane1.push_back(
              -999);
        }
        //////////////////////////////////////

        //Spatial + direction stuff from Pandora
        reco_daughter_allTrack_Theta.push_back(  pandora2Track->Theta() );
        reco_daughter_allTrack_Phi.push_back(  pandora2Track->Phi() );

        reco_daughter_allTrack_startDirX.push_back(pandora2Track->StartDirection().X());
        reco_daughter_allTrack_startDirY.push_back(pandora2Track->StartDirection().Y());
        reco_daughter_allTrack_startDirZ.push_back(pandora2Track->StartDirection().Z());

        reco_daughter_allTrack_len.push_back(    pandora2Track->Length() );
        reco_daughter_allTrack_startX.push_back( pandora2Track->Trajectory().Start().X() );
        reco_daughter_allTrack_startY.push_back( pandora2Track->Trajectory().Start().Y() );
        reco_daughter_allTrack_startZ.push_back( pandora2Track->Trajectory().Start().Z() );
        reco_daughter_allTrack_endX.push_back(   pandora2Track->Trajectory().End().X() );
        reco_daughter_allTrack_endY.push_back(   pandora2Track->Trajectory().End().Y() );
        reco_daughter_allTrack_endZ.push_back(   pandora2Track->Trajectory().End().Z() );

        //Using new michel tagging from Ajib
        if (!fSkipMVA) {
          std::pair<double, int> vertex_results =
              trackUtil.GetVertexMichelScore(
                  *pandora2Track, evt, "pandora2Track", fHitTag,
                  0., -500., 500., 0., 500., 0., false,
                  reco_beam_endX, reco_beam_endY, reco_beam_endZ);

          reco_daughter_allTrack_vertex_michel_score.push_back(
              vertex_results.first);
          reco_daughter_allTrack_vertex_nHits.push_back(
              vertex_results.second);
        }
        else {
          reco_daughter_allTrack_vertex_michel_score.push_back(
              -999.);
          reco_daughter_allTrack_vertex_nHits.push_back(
              -999);
        }

        if (fVerbose) std::cout << "pandora2Length " << pandora2Track->Length() << std::endl;
        //Another momentum by range calculation, but using the SCE-uncorrected info
        reco_daughter_allTrack_momByRange_proton.push_back( track_p_calc.GetTrackMomentum( pandora2Track->Length(), 2212 ) );
        reco_daughter_allTrack_momByRange_muon.push_back(   track_p_calc.GetTrackMomentum( pandora2Track->Length(), 13  ) );

      }
      else{  //Defaults
        reco_daughter_allTrack_ID.push_back( -999 );
        reco_daughter_allTrack_resRange_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_dEdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_dQdX_SCE.push_back( std::vector<double>() );
        reco_daughter_allTrack_calibrated_dEdX_SCE.push_back(std::vector<double>());
        reco_daughter_allTrack_calibrated_dQdX_SCE.push_back(std::vector<double>());
        reco_daughter_allTrack_EField_SCE.push_back(std::vector<double>());
        reco_daughter_allTrack_Chi2_proton.push_back( -999. );
        reco_daughter_allTrack_Chi2_ndof.push_back( -999 );

        reco_daughter_allTrack_Chi2_pion.push_back( -999. );
        reco_daughter_allTrack_Chi2_ndof_pion.push_back( -999 );

        reco_daughter_allTrack_Chi2_muon.push_back( -999. );
        reco_daughter_allTrack_Chi2_ndof_muon.push_back( -999 );

        reco_daughter_allTrack_calo_X.push_back( std::vector<double>() );
        reco_daughter_allTrack_calo_Y.push_back( std::vector<double>() );
        reco_daughter_allTrack_calo_Z.push_back( std::vector<double>() );

        //Calorimetry + chi2 for planes 0 and 1
        reco_daughter_allTrack_calibrated_dEdX_SCE_plane0.push_back(
            std::vector<double>());
        reco_daughter_allTrack_calibrated_dEdX_SCE_plane1.push_back(
            std::vector<double>());
        reco_daughter_allTrack_resRange_plane0.push_back(
            std::vector<double>());
        reco_daughter_allTrack_resRange_plane1.push_back(
            std::vector<double>());

        reco_daughter_allTrack_Chi2_proton_plane0.push_back( -999. );
        reco_daughter_allTrack_Chi2_ndof_plane0.push_back( -999 );
        reco_daughter_allTrack_Chi2_proton_plane1.push_back( -999. );
        reco_daughter_allTrack_Chi2_ndof_plane1.push_back( -999 );
        //////////////////////////////////

        reco_daughter_allTrack_Theta.push_back(-999. );
        reco_daughter_allTrack_Phi.push_back(-999.);
        reco_daughter_allTrack_len.push_back( -999. );
        reco_daughter_allTrack_alt_len.push_back(-999.);
        reco_daughter_allTrack_startX.push_back( -999. );
        reco_daughter_allTrack_startY.push_back( -999. );
        reco_daughter_allTrack_startZ.push_back( -999. );
        reco_daughter_allTrack_endX.push_back(   -999. );
        reco_daughter_allTrack_endY.push_back(   -999. );
        reco_daughter_allTrack_endZ.push_back(   -999. );

        reco_daughter_allTrack_startDirX.push_back(-999.); 
        reco_daughter_allTrack_startDirY.push_back(-999.); 
        reco_daughter_allTrack_startDirZ.push_back(-999.); 
        reco_daughter_allTrack_momByRange_proton.push_back(-999.);
        reco_daughter_allTrack_momByRange_muon.push_back(-999.);

        reco_daughter_allTrack_momByRange_alt_proton.push_back(-999.);
        reco_daughter_allTrack_momByRange_alt_muon.push_back(-999.);

        reco_daughter_allTrack_vertex_michel_score.push_back(-999.);
        reco_daughter_allTrack_vertex_nHits.push_back(-999);

      }
    }
    catch( const cet::exception &e ){
      MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Track object not found, moving on" << std::endl;
    }


    //If it exists (might not need this check anymore), get the forced shower object from pandora
    try{
      const recob::Shower* pandora2Shower = pfpUtil.GetPFParticleShower( *daughterPFP, evt, fPFParticleTag, "pandora2Shower" );
      if (fVerbose) std::cout << "pandora2 shower: " << pandora2Shower << std::endl;

      if( pandora2Shower ){
        reco_daughter_allShower_ID.push_back(     pandora2Shower->ID() );
        reco_daughter_allShower_len.push_back(    pandora2Shower->Length() );
        reco_daughter_allShower_startX.push_back( pandora2Shower->ShowerStart().X() );
        reco_daughter_allShower_startY.push_back( pandora2Shower->ShowerStart().Y() );
        reco_daughter_allShower_startZ.push_back( pandora2Shower->ShowerStart().Z() );

        reco_daughter_allShower_dirX.push_back( pandora2Shower->Direction().X() );
        reco_daughter_allShower_dirY.push_back( pandora2Shower->Direction().Y() );
        reco_daughter_allShower_dirZ.push_back( pandora2Shower->Direction().Z() );


        const std::vector<art::Ptr<recob::Hit>> hits =
            showerUtil.GetRecoShowerArtHits(
                *pandora2Shower, evt, "pandora2Shower");

        art::FindManyP<recob::SpacePoint> spFromHits(hits, evt, fHitTag);
        //double total_shower_energy = 0.;
        //need to get average y
        std::vector<double> x_vec, y_vec, z_vec;
        double total_y = 0.;
        int n_good_y = 0;
        std::vector<art::Ptr<recob::Hit>> good_hits;

        for (size_t iHit = 0; iHit < hits.size(); ++iHit) {
          auto theHit = hits[iHit];
          if (theHit->View() != 2) continue; //skip induction planes

          good_hits.push_back(theHit);

          double shower_hit_x = detProp.ConvertTicksToX(
              theHit->PeakTime(),
              theHit->WireID().Plane,
              theHit->WireID().TPC, 0);

          double shower_hit_z = wireReadout.Wire(theHit->WireID()).GetCenter().Z();

          x_vec.push_back(shower_hit_x);
          z_vec.push_back(shower_hit_z);

          std::vector<art::Ptr<recob::SpacePoint>> sps = spFromHits.at(iHit);
          //std::cout << shower_hit_x << " " << shower_hit_z << " ";
          if (!sps.empty()) {
            y_vec.push_back(sps[0]->XYZ()[1]);
            total_y += y_vec.back();
            ++n_good_y;
            //std::cout << shower_hit_y_vec.back();
          }
          else {
           y_vec.push_back(-999.);
          }
          //std::cout << std::endl;
        }

        if (fGetCalibratedShowerEnergy) {
          //std::cout << "Getting calibrated shower energy" << std::endl;
          auto calo = showerUtil.GetRecoShowerCalorimetry(
              *pandora2Shower, evt, "pandora2Shower", "pandora2cali");
          bool found_calo = false;
          size_t index = 0;
          for (index = 0; index < calo.size(); ++index) {
            if (calo[index].PlaneID().Plane == 2) {
              found_calo = true;
              break; 
            }
          }
          
          if (!found_calo) {
            reco_daughter_allShower_calibrated_energy.push_back(-999.);
          }
          else {
            reco_daughter_allShower_calibrated_energy.push_back(calo[index].KineticEnergy());
          }
        }

        if (fGetUncalibratedShowerEnergy && n_good_y > 1) {
          double total_shower_energy = 0.;
          for (size_t iHit = 0; iHit < good_hits.size(); ++iHit) {
            auto const& theHit = good_hits[iHit];
            if (theHit->View() != 2) continue; //skip induction planes

            if (y_vec[iHit] < -100.)
              y_vec[iHit] = total_y / n_good_y;

            total_shower_energy += calibration_SCE.HitToEnergy(
                *good_hits[iHit], x_vec[iHit], y_vec[iHit], z_vec[iHit]);
          }
          reco_daughter_allShower_energy.push_back(total_shower_energy);
        }
        else {
          reco_daughter_allShower_energy.push_back(-999.);
        }
        
      }
      else{
        reco_daughter_allShower_ID.push_back(-999);
        reco_daughter_allShower_len.push_back(-999.);
        reco_daughter_allShower_startX.push_back(-999.);
        reco_daughter_allShower_startY.push_back(-999.);
        reco_daughter_allShower_startZ.push_back(-999.);
        reco_daughter_allShower_dirX.push_back(-999.);
        reco_daughter_allShower_dirY.push_back(-999.);
        reco_daughter_allShower_dirZ.push_back(-999.);
        reco_daughter_allShower_energy.push_back(-999.);
        reco_daughter_allShower_calibrated_energy.push_back(-999.);
      }

    }
    catch( const cet::exception &e ){
      MF_LOG_WARNING("PDSPAnalyzer") << "pandora2Shower object not found, moving on" << std::endl;
    }
  }
}

void pduneana::PDSPAnalyzer::BeamMatchInfo(
    const art::Event & evt, const recob::PFParticle * particle,
    const simb::MCParticle * true_beam_particle,
    detinfo::DetectorClocksData const& clockData) {
  protoana::MCParticleSharedHits beam_match  = truthUtil.GetMCParticleByHits( clockData, *particle, evt, fPFParticleTag, fHitTag );
  if( beam_match.particle ){
    reco_beam_true_byHits_matched = ( beam_match.particle->TrackId() == true_beam_particle->TrackId() );
    reco_beam_true_byHits_PDG = beam_match.particle->PdgCode();
    reco_beam_true_byHits_ID = beam_match.particle->TrackId();

    reco_beam_true_byHits_process = beam_match.particle->Process();
    reco_beam_true_byHits_endProcess = beam_match.particle->EndProcess();
    reco_beam_true_byHits_origin = pi_serv->TrackIdToMCTruth_P(beam_match.particle->TrackId())->Origin();

    reco_beam_true_byHits_startPx = beam_match.particle->Px();
    reco_beam_true_byHits_startPy = beam_match.particle->Py();
    reco_beam_true_byHits_startPz = beam_match.particle->Pz();
    reco_beam_true_byHits_startP  = sqrt( reco_beam_true_byHits_startPx*reco_beam_true_byHits_startPx
                                   + reco_beam_true_byHits_startPy*reco_beam_true_byHits_startPy
                                   + reco_beam_true_byHits_startPz*reco_beam_true_byHits_startPz );
    reco_beam_true_byHits_startE = beam_match.particle->E();

    size_t np = beam_match.particle->NumberTrajectoryPoints();
    if( np > 1 ){
      reco_beam_true_byHits_endPx = beam_match.particle->Px( np - 2 );
      reco_beam_true_byHits_endPy = beam_match.particle->Py( np - 2 );
      reco_beam_true_byHits_endPz = beam_match.particle->Pz( np - 2 );
      reco_beam_true_byHits_endP  = sqrt( reco_beam_true_byHits_endPx*reco_beam_true_byHits_endPx
                                   + reco_beam_true_byHits_endPy*reco_beam_true_byHits_endPy
                                   + reco_beam_true_byHits_endPz*reco_beam_true_byHits_endPz );
      reco_beam_true_byHits_endE  = beam_match.particle->E( np - 2 );
    }

    auto list = truthUtil.GetMCParticleListByHits( clockData, *particle, evt, fPFParticleTag, fHitTag );
    double total = 0.;
    double matched_hits = 0.;
    for( size_t j = 0; j < list.size(); ++j ){
    //  std::cout << "Contrib " << j << " " << list[j].first->TrackId() << " " << list[j].second << std::endl;
      //std::cout << "Contrib " << j << " " << list[j].particle->TrackId() << " " << list[j].particle->PdgCode()
      //           << " " << pi_serv->TrackIdToMCTruth_P(list[j].particle->TrackId())->Origin()
      //           << " " << list[j].nSharedHits << " " << list[j].nSharedDeltaRayHits << std::endl;

      if( list[j].particle == beam_match.particle ){
         matched_hits = list[j].nSharedHits + list[j].nSharedDeltaRayHits;
      }

      total += list[j].nSharedHits + list[j].nSharedDeltaRayHits;
    }

    reco_beam_true_byHits_purity = ( matched_hits / total );

  }

  const simb::MCParticle* trueParticle = 0x0;
  trueParticle = truthUtil.GetMCParticleFromPFParticle(clockData, *particle, evt, fPFParticleTag);
  if( trueParticle ){

    //Check that this is the correct true particle
    if( trueParticle->TrackId() == true_beam_particle->TrackId() ){
      reco_beam_true_byE_matched = true;
    }

    reco_beam_true_byE_PDG = trueParticle->PdgCode();
    reco_beam_true_byE_ID = trueParticle->TrackId();

    reco_beam_true_byE_process = trueParticle->Process();
    reco_beam_true_byE_endProcess = trueParticle->EndProcess();
    reco_beam_true_byE_origin = pi_serv->TrackIdToMCTruth_P(trueParticle->TrackId())->Origin();

    reco_beam_true_byE_startPx = trueParticle->Px();
    reco_beam_true_byE_startPy = trueParticle->Py();
    reco_beam_true_byE_startPz = trueParticle->Pz();
    reco_beam_true_byE_startP  = sqrt( reco_beam_true_byE_startPx*reco_beam_true_byE_startPx
                                   + reco_beam_true_byE_startPy*reco_beam_true_byE_startPy
                                   + reco_beam_true_byE_startPz*reco_beam_true_byE_startPz );
    reco_beam_true_byE_startE = trueParticle->E();

    size_t np = trueParticle->NumberTrajectoryPoints();
    if( np > 1 ){
      reco_beam_true_byE_endPx = trueParticle->Px( np - 2 );
      reco_beam_true_byE_endPy = trueParticle->Py( np - 2 );
      reco_beam_true_byE_endPz = trueParticle->Pz( np - 2 );
      reco_beam_true_byE_endP  = sqrt( reco_beam_true_byE_endPx*reco_beam_true_byE_endPx
                                   + reco_beam_true_byE_endPy*reco_beam_true_byE_endPy
                                   + reco_beam_true_byE_endPz*reco_beam_true_byE_endPz );
      reco_beam_true_byE_endE  = trueParticle->E( np - 2 );
    }

  }
}


//Gets true info matched to the reconstructed daughter PFPs
void pduneana::PDSPAnalyzer::DaughterMatchInfo(
    const art::Event & evt, const recob::PFParticle * daughterPFP,
    detinfo::DetectorClocksData const& clockData) {

  const sim::ParticleList & plist = pi_serv->ParticleList();
  protoana::MCParticleSharedHits match = truthUtil.GetMCParticleByHits( clockData, *daughterPFP, evt, fPFParticleTag, fHitTag );
  if( match.particle ){

    reco_daughter_PFP_true_byHits_PDG.push_back( match.particle->PdgCode() );
    reco_daughter_PFP_true_byHits_ID.push_back( match.particle->TrackId() );
    reco_daughter_PFP_true_byHits_parID.push_back( match.particle->Mother() );
    reco_daughter_PFP_true_byHits_parPDG.push_back(
      ( (match.particle->Mother() > 0) ? plist[ match.particle->Mother() ]->PdgCode() : 0 )
    );

    reco_daughter_PFP_true_byHits_process.push_back( match.particle->Process() );
    reco_daughter_PFP_true_byHits_origin.push_back(
      pi_serv->TrackIdToMCTruth_P(match.particle->TrackId())->Origin()
    );
    reco_daughter_PFP_true_byHits_sharedHits.push_back( match.nSharedHits );
    reco_daughter_PFP_true_byHits_emHits.push_back( match.nSharedDeltaRayHits );

    reco_daughter_PFP_true_byHits_len.push_back( match.particle->Trajectory().TotalLength() );
    reco_daughter_PFP_true_byHits_startX.push_back( match.particle->Position(0).X() );
    reco_daughter_PFP_true_byHits_startY.push_back( match.particle->Position(0).Y() );
    reco_daughter_PFP_true_byHits_startZ.push_back( match.particle->Position(0).Z() );

    reco_daughter_PFP_true_byHits_endX.push_back( match.particle->EndPosition().X() );
    reco_daughter_PFP_true_byHits_endY.push_back( match.particle->EndPosition().Y() );
    reco_daughter_PFP_true_byHits_endZ.push_back( match.particle->EndPosition().Z() );

    reco_daughter_PFP_true_byHits_startPx.push_back( match.particle->Px() );
    reco_daughter_PFP_true_byHits_startPy.push_back( match.particle->Py() );
    reco_daughter_PFP_true_byHits_startPz.push_back( match.particle->Pz() );
    reco_daughter_PFP_true_byHits_startE.push_back( match.particle->E() );
    reco_daughter_PFP_true_byHits_startP.push_back(
                    sqrt(match.particle->Px()*match.particle->Px() +
                            match.particle->Py()*match.particle->Py() +
                            match.particle->Pz()*match.particle->Pz()) );
    reco_daughter_PFP_true_byHits_endProcess.push_back( match.particle->EndProcess());

    auto list = truthUtil.GetMCParticleListByHits( clockData, *daughterPFP, evt, fPFParticleTag, fHitTag );
    double total = 0.;
    double matched_hits = 0.;
    for( size_t j = 0; j < list.size(); ++j ){
      if( list[j].particle == match.particle ){
         matched_hits = list[j].nSharedHits + list[j].nSharedDeltaRayHits;
      }
      total += list[j].nSharedHits + list[j].nSharedDeltaRayHits;
    }

    reco_daughter_PFP_true_byHits_purity.push_back( matched_hits / total );

    double totalTruth = truthUtil.GetMCParticleHits( clockData, *match.particle, evt, fHitTag).size();
    double sharedHits = truthUtil.GetSharedHits( clockData, *match.particle, *daughterPFP, evt, fPFParticleTag).size();
    reco_daughter_PFP_true_byHits_completeness.push_back( sharedHits/totalTruth );
  }
  else{
    reco_daughter_PFP_true_byHits_PDG.push_back( -999 );
    reco_daughter_PFP_true_byHits_ID.push_back( -999 );
    reco_daughter_PFP_true_byHits_origin.push_back( -999 );
    reco_daughter_PFP_true_byHits_parID.push_back( -999 );
    reco_daughter_PFP_true_byHits_parPDG.push_back( -999 );
    reco_daughter_PFP_true_byHits_process.push_back( "empty" );
    reco_daughter_PFP_true_byHits_sharedHits.push_back( -999 );
    reco_daughter_PFP_true_byHits_emHits.push_back( -999 );

    reco_daughter_PFP_true_byHits_len.push_back( -999. );
    reco_daughter_PFP_true_byHits_startX.push_back( -999. );
    reco_daughter_PFP_true_byHits_startY.push_back( -999. );
    reco_daughter_PFP_true_byHits_startZ.push_back( -999. );
    reco_daughter_PFP_true_byHits_endX.push_back( -999. );
    reco_daughter_PFP_true_byHits_endY.push_back( -999. );
    reco_daughter_PFP_true_byHits_endZ.push_back( -999. );
    reco_daughter_PFP_true_byHits_startPx.push_back( -999. );
    reco_daughter_PFP_true_byHits_startPy.push_back( -999. );
    reco_daughter_PFP_true_byHits_startPz.push_back( -999. );
    reco_daughter_PFP_true_byHits_startP.push_back( -999. );
    reco_daughter_PFP_true_byHits_startE.push_back( -999. );
    reco_daughter_PFP_true_byHits_endProcess.push_back("empty");
    reco_daughter_PFP_true_byHits_purity.push_back( -999. );
    reco_daughter_PFP_true_byHits_completeness.push_back( -999. );
  }

  //Matching by energy
  const simb::MCParticle* true_daughter_byE = truthUtil.GetMCParticleFromPFParticle(clockData, *daughterPFP, evt, fPFParticleTag);
  if( true_daughter_byE ){
    reco_daughter_PFP_true_byE_PDG.push_back( true_daughter_byE->PdgCode() );
    reco_daughter_PFP_true_byE_len.push_back( true_daughter_byE->Trajectory().TotalLength() );
    double purity = truthUtil.GetPurity( clockData, *daughterPFP, evt, fPFParticleTag);
    double completeness = truthUtil.GetCompleteness( clockData, *daughterPFP, evt, fPFParticleTag, fHitTag );
    reco_daughter_PFP_true_byE_purity.push_back( purity );
    reco_daughter_PFP_true_byE_completeness.push_back( completeness );
  }
  else {
    reco_daughter_PFP_true_byE_PDG.push_back( -999 );
    reco_daughter_PFP_true_byE_len.push_back( -999. );
    reco_daughter_PFP_true_byE_purity.push_back( -999. );
    reco_daughter_PFP_true_byE_completeness.push_back( -999. );
  }

}

//Info from the forced pandora tracking applied to the beam PFParticle
void pduneana::PDSPAnalyzer::BeamForcedTrackInfo(
    const art::Event & evt, const recob::PFParticle * particle) {
  try{
    const recob::Track* pandora2Track = pfpUtil.GetPFParticleTrack( *particle, evt, fPFParticleTag, "pandora2Track" );
    if (fVerbose) std::cout << "pandora2 track: " << pandora2Track << std::endl;

    if( pandora2Track ){
      // michle score
      if (!fSkipMVA) {
        std::pair<double, int> vertex_michel_score = trackUtil.GetVertexMichelScore(*pandora2Track, evt, "pandora2Track", fHitTag);
        reco_beam_vertex_nHits_allTrack = vertex_michel_score.second;
        reco_beam_vertex_michel_score_allTrack = vertex_michel_score.first;
        
        std::pair<double, double> vertex_michel_score_weight_by_charge =
            trackUtil.GetVertexMichelScore_weight_by_charge(*pandora2Track, evt, "pandora2Track", fHitTag);
        if (vertex_michel_score_weight_by_charge.second != 0)
          reco_beam_vertex_michel_score_weight_by_charge_allTrack = vertex_michel_score_weight_by_charge.first/vertex_michel_score_weight_by_charge.second;
        else
          reco_beam_vertex_michel_score_weight_by_charge_allTrack = -999.;
      }
      reco_beam_allTrack_ID = pandora2Track->ID();
      reco_beam_allTrack_beam_cuts = beam_cuts.IsBeamlike( *pandora2Track, evt, "1" );
      reco_beam_allTrack_startX = pandora2Track->Trajectory().Start().X();
      reco_beam_allTrack_startY = pandora2Track->Trajectory().Start().Y();
      reco_beam_allTrack_startZ = pandora2Track->Trajectory().Start().Z();
      reco_beam_allTrack_endX = pandora2Track->Trajectory().End().X();
      reco_beam_allTrack_endY = pandora2Track->Trajectory().End().Y();
      reco_beam_allTrack_endZ = pandora2Track->Trajectory().End().Z();

      auto startDir = pandora2Track->StartDirection();
      auto endDir   = pandora2Track->EndDirection();

      //try flipping
      if( reco_beam_allTrack_startZ > reco_beam_allTrack_endZ ){
        reco_beam_allTrack_flipped = true;
        reco_beam_allTrack_endX = pandora2Track->Trajectory().Start().X();
        reco_beam_allTrack_endY = pandora2Track->Trajectory().Start().Y();
        reco_beam_allTrack_endZ = pandora2Track->Trajectory().Start().Z();
        reco_beam_allTrack_startX = pandora2Track->Trajectory().End().X();
        reco_beam_allTrack_startY = pandora2Track->Trajectory().End().Y();
        reco_beam_allTrack_startZ = pandora2Track->Trajectory().End().Z();

        reco_beam_allTrack_trackDirX =  -1. * endDir.X();
        reco_beam_allTrack_trackDirY =  -1. * endDir.Y();
        reco_beam_allTrack_trackDirZ =  -1. * endDir.Z();

        reco_beam_allTrack_trackEndDirX =  -1. * startDir.X();
        reco_beam_allTrack_trackEndDirY =  -1. * startDir.Y();
        reco_beam_allTrack_trackEndDirZ =  -1. * startDir.Z();
      }
      else{
        reco_beam_allTrack_flipped = false;
        reco_beam_allTrack_trackDirX    =  startDir.X();
        reco_beam_allTrack_trackDirY    =  startDir.Y();
        reco_beam_allTrack_trackDirZ    =  startDir.Z();
        reco_beam_allTrack_trackEndDirX =  endDir.X();
        reco_beam_allTrack_trackEndDirY =  endDir.Y();
        reco_beam_allTrack_trackEndDirZ =  endDir.Z();
      }

      reco_beam_allTrack_len  = pandora2Track->Length();

      auto allHits = evt.getValidHandle<std::vector<recob::Hit> >(fHitTag);
      auto calo = trackUtil.GetRecoTrackCalorimetry(*pandora2Track, evt, "pandora2Track", fPandora2CaloSCE);
      size_t index = 0;
      bool found_plane = false;
      for (index = 0; index < calo.size(); ++index) {
        if (calo[index].PlaneID().Plane == 2) {
          found_plane = true;
          break; 
        }
      }

      if (found_plane) {
        reco_beam_alt_len_allTrack = calo[index].Range();
        auto calo_range = calo[index].ResidualRange();
        auto calo_dEdX = calo[index].dEdx();
        auto calo_dQdX = calo[index].dQdx();
        auto TpIndices = calo[index].TpIndices();

        for( size_t i = 0; i < calo_range.size(); ++i ){
          reco_beam_allTrack_resRange.push_back( calo_range[i] );
          //reco_beam_dQdX_SCE_allTrack.push_back( calo_dQdX[i] );
          //reco_beam_dEdX_SCE_allTrack.push_back( calo_dEdX[i] );
          reco_beam_TrkPitch_SCE_allTrack.push_back( calo[index].TrkPitchVec()[i] );
        }
        
        auto theXYZPoints = calo[index].XYZ();
        for( size_t i = 0; i < calo_dQdX.size(); ++i ){
          const recob::Hit & theHit = (*allHits)[ TpIndices[i] ];
          //reco_beam_dQ_allTrack.push_back(theHit.Integral());
          //reco_beam_calo_TPC_allTrack.push_back(theHit.WireID().TPC);
          if (theHit.WireID().TPC == 1) {
            reco_beam_calo_wire_allTrack.push_back( theHit.WireID().Wire );
          }
          else if (theHit.WireID().TPC == 5) {
            reco_beam_calo_wire_allTrack.push_back( theHit.WireID().Wire + 479);
          }
          //Need other TPCs?
          else {
            reco_beam_calo_wire_allTrack.push_back(theHit.WireID().Wire );
          }
          //reco_beam_calo_tick_allTrack.push_back( theHit.PeakTime() );
          //calo_hit_indices_allTrack.push_back( TpIndices[i] );
          //reco_beam_calo_wire_z_allTrack.push_back(wireReadout.Wire(theHit.WireID()).GetCenter().Z());
          
          reco_beam_calo_X_allTrack.push_back(theXYZPoints[i].X());
          reco_beam_calo_Y_allTrack.push_back(theXYZPoints[i].Y());
          reco_beam_calo_Z_allTrack.push_back(theXYZPoints[i].Z());
        }
        
        //Getting the SCE corrected start/end positions & directions
        std::sort(theXYZPoints.begin(), theXYZPoints.end(), [](auto a, auto b)
                  {return (a.Z() < b.Z());});
        if (theXYZPoints.size()) {
          reco_beam_calo_startX_allTrack = theXYZPoints[0].X();
          reco_beam_calo_startY_allTrack = theXYZPoints[0].Y();
          reco_beam_calo_startZ_allTrack = theXYZPoints[0].Z();
          reco_beam_calo_endX_allTrack = theXYZPoints.back().X();
          reco_beam_calo_endY_allTrack = theXYZPoints.back().Y();
          reco_beam_calo_endZ_allTrack = theXYZPoints.back().Z();
        }

        //New Calibration
        if (fRecalibrate){
          std::vector< float > new_dEdX = calibration_SCE.GetCalibratedCalorimetry( *pandora2Track, evt, "pandora2Track", fPandora2CaloSCE, 2);
          for( size_t i = 0; i < new_dEdX.size(); ++i ) {
            reco_beam_allTrack_calibrated_dEdX.push_back( new_dEdX[i] );
          }
          /*std::vector<double> new_dQdX = calibration_SCE.CalibratedQdX( *pandora2Track, evt, "pandora2Track", fCalorimetryTagSCE, 2, -10.);
           for (auto dqdx : new_dQdX)  reco_beam_calibrated_dQdX_SCE_allTrack.push_back(dqdx);*/
          
          /*std::vector<double> efield = calibration_SCE.GetEFieldVector( *pandora2Track, evt, "pandora2Track", fCalorimetryTagSCE, 2, -10.);
           for (auto ef : efield)  reco_beam_EField_SCE_allTrack.push_back(ef);*/
        }
        else{
          for (size_t i = 0; i < calo_dEdX.size(); ++i){
            reco_beam_allTrack_calibrated_dEdX.push_back(calo_dEdX[i]);
            //reco_beam_calibrated_dQdX_SCE_allTrack.push_back(calo_dQdX[i]);
            
            /*double E_field_nominal = detProp.Efield();   // Electric Field in the drift region in KV/cm
             geo::Vector_t E_field_offsets = {0., 0., 0.};
             if(sce->EnableCalEfieldSCE()&&fSCE) E_field_offsets = sce->GetCalEfieldOffsets(geo::Point_t{theXYZPoints[i].X(), theXYZPoints[i].Y(), theXYZPoints[i].Z()},calo[index].PlaneID().TPC);
             TVector3 E_field_vector = {E_field_nominal*(1 + E_field_offsets.X()), E_field_nominal*E_field_offsets.Y(), E_field_nominal*E_field_offsets.Z()};
             double E_field = E_field_vector.Mag();
             reco_beam_EField_SCE.push_back(E_field);*/
          }
        }
        double efield_placeholder = 1.;
        ////////////////////////////////////////////

        std::pair< double, int > pid_chi2_ndof = trackUtil.Chi2PID( reco_beam_allTrack_calibrated_dEdX, reco_beam_allTrack_resRange, templates[ 2212 ] );
        reco_beam_allTrack_Chi2_proton = pid_chi2_ndof.first;
        reco_beam_allTrack_Chi2_ndof = pid_chi2_ndof.second;
        
        std::vector< calo_point > reco_beam_calo_points;
        //Doing thin slice
        if (reco_beam_allTrack_calibrated_dEdX.size() &&
            reco_beam_allTrack_calibrated_dEdX.size() == reco_beam_TrkPitch_SCE_allTrack.size() &&
            reco_beam_allTrack_calibrated_dEdX.size() == reco_beam_calo_wire_allTrack.size()) {
          
          for( size_t i = 0; i < reco_beam_allTrack_calibrated_dEdX.size(); ++i ){
            reco_beam_calo_points.push_back(
                                            calo_point(reco_beam_calo_wire_allTrack[i], 1.,
                                                       reco_beam_TrkPitch_SCE_allTrack[i],
                                                       1., 1.,
                                                       1.,
                                                       1.,
                                                       reco_beam_allTrack_calibrated_dEdX[i],
                                                       reco_beam_allTrack_resRange[i], 0,
                                                       1., 0,
                                                       efield_placeholder, reco_beam_calo_X_allTrack[i],
                                                       reco_beam_calo_Y_allTrack[i], reco_beam_calo_Z_allTrack[i]));
          }
          
          //std::cout << "N Calo points: " << reco_beam_calo_points.size() << std::endl;
          //Sort
          std::sort( reco_beam_calo_points.begin(), reco_beam_calo_points.end(), [](calo_point a, calo_point b) {return ( a.z < b.z );} );
          
          //And also put these in the right order
          for( size_t i = 0; i < reco_beam_calo_points.size(); ++i ){
            calo_point thePoint = reco_beam_calo_points[i];
            reco_beam_calo_wire_allTrack[i] = thePoint.wire;
            reco_beam_TrkPitch_SCE_allTrack[i] = thePoint.pitch;
            reco_beam_allTrack_calibrated_dEdX[i] = thePoint.calibrated_dEdX;
            reco_beam_allTrack_resRange[i] = thePoint.res_range;
            reco_beam_calo_X_allTrack[i] = thePoint.x;
            reco_beam_calo_Y_allTrack[i] = thePoint.y;
            reco_beam_calo_Z_allTrack[i] = thePoint.z;
          }
          
          //Get the initial Energy KE
          //double mass = 0.;
          /*
          double init_KE = 0.;
          //std::cout << "Has BI? " << fMCHasBI << " " << evt.isRealData() << std::endl;
          if (evt.isRealData() || fMCHasBI) {
            double mass = 139.57;
            
            init_KE =  sqrt( 1.e6*beam_inst_P*beam_inst_P + mass*mass ) - mass;
            // std::cout << "MC has BI: " << init_KE << std::endl;
          }
          else{
            init_KE = sqrt(1.e6*true_beam_startP*true_beam_startP + true_beam_mass*true_beam_mass) - true_beam_mass;
          }
          
          reco_beam_incidentEnergies_allTrack.push_back( init_KE );
          for( size_t i = 0; i < reco_beam_calo_points.size() - 1; ++i ){ //-1 to not count the last slice
            //use dedx * pitch or new hit calculation?
            if (reco_beam_calo_points[i].calibrated_dEdX < 0.) continue;
            double this_energy = reco_beam_incidentEnergies_allTrack.back() - ( reco_beam_calo_points[i].calibrated_dEdX * reco_beam_calo_points[i].pitch );
            reco_beam_incidentEnergies_allTrack.push_back( this_energy );
          }
          if( reco_beam_incidentEnergies_allTrack.size() ) reco_beam_interactingEnergy_allTrack = reco_beam_incidentEnergies_allTrack.back();
          */
        }

      }
      else{
        reco_beam_allTrack_Chi2_proton = -999;
        reco_beam_allTrack_Chi2_ndof = -999;
      }
    }
    else{
      reco_beam_allTrack_ID = -999;
      reco_beam_allTrack_beam_cuts = -999;
      reco_beam_allTrack_startX = -999;
      reco_beam_allTrack_startY = -999;
      reco_beam_allTrack_startZ = -999;
      reco_beam_allTrack_endX = -999;
      reco_beam_allTrack_endY = -999;
      reco_beam_allTrack_endZ = -999;
      reco_beam_allTrack_Chi2_proton = -999;
      reco_beam_allTrack_Chi2_ndof = -999;
      reco_beam_allTrack_trackDirX    =  -999;
      reco_beam_allTrack_trackDirY    =  -999;
      reco_beam_allTrack_trackDirZ    =  -999;
      reco_beam_allTrack_trackEndDirX =  -999;
      reco_beam_allTrack_trackEndDirY =  -999;
      reco_beam_allTrack_trackEndDirZ =  -999;

    }
  }
  catch( const cet::exception &e ){
    MF_LOG_WARNING("PDSPAnalyzer") << "beam pandora2Track object not found, moving on" << std::endl;
  }
}

void pduneana::PDSPAnalyzer::G4RWGridWeights(
    std::vector<std::vector<G4ReweightTraj *>> & hierarchy,
    std::vector<fhicl::ParameterSet> & pars,
    std::vector<std::vector<double>> & weights,
    G4MultiReweighter * multi_rw) {
  std::vector<double> input(pars.size(), 1.);
  for (size_t i = 0; i < pars.size(); ++i) {
    weights.push_back(std::vector<double>());
    for (size_t j = 0; j < fGridPoints.size(); ++j) {
      input[i] = fGridPoints[j];
      bool set_values = multi_rw->SetAllParameterValues(input);

      if (set_values) {
        if (hierarchy.size()) {
          std::vector<G4ReweightTraj *> & init_trajs = hierarchy[0];
          weights.back().push_back(
              GetNTrajWeightFromSetPars(init_trajs, *multi_rw)); 
          for (size_t k = 1; k < hierarchy.size(); ++k) {
            std::vector<G4ReweightTraj *> & temp_trajs = hierarchy[k];
            weights.back().back()
                *= GetNTrajWeightFromSetPars(temp_trajs, *multi_rw);
          }
        }
        else {
          weights.back().push_back(1.);
        }
      }
      else {
        std::string message = "Could not Get N Traj Weight from set pars";
        throw std::runtime_error(message);
      }
    }

    //Reset to 1.
    input[i] = 1.;
  }
}

/////////////Below: To get the fitted lines for start/end directions
// define the parameteric line equation
void pduneana::line(double t, double *p,
                                 double &x, double &y, double &z) {
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

// calculate distance line-point
double pduneana::distance2(double x,double y,double z, double *p) {
   // distance line point is D= | (xp-x0) cross  ux |
   // where ux is direction of line and x0 is a point in the line (like t = 0)
   ROOT::Math::XYZVector xp(x,y,z);
   ROOT::Math::XYZVector x0(p[0], p[2], 0.);
   ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1.);
   ROOT::Math::XYZVector u = (x1-x0).Unit();
   double d2 = ((xp-x0).Cross(u)) .Mag2();
   return d2;
}

// function to be minimized
void pduneana::SumDistance2(int &, double *, double & sum,
                                         double * par, int) {
   // the TGraph must be a global variable
   TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
   assert(gr != 0);
   double * x = gr->GetX();
   double * y = gr->GetY();
   double * z = gr->GetZ();
   int npoints = gr->GetN();
   sum = 0;
   for (int i  = 0; i < npoints; ++i) {
      double d = distance2(x[i],y[i],z[i],par);
      sum += d;
   }
}

TVector3 pduneana::FitLine(const std::vector<TVector3> & input) {
  TGraph2D * gr = new TGraph2D();
  for (size_t i = 0; i < input.size(); ++i) {
    gr->SetPoint(i, input[i].X(), input[i].Y(), input[i].Z());
  }

  TVirtualFitter * min = TVirtualFitter::Fitter(0,4);
  min->SetObjectFit(gr);
  min->SetFCN(SumDistance2);

  double arglist[10];
  
  arglist[0] = -1;
  min->ExecuteCommand("SET PRINT",arglist,1);
  

  double pStart[4] = {1,1,1,1};
  min->SetParameter(0,"x0",pStart[0],0.01,0,0);
  min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
  min->SetParameter(2,"y0",pStart[2],0.01,0,0);
  min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

  arglist[0] = 1000; // number of function calls 
  arglist[1] = 0.001; // tolerance 
  min->ExecuteCommand("MIGRAD", arglist, 2);

  // get fit parameters
  double parFit[4];
  for (int i = 0; i < 4; ++i) {
     parFit[i] = min->GetParameter(i);
  }
  double startX1, startY1, startZ1;
  double startX2, startY2, startZ2;
  line(0, parFit, startX1, startY1, startZ1);
  line(1, parFit, startX2, startY2, startZ2);
  
  TVector3 diff(startX2 - startX1,
                startY2 - startY1,
                startZ2 - startZ1);
  delete gr;
  delete min;
  return diff;
}

void pduneana::PDSPAnalyzer::GetG4RWCoeffs(std::vector<double> & weights, std::vector<double> & coeffs) {
  //Check if there are no weights
  if (weights.empty()) return;

  TGraph gr(fGridPoints.size(), &fGridPoints[0], &weights[0]);
  gr.Fit("pol9", "Q");
  TF1 * fit_func = (TF1*)gr.GetFunction("pol9");
  for (int j = 0; j < fit_func->GetNpar(); ++j) {
    coeffs.push_back(fit_func->GetParameter(j));
  }
}

double pduneana::PDSPAnalyzer::GetG4RWExpCoeffs(std::vector<double> & weights, std::vector<double> & coeffs) {
  //Check if there are no weights
  if (weights.empty()) return -999.;

  TGraph gr(fGridPoints.size(), &fGridPoints[0], &weights[0]);
  TF1 fit_func("weight_func", "[0]*exp([1]*x)", fGridPoints[0], fGridPoints.back());
  TFitResultPtr result = gr.Fit(&fit_func, "QS");
  for (int j = 0; j < fit_func.GetNpar(); ++j) {
    coeffs.push_back(fit_func.GetParameter(j));
  }

  return result->Chi2();
}



std::vector<int> pduneana::PDSPAnalyzer::PrimaryHierarchy(
        int ID, const sim::ParticleList & plist, bool verbose) {

  std::deque<int> to_create = {ID};
  std::vector<int> results;

  while (to_create.size()) {
    auto part = plist[to_create[0]];
    results.push_back(to_create[0]);
    for (int i = 0; i < part->NumberDaughters(); ++i) {
      int daughter_ID = part->Daughter(i);
      auto d_part = plist[daughter_ID];
      if (abs(d_part->PdgCode()) != 11) {
        to_create.push_back(daughter_ID);
        if (verbose) {
          std::cout << "Adding daughter " << to_create.back() << " " <<
                       d_part->PdgCode() << std::endl;
        }
      }
    }
  
    to_create.pop_front();
  }

  return results;
}
DEFINE_ART_MODULE(pduneana::PDSPAnalyzer)
