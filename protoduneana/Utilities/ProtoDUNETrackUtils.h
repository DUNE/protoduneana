#ifndef PROTODUNE_TRACK_UTILS_H
#define PROTODUNE_TRACK_UTILS_H

///////////////////////////////////////////////////////////////
// ProtoDUNETrackUtils
//  - Class to help analysers access useful track information
// 
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "art/Framework/Principal/Event.h"
#include "TProfile.h"

namespace protoana {

  struct BrokenTrack {
    const recob::Track * firstTrack;
    const recob::Track * secondTrack;

    double CosTheta;
    std::vector< float > Combined_ResidualRange;
    std::vector< float > Combined_dQdx;
    std::vector< float > Combined_dEdx;

    bool Valid;

  };

  class ProtoDUNETrackUtils {

  public:

    ProtoDUNETrackUtils();
    ~ProtoDUNETrackUtils();

    /// Get the cosmic tag(s) from a given reco track
    std::vector<anab::CosmicTag> GetRecoTrackCosmicTag(const recob::Track &track, art::Event const &evt, std::string trackModule) const;
    /// Get the T0(s) from a given reco track
    std::vector<anab::T0> GetRecoTrackT0(const recob::Track &track, art::Event const &evt, std::string trackModule) const;
    /// Get the Calorimetry(s) from a given reco track
    std::vector<anab::Calorimetry> GetRecoTrackCalorimetry(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule) const;
    /// Calibrate a Calorimetry object for a given plane from a given track 
    std::vector< float > CalibrateCalorimetry(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule, const fhicl::ParameterSet & ps );
    float calc_dEdX( double, double, double, double, double, double );
    /// Get the hits from a given reco track
    const std::vector<const recob::Hit*> GetRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const;
    /// Get the hits from a given reco track from a specific plane
    const std::vector<const recob::Hit*> GetRecoTrackHitsFromPlane(const recob::Track &track, art::Event const &evt, const std::string trackModule, unsigned int planeID) const;
    /// Get the number of hits from a given reco track
    unsigned int GetNumberRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const;
    /// Get the PID from a given track
    std::vector<anab::ParticleID> GetRecoTrackPID(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string pidModule) const;
    /// Try to determine if it's a broken track
    BrokenTrack IsBrokenTrack( const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule, const fhicl::ParameterSet & BrokenTrackPars, const fhicl::ParameterSet & CalorimetryPars );
   /// Until we have fully calibrated calorimetry, use this PID algo
   std::pair< double, int > Chi2PIDFromTrack_MC( const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule, TProfile * profile );
   std::pair< double, int > Chi2PID( const std::vector< double > & track_dedx, const std::vector< double > & range, TProfile * profile ) const;

   // Get the hits from the track, organized by the Track Trajectory points
   //std::map< size_t, std::vector< const recob::Hit * > > GetRecoHitsFromTrajPoints( const recob::Track & track, art::Event const & evt, std::string trackModule );
   std::map< size_t, const recob::Hit * > GetRecoHitsFromTrajPoints( const recob::Track & track, art::Event const & evt, std::string trackModule );

   bool IsBeamlike( const recob::Track & track, art::Event const & evt, const fhicl::ParameterSet & BeamPars, bool flip=false );


  std::pair<double, int> GetVertexMichelScore(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    double min_length = 5., double min_x = -200.,
    double max_x = 0., double min_y = 200., double max_y = 500.,
    double min_z = 25., bool check_wire = true, double check_x = 0,
    double check_y = 0., double check_z = 0.);
  std::pair<double, double> GetVertexMichelScore_weight_by_charge(
      const recob::Track & track, const art::Event & evt,
      const std::string trackModule, const std::string hitModule,
      double min_length = 5., double min_x = -200.,
      double max_x = 0., double min_y = 200., double max_y = 500.,
      double min_z = 25., bool check_wire = true, double check_x = 0,
      double check_y = 0., double check_z = 0.);
  std::pair<double, int> GetVertexMichelScoreAlt(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    double min_length = 5., double min_x = -200.,
    double max_x = 0., double min_y = 200., double max_y = 500.,
    double min_z = 25., bool check_wire = true, double check_x = 0,
    double check_y = 0., double check_z = 0.);

  private:


  };

}

#endif

