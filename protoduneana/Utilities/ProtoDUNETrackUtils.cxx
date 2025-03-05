#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

#include "larcore/Geometry/WireReadout.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "lardata/ArtDataHelper/MVAReader.h"

#include "TFile.h"
#include "TH1F.h"

#include <string>

protoana::ProtoDUNETrackUtils::ProtoDUNETrackUtils(){

}

protoana::ProtoDUNETrackUtils::~ProtoDUNETrackUtils(){

}

std::vector<anab::CosmicTag> protoana::ProtoDUNETrackUtils::GetRecoTrackCosmicTag(const recob::Track &track, art::Event const &evt, std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  unsigned int trackIndex = track.ID();

  // Convert to std::vector<anab::CosmicTag> from std::vector<art::Ptr<anab::CosmicTag>>
  std::vector<anab::CosmicTag> trackTags;

  try{
    const art::FindManyP<anab::CosmicTag> findCosmicTags(recoTracks,evt,trackModule);
    for(unsigned int t = 0; t < findCosmicTags.at(trackIndex).size(); ++t){
      trackTags.push_back((*(findCosmicTags.at(trackIndex)[t])));
    }
  }
  catch(...){
//    std::cerr << "Product not found - returning empty vector" << std::endl;
  }

  return trackTags;    
}

std::vector<anab::T0> protoana::ProtoDUNETrackUtils::GetRecoTrackT0(const recob::Track &track, art::Event const &evt, std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  unsigned int trackIndex = track.ID();

  // Convert to std::vector<anab::T0> from std::vector<art::Ptr<anab::T0>>
  std::vector<anab::T0> trackT0s;
  
  try{
    const art::FindManyP<anab::T0> findTrackT0s(recoTracks,evt,trackModule);
    for(unsigned int t = 0; t < findTrackT0s.at(trackIndex).size(); ++t){
      trackT0s.push_back((*(findTrackT0s.at(trackIndex)[t])));
    }
  }
  catch(...){
//    std::cerr << "Product not found - returning empty vector" << std::endl;
  }
  
  return trackT0s;

}

// Get the Calorimetry(s) from a given reco track
std::vector<anab::Calorimetry> protoana::ProtoDUNETrackUtils::GetRecoTrackCalorimetry(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  std::vector<anab::Calorimetry> caloInfo;
  
  try{
    const art::FindManyP<anab::Calorimetry> findCalorimetry(recoTracks,evt,caloModule);
    std::vector<art::Ptr<anab::Calorimetry>> theseCalos = findCalorimetry.at(track.ID());

    for( auto calo : theseCalos){
      caloInfo.push_back(*calo);
    }
  }
  catch(...){
    std::cerr << "No track calorimetry object found... returning empty vector" << std::endl;
  }

  return caloInfo;
}

// Get the hits from a given reco track
const std::vector<const recob::Hit*> protoana::ProtoDUNETrackUtils::GetRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,trackModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(track.ID());

  std::vector<const recob::Hit*> trackHits;

  for(const art::Ptr<recob::Hit> hit : inputHits){

    trackHits.push_back(hit.get());

  }

  return trackHits;  

}

const std::vector<const recob::Hit*> protoana::ProtoDUNETrackUtils::GetRecoTrackHitsFromPlane(const recob::Track &track, art::Event const &evt, const std::string trackModule, unsigned int planeID ) const{

  std::vector<const recob::Hit*> trackHits;
  if( planeID > 2 ){
    std::cout << "Please input plane 0, 1, or 2" << std::endl;
    return trackHits;
  }

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,trackModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(track.ID());

  for(const art::Ptr<recob::Hit> hit : inputHits){
    unsigned int thePlane = hit.get()->WireID().asPlaneID().Plane;
    if( thePlane != planeID ) continue;
       
    trackHits.push_back(hit.get());

  }

  return trackHits;  

}

std::vector< float >  protoana::ProtoDUNETrackUtils::CalibrateCalorimetry(  const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule, const fhicl::ParameterSet &ps ) {


  unsigned int planeID = ps.get< unsigned int >( "PlaneID" );
   
  double betap  = ps.get< double >( "betap"  );
  double Rho    = ps.get< double >( "Rho"    );
  double Efield = ps.get< double >( "Efield" );
  double Wion   = ps.get< double >( "Wion"   );
  double alpha  = ps.get< double >( "alpha"  );
  double norm_factor = ps.get< double >( "norm_factor" );
  double calib_factor = ps.get< double >( "calib_factor" );
  std::string X_correction_name = ps.get< std::string >( "X_correction" );
  TFile X_correction_file = TFile( X_correction_name.c_str(), "OPEN" );
  TH1F * X_correction_hist = NULL;

  bool UseNewVersion = ps.get< bool >( "UseNewVersion", false );
  if( UseNewVersion ){
    std::string hist_name = "dqdx_X_correction_hist_" + std::to_string(planeID);
    X_correction_hist = (TH1F*)X_correction_file.Get( hist_name.c_str() );
    
  }
  else{
    X_correction_hist = (TH1F*)X_correction_file.Get( "dqdx_X_correction_hist" );
  }


  std::vector< float > calibrated_dEdx;

  //Get the Calorimetry vector from the track
  /*std::vector< anab::Calorimetry >*/ auto caloVector = GetRecoTrackCalorimetry( track, evt, trackModule, caloModule ); 
  
  size_t calo_position;
  bool found_plane = false;
  for( size_t i = 0; i < caloVector.size(); ++i ){
     unsigned int thePlane = caloVector.at(i).PlaneID().Plane;
     if( thePlane == planeID ){
       calo_position = i;
       found_plane = true;
       break;
     }
  }

  if( !found_plane ){
    std::cout << "Could not find the correct plane in the calorimetry vector" << std::endl;
    return calibrated_dEdx;
  }

  std::vector< float > dQdX = caloVector.at( calo_position).dQdx();
  auto theXYZPoints = caloVector.at( calo_position).XYZ();

  //Get the hits from the track from a specific plane
  const std::vector< const recob::Hit* > hits = GetRecoTrackHitsFromPlane( track, evt, trackModule, planeID ); 
  if( hits.size() == 0 ){
    std::cout << "Got empty hits vector" << std::endl;
    return calibrated_dEdx;
  }

  //Do Ajib's correction 
  for( size_t i = 0; i < dQdX.size(); ++i ){ 
    float hit_x = theXYZPoints[i].X();
    int X_bin = X_correction_hist->FindBin( hit_x );
    float X_correction = X_correction_hist->GetBinContent(X_bin);

    float corrected_dq_dx = dQdX[i] * X_correction * norm_factor;
    float scaled_corrected_dq_dx = corrected_dq_dx / calib_factor;
    float cal_de_dx = calc_dEdX( scaled_corrected_dq_dx,  betap,  Rho,  Efield,  Wion,  alpha );
 
    calibrated_dEdx.push_back( cal_de_dx );
  }


  return calibrated_dEdx;
}

float protoana::ProtoDUNETrackUtils::calc_dEdX(double dqdx, double betap, double Rho, double Efield, double Wion, double alpha){
  return (exp(dqdx*(betap/(Rho*Efield)*Wion))-alpha)/(betap/(Rho*Efield));  
}


// Get the hits from a given reco track
unsigned int protoana::ProtoDUNETrackUtils::GetNumberRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const{

  return GetRecoTrackHits(track,evt,trackModule).size();

}

// Get the PID from a given track
std::vector<anab::ParticleID> protoana::ProtoDUNETrackUtils::GetRecoTrackPID(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string pidModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  std::vector<anab::ParticleID> pidvec;

  try{
    const art::FindManyP<anab::ParticleID> findPID(recoTracks,evt,pidModule);
    std::vector<art::Ptr<anab::ParticleID>> thePID = findPID.at(track.ID());
  
    for( auto pid : thePID){
      pidvec.push_back(*pid);
    }
  }
  catch(...){
    std::cerr << "No track PID object found... returning empty vector" << std::endl;
  }

  return pidvec;

}

protoana::BrokenTrack protoana::ProtoDUNETrackUtils::IsBrokenTrack( const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule, const fhicl::ParameterSet & BrokenTrackPars, const fhicl::ParameterSet & CalorimetryPars ){
  
  BrokenTrack theBrokenTrack;
  theBrokenTrack.Valid = false;
  
  double fBrokenTrackZ_low  = BrokenTrackPars.get< double >( "BrokenTrackZ_low" );
  double fBrokenTrackZ_high = BrokenTrackPars.get< double >( "BrokenTrackZ_high" );

  double fStitchTrackZ_low  = BrokenTrackPars.get< double >( "StitchTrackZ_low" );
  double fStitchTrackZ_high = BrokenTrackPars.get< double >( "StitchTrackZ_high" );
  
  double fStitchXTol = BrokenTrackPars.get< double >( "StitchXTol" );
  double fStitchYTol = BrokenTrackPars.get< double >( "StitchYTol" );

  ////Check the end of the track
  double endZ = track.Trajectory().End().Z();
  double endX = track.Trajectory().End().X();
  double endY = track.Trajectory().End().Y();
  if( fBrokenTrackZ_low < endZ && endZ < fBrokenTrackZ_high ){

    const auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
    for( auto const & tr : *recoTracks ){          

      //Skip the track in question 
      if( tr.ID() == track.ID() ) continue;

      double stitchStartZ = tr.Trajectory().Start().Z();
      if( fStitchTrackZ_low < stitchStartZ && stitchStartZ < fStitchTrackZ_high ){
        double deltaX = fabs(endX - tr.Trajectory().Start().X());
        double deltaY = fabs(endY - tr.Trajectory().Start().Y());

        std::cout << "Possible stitching track: " << stitchStartZ << " " << deltaX << " " << deltaY << std::endl;

        if( deltaX < fStitchXTol && deltaY < fStitchYTol ){


          //Get the cosine of the angle between them
          auto stitchDir = tr.StartDirection();
          double stitch_cos_theta =  stitchDir.X()*track.EndDirection().X() + stitchDir.Y()*track.EndDirection().Y() + stitchDir.Z()*track.EndDirection().Z() ;
          std::cout << "Cos_theta " << stitch_cos_theta << std::endl;


          unsigned int planeID = CalorimetryPars.get< unsigned int >( "PlaneID");
          //Get the calorimetries, calibrate, and combine
          /*std::vector< anab::Calorimetry >*/ auto broken_calos = GetRecoTrackCalorimetry(track, evt, trackModule, caloModule);

          bool found_broken_calo = false; 
          size_t calo_position= 0;
          for( size_t i = 0; i < broken_calos.size(); ++i ){
            unsigned int thePlane = broken_calos.at(i).PlaneID().Plane;
            if( thePlane == planeID ){
              found_broken_calo = true;
              calo_position = i;
              break;
            }
          }
          
          //If no calorimetry object is found for this track, return the default BrokenTrack.
          if( !found_broken_calo ) return theBrokenTrack; 

          auto broken_range = broken_calos.at( calo_position ).ResidualRange();
          auto broken_dQdx  = broken_calos.at( calo_position ).dQdx();
          std::vector< float > broken_cal_dEdx = CalibrateCalorimetry(  track, evt, trackModule, caloModule, CalorimetryPars );

          calo_position = 0;
          /*std::vector< anab::Calorimetry >*/ auto stitch_calos = GetRecoTrackCalorimetry(tr, evt, trackModule, caloModule);

          bool found_stitch_calo = false;
          for( size_t i = 0; i < stitch_calos.size(); ++i ){
            unsigned int thePlane = stitch_calos.at(i).PlaneID().Plane;
            if( thePlane == planeID ){
              found_stitch_calo = true;
              calo_position = i;
              break;
            }
          }

          //If no calorimetry object is found for this track, try the next 
          if( !found_stitch_calo ) continue; 

          auto stitch_range = stitch_calos.at( calo_position ).ResidualRange();
          auto stitch_dQdx  = stitch_calos.at( calo_position ).dQdx();
          std::vector< float > stitch_cal_dEdx = CalibrateCalorimetry(  tr, evt, trackModule, caloModule, CalorimetryPars );

          //piece them together in order       
          std::vector< float > combined_range, combined_dQdx, combined_dEdx;
          

          combined_range = stitch_range;
          if( stitch_range[0] > stitch_range.back() ){
            std::cout << "Adding range: " << stitch_range[0] << std::endl;
            for( size_t i = 0 ; i < broken_range.size(); ++i ){
              combined_range.push_back( broken_range[i] + stitch_range[0] );
            }
          }
          else{
            std::cout << "Adding range: " << stitch_range[0] << std::endl;
            for( size_t i = 0 ; i < broken_range.size(); ++i ){
              combined_range.push_back( broken_range[i] + stitch_range.back() );
            }
          }

          for( size_t i = 0; i < combined_range.size(); ++i ){
            std::cout << combined_range[i] << std::endl;
          }

          combined_dQdx = stitch_dQdx;
          combined_dQdx.insert( combined_dQdx.end(), broken_dQdx.begin(), broken_dQdx.end() );
          combined_dEdx = stitch_cal_dEdx;
          combined_dEdx.insert( combined_dEdx.end(), broken_cal_dEdx.begin(), broken_cal_dEdx.end() );
/*          float total_stitch_range;
          if( stitch_range[0] > stitch_range.back() ){
            total_stitch_range = stitch_range[0]; 
          }
          else{
            total_stitch_range = stitch_range.back();
          }

          std::cout << stitch_range[0] << " " << stitch_range.back() << " " << total_stitch_range << std::endl;

          if( broken_range[0] > broken_range.back() ){
            for( size_t i = 0; i < broken_range.size(); ++i ){
              combined_range.push_back( broken_range.at(i) + total_stitch_range );
              combined_dQdx.push_back(  broken_dQdx.at(i) );
              combined_dEdx.push_back(  broken_cal_dEdx.at(i) );
            }
          }
          else{
            for( size_t i = broken_range.size() - 1; i >= 0; --i ){
              combined_range.push_back( broken_range.at(i) + total_stitch_range );
              combined_dQdx.push_back(  broken_dQdx.at(i) );
              combined_dEdx.push_back(  broken_cal_dEdx.at(i) );
            }
          }

          if( stitch_range[0] > stitch_range.back() ){
            combined_range.insert( combined_range.end(), stitch_range.begin(), stitch_range.end() );
            combined_dQdx.insert( combined_dQdx.end(), stitch_dQdx.begin(), stitch_dQdx.end() );
            combined_dEdx.insert( combined_dEdx.end(), stitch_cal_dEdx.begin(), stitch_cal_dEdx.end() );
          }
          else{
            combined_range.insert( combined_range.end(), stitch_range.rbegin(), stitch_range.rend() );
            combined_dQdx.insert( combined_dQdx.end(), stitch_dQdx.rbegin(), stitch_dQdx.rend() );
            combined_dEdx.insert( combined_dEdx.end(), stitch_cal_dEdx.rbegin(), stitch_cal_dEdx.rend() );
          }
*/          

          theBrokenTrack.firstTrack = &track;
          theBrokenTrack.secondTrack = &tr;
          theBrokenTrack.CosTheta = stitch_cos_theta; 
          theBrokenTrack.Combined_ResidualRange = combined_range;
          theBrokenTrack.Combined_dQdx = combined_dQdx;
          theBrokenTrack.Combined_dEdx = combined_dEdx;
          theBrokenTrack.Valid = true;

          return theBrokenTrack;
        }
      }
    }
  }
  return theBrokenTrack;
}


std::pair< double, int > protoana::ProtoDUNETrackUtils::Chi2PIDFromTrack_MC( const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule, TProfile * profile ){

  /*std::vector< anab::Calorimetry >*/ auto calo = GetRecoTrackCalorimetry(track, evt, trackModule, caloModule);
  std::vector< double > calo_dEdX;
  std::vector< double > calo_range;
  for( size_t i = 0; i < calo[0].dEdx().size(); ++i ){
    calo_dEdX.push_back( calo[0].dEdx()[i] );
    calo_range.push_back( calo[0].ResidualRange()[i] );
  }

  return Chi2PID( calo_dEdX, calo_range, profile );

}

std::pair< double, int > protoana::ProtoDUNETrackUtils::Chi2PID( const std::vector< double > & track_dedx, const std::vector< double > & range, TProfile * profile ){

  double pid_chi2 = 0.; 
  int npt = 0;

  if( track_dedx.size() < 1 || range.size() < 1 )
    return std::make_pair(9999., -1);

  //Ignore first and last point
  for( size_t i = 1; i < track_dedx.size()-1; ++i ){

    //Skip large pulse heights
    if( track_dedx[i] > 1000. || track_dedx[i] < 0.) {
      continue;
    }

    int bin = profile->FindBin( range[i] );
    if( bin >= 1 && bin <= profile->GetNbinsX() ){

      double template_dedx = profile->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( profile->GetBinContent( bin - 1 ) + profile->GetBinContent( bin + 1 ) ) / 2.;        
      }

      double template_dedx_err = profile->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( profile->GetBinError( bin - 1 ) + profile->GetBinError( bin + 1 ) ) / 2.;        
      }


      double dedx_res = 0.04231 + 0.0001783 * track_dedx[i] * track_dedx[i];
      dedx_res *= track_dedx[i]; 

      //Chi2 += ( track_dedx - template_dedx )^2  / ( (template_dedx_err)^2 + (dedx_res)^2 )
      pid_chi2 += ( pow( (track_dedx[i] - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) ); 

      ++npt;
    }
  }

  if( npt == 0 )
    return std::make_pair(9999., -1);
  
    

  return std::make_pair(pid_chi2, npt); 
}

//std::map< size_t, std::vector< const recob::Hit * > > protoana::ProtoDUNETrackUtils::GetRecoHitsFromTrajPoints(const recob::Track & track, art::Event const & evt, std::string trackModule){
std::map< size_t, const recob::Hit * > protoana::ProtoDUNETrackUtils::GetRecoHitsFromTrajPoints(const recob::Track & track, art::Event const & evt, std::string trackModule){

   auto recoTracks = evt.getValidHandle< std::vector< recob::Track > >(trackModule);
   art::FindManyP< recob::Hit, recob::TrackHitMeta >  trackHitMetas(recoTracks,evt,trackModule);
   art::FindManyP< recob::Hit > findHits(recoTracks,evt,trackModule);

   std::vector< art::Ptr< recob::Hit > > track_hits = findHits.at( track.ID() );


   //First, find the location of the beam track in the track list
   size_t beam_index = 0;
   for( size_t i = 0; i < recoTracks->size(); ++i ){
     if( (*recoTracks)[i].ID() == track.ID() ){
       beam_index = i;
       break;
     }        
   }


   //std::map< size_t, std::vector< const recob::Hit * > > results;
   std::map< size_t, const recob::Hit * > results;
   if( trackHitMetas.isValid() ){

     auto beamHits  = trackHitMetas.at( beam_index );
     auto beamMetas = trackHitMetas.data( beam_index );    

     for( size_t i = 0; i < beamHits.size(); ++i ){

       if( beamMetas[i]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max()) )
         continue;

       if( !track.HasValidPoint( beamMetas[i]->Index() ) ){
         std::cout << "Has no valid hit: " << beamMetas[i]->Index() << std::endl;
         continue;
       }

       //results[ beamMetas[i]->Index() ] = std::vector< const recob::Hit * >();

       for( size_t j = 0; j < track_hits.size(); ++j ){

         //if( track_hits[j]->WireID().Plane == 2 ){//Look at just the collection plane

           if( beamHits[i].key() == track_hits[j].key() ){

             if( beamMetas[i]->Index() >= track.NumberTrajectoryPoints() ){
               throw cet::exception("ProtoDUNETrackUtils.cxx") 
                     << "Requested track trajectory index " << beamMetas[i]->Index() 
                     << " exceeds the total number of trajectory points "<< track.NumberTrajectoryPoints() 
                     << " for track index " << beam_index 
                     << ". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov";
             }

             //If we've reached here, it's a good hit within the track. Connect to the trajectory point
             //results[ beamMetas[i]->Index() ].push_back( track_hits[j].get() ); 
             results[ beamMetas[i]->Index() ] = track_hits[j].get(); 
           }
         //}
       }
     }

   }
   return results;
   
}

bool protoana::ProtoDUNETrackUtils::IsBeamlike( const recob::Track & track, art::Event const & evt, const fhicl::ParameterSet & BeamPars, bool flip ){

   double startX = track.Trajectory().Start().X();
   double startY = track.Trajectory().Start().Y();
   double startZ = track.Trajectory().Start().Z();

   double endX = track.Trajectory().End().X();
   double endY = track.Trajectory().End().Y();
   double endZ = track.Trajectory().End().Z();

   auto startDir = track.StartDirection();
   auto endDir   = track.EndDirection();
   double trackDirX = 0.;
   double trackDirY = 0.;
   double trackDirZ = 0.;

   //'Flip' the track if endZ < startZ
   if( flip && ( endZ < startZ ) ){
     startX = endX;
     startY = endY;
     startZ = endZ;

     trackDirX = -1. * endDir.X();
     trackDirY = -1. * endDir.Y();
     trackDirZ = -1. * endDir.Z();
   }
   else{
     trackDirX = startDir.X();
     trackDirY = startDir.Y();
     trackDirZ = startDir.Z();
   }

   double beamX = 0.;
   double beamY = 0.;

   double beamDirX = 0.;
   double beamDirY = 0.;
   double beamDirZ = 0.;
   
   //Real data: compare to the reconstructed beam particle
   if( evt.isRealData() ){
     std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
     auto beamHandle = evt.getValidHandle< std::vector< beam::ProtoDUNEBeamEvent > >("beamevent");
                           
     if( beamHandle.isValid()){
       art::fill_ptr_vector(beamVec, beamHandle);
     }
     //Should just have one
     const beam::ProtoDUNEBeamEvent & beamEvent = *(beamVec.at(0));

     const std::vector< recob::Track > & beamTracks = beamEvent.GetBeamTracks();
     if( beamTracks.size() == 0 ){
       std::cout << "Warning: no tracks associated to beam data" << std::endl;
       return false;  
     }
     else if( beamTracks.size() > 1 ){
       std::cout << "Warning: mutiple tracks associated to beam data" << std::endl;
       return false;  
     }
     
     beamX = beamTracks.at(0).Trajectory().End().X();
     beamY = beamTracks.at(0).Trajectory().End().Y();

     beamDirX = beamTracks.at(0).EndDirection().X();
     beamDirY = beamTracks.at(0).EndDirection().Y();
     beamDirZ = beamTracks.at(0).EndDirection().Z();

   }
   //MC: compare to the projected particle from the particle gun
   else{
     protoana::ProtoDUNETruthUtils truthUtil;
     auto mcTruths = evt.getValidHandle< std::vector< simb::MCTruth > >("generator");
                         
     const simb::MCParticle* true_beam_particle = truthUtil.GetGeantGoodParticle((*mcTruths)[0],evt);

     if( !true_beam_particle ){
       std::cout << "No true beam particle" << std::endl;       
       return false;
     }


     beamDirX = true_beam_particle->Px() / true_beam_particle->P();
     beamDirY = true_beam_particle->Py() / true_beam_particle->P();
     beamDirZ = true_beam_particle->Pz() / true_beam_particle->P();

     //Project the beam to Z = 0
     beamX = true_beam_particle->Position(0).X() + (-1.*true_beam_particle->Position(0).Z())*(beamDirX / beamDirZ);
     beamY = true_beam_particle->Position(0).Y() + (-1.*true_beam_particle->Position(0).Z())*(beamDirY / beamDirZ);

   }

   double deltaX = startX - beamX;
   double deltaY = startY - beamY;

   double costheta = beamDirX*trackDirX + beamDirY*trackDirY + beamDirZ*trackDirZ;

   std::pair< double, double > startX_cut = BeamPars.get< std::pair< double, double > >("TrackStartXCut");
   std::pair< double, double > startY_cut = BeamPars.get< std::pair< double, double > >("TrackStartYCut");
   std::pair< double, double > startZ_cut = BeamPars.get< std::pair< double, double > >("TrackStartZCut");
   double costheta_cut = BeamPars.get< double >("TrackDirCut");

   if( deltaX < startX_cut.first || deltaX > startX_cut.second )
     return false;
   if( deltaY < startY_cut.first || deltaY > startY_cut.second )
     return false;
   if( startZ < startZ_cut.first || startZ > startZ_cut.second )
     return false;
   if( costheta < costheta_cut )
     return false;   

   //If here, the track is in the good beam region
   return true;
}

//Jake's implementation
std::pair<double, int> protoana::ProtoDUNETrackUtils::GetVertexMichelScore(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    double min_length, double min_x,
    double max_x, double min_y, double max_y, double min_z, bool check_wire,
    double check_x, double check_y, double check_z) {

  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  anab::MVAReader<recob::Hit, 4> hitResults(evt, "emtrkmichelid:emtrkmichel");

  //Skip short tracks
  //Skip tracks that start/end outside of interesting volume
  auto start = track.Vertex();
  auto end = track.End();
  if ((TMath::Max(start.X(), end.X()) > max_x) ||
      (TMath::Min(start.X(), end.X()) < min_x) ||
      (TMath::Max(start.Y(), end.Y()) > max_y) ||
      (TMath::Min(start.Y(), end.Y()) < min_y) ||
      (TMath::Min(start.Z(), end.Z()) < min_z) ||
      (track.Length() < min_length)) {
    return {-1., 0};
  }


  //Get the hits from the TrackHitMetas and only for view 2 (collection plane)
  std::map<size_t, const recob::Hit *>
      hits_from_traj = GetRecoHitsFromTrajPoints(track, evt, trackModule);
  std::vector<const recob::Hit *> hits_from_traj_view2;
  std::vector<size_t> index_from_traj_view2;

  for (auto it = hits_from_traj.begin(); it != hits_from_traj.end(); ++it) {
    if (it->second->View() != 2) continue;
    hits_from_traj_view2.push_back(it->second); 
    index_from_traj_view2.push_back(it->first);
  }

  //Find the vertex hit & info to compare to later
  double highest_z = -100.; 
  int vertex_tpc = -1;
  int vertex_wire = -1;
  float vertex_peak_time = -1.;

  if (check_z) {
    for (const auto * hit : hits_from_traj_view2) {
      double wire_z = wireReadout.Wire(hit->WireID()).GetCenter().Z();
      if (wire_z > highest_z) {
        highest_z = wire_z;
        vertex_tpc = hit->WireID().TPC;
        vertex_peak_time = hit->PeakTime();
        vertex_wire = hit->WireID().Wire;
      }
    }
  }
  else {
    const recob::TrackTrajectory & traj = track.Trajectory();
    double highest_diff = -1.;
    for (size_t i = 0; i < hits_from_traj_view2.size(); ++i) {
      const recob::Hit * hit = hits_from_traj_view2[i];
      size_t traj_index = index_from_traj_view2[i];
      
      double traj_x = traj.LocationAtPoint(traj_index).X();
      double traj_y = traj.LocationAtPoint(traj_index).Y();
      double traj_z = traj.LocationAtPoint(traj_index).Z();

      double diff = sqrt(std::pow((traj_x - check_x), 2) + 
                         std::pow((traj_y - check_y), 2) + 
                         std::pow((traj_z - check_z), 2));
      if (diff > highest_diff) {
        highest_diff = diff;
        vertex_tpc = hit->WireID().TPC;
        vertex_peak_time = hit->PeakTime();
        vertex_wire = hit->WireID().Wire;
      }
    }
  }
  
  std::pair<double, int> results = {0., 0};

  //Go through all hits in the event.
  auto allHits = evt.getValidHandle<std::vector<recob::Hit>>(hitModule);
  std::vector<art::Ptr<recob::Hit>> hit_vector;
  art::fill_ptr_vector(hit_vector, allHits);
  art::FindManyP<recob::Track> tracks_from_all_hits(allHits, evt, trackModule);
  for (size_t i = 0; i < hit_vector.size(); ++i) {

    //If this hit is in the trajectory hits vector, skip
    const recob::Hit * theHit = hit_vector[i].get();
    if (std::find(hits_from_traj_view2.begin(),
        hits_from_traj_view2.end(),
        theHit) != hits_from_traj_view2.end()) {
      continue;
    }

    //Skip hits that are outside of our TPC/plane or window of interest
    int wire = theHit->WireID().Wire;
    float peak_time = theHit->PeakTime();
    int tpc = theHit->WireID().TPC; 
    int plane = theHit->View();
    if ((abs(wire - vertex_wire) > 15) ||
        (abs(peak_time - vertex_peak_time) > 100.) ||
        (tpc != vertex_tpc) || (plane != 2)) {
      continue;
    }

    //It's ok if the hits don't come from any track
    //or if that track is the primary one, because sometimes the michel hits
    //are associated to it.
    //
    //It's not ok if the hits come from another track that is long
    //(i.e. an actual daughter). Skip these
    auto & tracks_from_hit = tracks_from_all_hits.at(hit_vector[i].key());
    if (!tracks_from_hit.empty() &&
        (tracks_from_hit[0].get()->ID() != track.ID()) &&
        (tracks_from_hit[0].get()->Length() > 25.))
      continue;

    //add up the CNN results 
    std::array<float, 4> cnn_out = hitResults.getOutput(hit_vector[i]);
    results.first += cnn_out[hitResults.getIndex("michel")];
    results.second += 1;
  }

  return results;
}
//Yinrui's modification
std::pair<double, double> protoana::ProtoDUNETrackUtils::GetVertexMichelScore_weight_by_charge(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    double min_length, double min_x,
    double max_x, double min_y, double max_y, double min_z, bool check_wire,
    double check_x, double check_y, double check_z) {

  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  anab::MVAReader<recob::Hit, 4> hitResults(evt, "emtrkmichelid:emtrkmichel");

  //Skip short tracks
  //Skip tracks that start/end outside of interesting volume
  auto start = track.Vertex();
  auto end = track.End();
  if ((TMath::Max(start.X(), end.X()) > max_x) ||
      (TMath::Min(start.X(), end.X()) < min_x) ||
      (TMath::Max(start.Y(), end.Y()) > max_y) ||
      (TMath::Min(start.Y(), end.Y()) < min_y) ||
      (TMath::Min(start.Z(), end.Z()) < min_z) ||
      (track.Length() < min_length)) {
    return {-1., 0};
  }


  //Get the hits from the TrackHitMetas and only for view 2 (collection plane)
  std::map<size_t, const recob::Hit *>
      hits_from_traj = GetRecoHitsFromTrajPoints(track, evt, trackModule);
  std::vector<const recob::Hit *> hits_from_traj_view2;
  std::vector<size_t> index_from_traj_view2;

  for (auto it = hits_from_traj.begin(); it != hits_from_traj.end(); ++it) {
    if (it->second->View() != 2) continue;
    hits_from_traj_view2.push_back(it->second);
    index_from_traj_view2.push_back(it->first);
  }

  //Find the vertex hit & info to compare to later
  double highest_z = -100.;
  int vertex_tpc = -1;
  int vertex_wire = -1;
  float vertex_peak_time = -1.;

  if (check_z) {
    for (const auto * hit : hits_from_traj_view2) {
      double wire_z = wireReadout.Wire(hit->WireID()).GetCenter().Z();
      if (wire_z > highest_z) {
        highest_z = wire_z;
        vertex_tpc = hit->WireID().TPC;
        vertex_peak_time = hit->PeakTime();
        vertex_wire = hit->WireID().Wire;
      }
    }
  }
  else {
    const recob::TrackTrajectory & traj = track.Trajectory();
    double highest_diff = -1.;
    for (size_t i = 0; i < hits_from_traj_view2.size(); ++i) {
      const recob::Hit * hit = hits_from_traj_view2[i];
      size_t traj_index = index_from_traj_view2[i];
      
      double traj_x = traj.LocationAtPoint(traj_index).X();
      double traj_y = traj.LocationAtPoint(traj_index).Y();
      double traj_z = traj.LocationAtPoint(traj_index).Z();

      double diff = sqrt(std::pow((traj_x - check_x), 2) +
                         std::pow((traj_y - check_y), 2) +
                         std::pow((traj_z - check_z), 2));
      if (diff > highest_diff) {
        highest_diff = diff;
        vertex_tpc = hit->WireID().TPC;
        vertex_peak_time = hit->PeakTime();
        vertex_wire = hit->WireID().Wire;
      }
    }
  }
  
  std::pair<double, double> results = {0., 0.};

  //Go through all hits in the event.
  auto allHits = evt.getValidHandle<std::vector<recob::Hit>>(hitModule);
  std::vector<art::Ptr<recob::Hit>> hit_vector;
  art::fill_ptr_vector(hit_vector, allHits);
  art::FindManyP<recob::Track> tracks_from_all_hits(allHits, evt, trackModule);
  for (size_t i = 0; i < hit_vector.size(); ++i) {

    //If this hit is in the trajectory hits vector, skip
    const recob::Hit * theHit = hit_vector[i].get();
    if (std::find(hits_from_traj_view2.begin(),
        hits_from_traj_view2.end(),
        theHit) != hits_from_traj_view2.end()) {
      continue;
    }

    //Skip hits that are outside of our TPC/plane or window of interest
    int wire = theHit->WireID().Wire;
    float peak_time = theHit->PeakTime();
    int tpc = theHit->WireID().TPC;
    int plane = theHit->View();
    if ((abs(wire - vertex_wire) > 15) ||
        (abs(peak_time - vertex_peak_time) > 100.) ||
        (tpc != vertex_tpc) || (plane != 2)) {
      continue;
    }

    //It's ok if the hits don't come from any track
    //or if that track is the primary one, because sometimes the michel hits
    //are associated to it.
    //
    //It's not ok if the hits come from another track that is long
    //(i.e. an actual daughter). Skip these
    auto & tracks_from_hit = tracks_from_all_hits.at(hit_vector[i].key());
    if (!tracks_from_hit.empty() &&
        (tracks_from_hit[0].get()->ID() != track.ID()) &&
        (tracks_from_hit[0].get()->Length() > 25.))
      continue;

    //add up the CNN results
    std::array<float, 4> cnn_out = hitResults.getOutput(hit_vector[i]);
    double hitcharge = hit_vector[i]->Integral();
    results.first += hitcharge*cnn_out[hitResults.getIndex("michel")];
    results.second += hitcharge;
  }

  return results;
}
//Ajib's implementation
std::pair<double, int> protoana::ProtoDUNETrackUtils::GetVertexMichelScoreAlt(
    const recob::Track & track, const art::Event & evt,
    const std::string trackModule, const std::string hitModule,
    double min_length, double min_x,
    double max_x, double min_y, double max_y, double min_z, bool check_wire,
    double check_x, double check_y, double check_z) {

  // Get all tracks
  std::vector < art::Ptr < recob::Track > > trkList;
  auto trkListHandle = evt.getHandle < std::vector < recob::Track > >(trackModule);
  if (trkListHandle) {
    art::fill_ptr_vector(trkList, trkListHandle);
  }

  // Get all hits
  std::vector < art::Ptr < recob::Hit > > hitList;
  auto hitListHandle = evt.getHandle < std::vector < recob::Hit > >(hitModule);
  if (hitListHandle) {
    art::fill_ptr_vector(hitList, hitListHandle);
  }

  // Get track-hit association
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trkListHandle, evt,trackModule); // to associate tracks and hits

  // Get hit-track association
  art::FindManyP<recob::Track> thass(hitListHandle, evt, trackModule); //to associate hit just trying

  // Get CNN scores
  anab::MVAReader<recob::Hit, 4> hitResults(evt, "emtrkmichelid:emtrkmichel");

  int endwire = -1;
  int endtpc = -1;
  double endpeakt = -1;
  std::vector<int> wirekeys;

  if (fmthm.isValid()){
    float zlast0=-99999;
    auto vhit=fmthm.at(track.ID());
    auto vmeta=fmthm.data(track.ID());
    for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
      bool fBadhit = false;
      if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
        fBadhit = true;
        continue;
      }
      if (vmeta[ii]->Index()>=track.NumberTrajectoryPoints()){
        throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<track.NumberTrajectoryPoints()<<" for track index "<<track.ID()<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
      }
      if (!track.HasValidPoint(vmeta[ii]->Index())){
        fBadhit = true;
        continue;
      }
      auto loc = track.LocationAtPoint(vmeta[ii]->Index());
      if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
      if (loc.Z()<-100) continue; //hit not on track
      if(vhit[ii]->WireID().Plane==2){
        wirekeys.push_back(vhit[ii].key());
        float zlast=loc.Z();
        if(zlast>zlast0){
          zlast0=zlast;
          endwire=vhit[ii]->WireID().Wire;
          endpeakt=vhit[ii]->PeakTime();
          endtpc=vhit[ii]->WireID().TPC;
        }
      }
    }
  }

  int ndaughterhits = 0;
  double average_daughter_score_mic = 0;
  
  for(size_t hitl=0;hitl<hitList.size();hitl++){
    std::array<float,4> cnn_out=hitResults.getOutput(hitList[hitl]);
    auto & tracks = thass.at(hitList[hitl].key());
    // hit not on the track
    if (std::find(wirekeys.begin(), wirekeys.end(), hitl) != wirekeys.end()) continue;
    // hit not on a long track
    if (!tracks.empty() && int(tracks[0].key()) != track.ID() && trkList[tracks[0].key()]->Length()>25) continue;
    int planeid=hitList[hitl]->WireID().Plane;
    if (planeid!=2) continue;
    int tpcid=hitList[hitl]->WireID().TPC;
    if (tpcid!=endtpc) continue;
    float peakth1=hitList[hitl]->PeakTime();
    int wireh1=hitList[hitl]->WireID().Wire;
    if(std::abs(wireh1-endwire)<8 && std::abs(peakth1-endpeakt)<150 && tpcid==endtpc){
      ++ndaughterhits;
      average_daughter_score_mic += cnn_out[hitResults.getIndex("michel")];
      //std::cout<<hitList[hitl]->WireID().Wire<<" "<<hitList[hitl]->PeakTime()<<" "<<hitList[hitl]->Integral()<<" "<<cnn_out[hitResults.getIndex("michel")]<<std::endl;
    }
  }

  if (ndaughterhits) average_daughter_score_mic /= ndaughterhits;
  
  return {average_daughter_score_mic, ndaughterhits};
}
