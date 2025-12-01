#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "lardataalg/Utilities/TrackingTypeUtils.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>
#include "TVector3.h"

protoana::ProtoDUNEBeamlineUtils::ProtoDUNEBeamlineUtils(fhicl::ParameterSet const& p){
  this->reconfigure(p);
}

protoana::ProtoDUNEBeamlineUtils::~ProtoDUNEBeamlineUtils(){

}

const beam::ProtoDUNEBeamEvent& protoana::ProtoDUNEBeamlineUtils::GetBeamEvent(art::Event const & evt) const{
  //Should just have one
  return evt.getProduct< std::vector< beam::ProtoDUNEBeamEvent> >(fBeamEventTag).front();
}

void protoana::ProtoDUNEBeamlineUtils::GetCurrent( art::Event const & evt ){
  Current = GetBeamEvent(evt).GetMagnetCurrent();
}

void protoana::ProtoDUNEBeamlineUtils::GetFibers( art::Event const & evt ){
 
  const beam::ProtoDUNEBeamEvent & beamEvent = GetBeamEvent(evt);

  for( size_t i = 0; i < AllDevices.size(); ++i ){
    std::string monitor = AllDevices[i];
    ActiveFibers[ monitor ] = beamEvent.GetActiveFibers( monitor );

    
    if( ActiveFibers.at(monitor).size() > 0){
      std::cout << monitor << " has active fibers: " << std::endl;
      for( size_t j = 0; j < ActiveFibers.at(monitor).size(); ++j ){
        std::cout << ActiveFibers.at(monitor).at(j) << " ";    
      }
      std::cout << std::endl << std::endl;
    }
    else{ 
      std::cout << monitor << " has no active fibers" << std::endl;
    }
  }
}

std::vector< recob::Track > protoana::ProtoDUNEBeamlineUtils::MakeTracks( art::Event const & evt ){

  std::vector< recob::Track > tracks;
 
  //Load fibers 
  GetFibers( evt ); 

  std::vector<short> HorizUpstreamFibers   = ActiveFibers.at(HorizUpstream);
  std::vector<short> HorizDownstreamFibers = ActiveFibers.at(HorizDownstream);
  std::vector<short> VertUpstreamFibers    = ActiveFibers.at(VertUpstream);
  std::vector<short> VertDownstreamFibers  = ActiveFibers.at(VertDownstream);

  //Check that the monitors all have at least one hit each
  bool all_good ;
  all_good = ( ( HorizUpstreamFibers.size() > 0 )   && ( VertUpstreamFibers.size() > 0 ) 
            && ( HorizDownstreamFibers.size() > 0 ) && ( VertDownstreamFibers.size() > 0 ) );

  if( !all_good ){
    std::cout << "At least one empty monitor. Producing no track" << std::endl;
    return tracks;
  }

  //Convention: (Horiz, Vert)
  std::vector< std::pair<short, short> > UpstreamPairedFibers;
  std::vector< std::pair<short, short> > DownstreamPairedFibers;

  std::vector< TVector3 > UpstreamPositions;
  std::vector< TVector3 > DownstreamPositions;

  for(size_t iH = 0; iH < HorizUpstreamFibers.size(); ++iH){

    size_t HorizFiber = HorizUpstreamFibers[iH];

    for(size_t iV = 0; iV < VertUpstreamFibers.size(); ++iV){
      size_t VertFiber = VertUpstreamFibers[iV];

      //MF_LOG_DEBUG("BeamEvent") << "Paired: " << HorizFiber << " " << VertFiber << "\n"; 
      std::cout << "Paired: " << HorizFiber << " " << VertFiber << std::endl; 
      UpstreamPairedFibers.push_back(std::make_pair(HorizFiber, VertFiber));

      //If there's 2 adjacent fibers. Skip the next. Could replace this with averaging the position
      //todo later
      if (iV < VertUpstreamFibers.size() - 1){
        if (VertUpstreamFibers[iV] == (VertUpstreamFibers[iV + 1] - 1)) ++iV;
      }    
    }    

    if (iH < HorizUpstreamFibers.size() - 1){
      if (HorizUpstreamFibers[iH] == (HorizUpstreamFibers[iH + 1] - 1)) ++iH;
    }    
  } 
  
  std::cout << "Upstream " << std::endl;

  for(size_t i = 0; i < UpstreamPairedFibers.size(); ++i){
    
    std::pair<short,short> thePair = UpstreamPairedFibers.at(i);
 
    double xPos = GetPosition(thePair.first);
    double yPos = GetPosition(thePair.second);
    
    std::cout << "normal " << xPos << " " << yPos << std::endl;
    TVector3 posInDet = ConvertMonitorCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
    std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
    UpstreamPositions.push_back( posInDet );
  }

  for(size_t iH = 0; iH < HorizDownstreamFibers.size(); ++iH){

    size_t HorizFiber = HorizDownstreamFibers[iH];

    for(size_t iV = 0; iV < VertDownstreamFibers.size(); ++iV){
      size_t VertFiber = VertDownstreamFibers[iV];

      //MF_LOG_DEBUG("BeamEvent") << "Paired: " << HorizFiber << " " << VertFiber << "\n"; 
      std::cout << "Paired: " << HorizFiber << " " << VertFiber << std::endl; 
      DownstreamPairedFibers.push_back(std::make_pair(HorizFiber, VertFiber));

      //If there's 2 adjacent fibers. Skip the next. Could replace this with averaging the position
      //todo later
      if (iV < VertDownstreamFibers.size() - 1){
        if (VertDownstreamFibers[iV] == (VertDownstreamFibers[iV + 1] - 1)) ++iV;
      }    
    }    

    if (iH < HorizDownstreamFibers.size() - 1){
      if (HorizDownstreamFibers[iH] == (HorizDownstreamFibers[iH + 1] - 1)) ++iH;
    }    
  } 
  
  std::cout << "Downstream " << std::endl;

  for(size_t i = 0; i < DownstreamPairedFibers.size(); ++i){
    
    std::pair<short,short> thePair = DownstreamPairedFibers.at(i);
 
    double xPos = GetPosition(thePair.first);
    double yPos = GetPosition(thePair.second);
    
    std::cout << "normal " << xPos << " " << yPos << std::endl;
    TVector3 posInDet = ConvertMonitorCoordinates(xPos,yPos,0.,fSecondTrackingProfZ);
    std::cout << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << std::endl;
    DownstreamPositions.push_back( posInDet );
  }

  for(size_t iU = 0; iU < UpstreamPositions.size(); ++iU){
    for(size_t iD = 0; iD < DownstreamPositions.size(); ++iD){
      std::vector<TVector3> thePoints;
      thePoints.push_back(UpstreamPositions.at(iU));
      thePoints.push_back(DownstreamPositions.at(iD));

      //Now project the last point to the TPC face
      thePoints.push_back( ProjectToTPC(thePoints[0],thePoints[1]) );    
      std::cout << "Projected: " << thePoints.back().X() << " " << thePoints.back().Y() << " " << thePoints.back().Z() << std::endl;

     
      std::vector<TVector3> theMomenta;
      //Just push back the unit vector for each point 
      theMomenta.push_back( ( DownstreamPositions.at(iD) - UpstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( DownstreamPositions.at(iD) - UpstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( DownstreamPositions.at(iD) - UpstreamPositions.at(iU) ).Unit() );

      recob::Track tempTrack(recob::TrackTrajectory(recob::tracking::convertCollToPoint(thePoints),
						    recob::tracking::convertCollToVector(theMomenta),
						    recob::Track::Flags_t(thePoints.size()), false),
			     0, -1., 0, recob::tracking::SMatrixSym55(), recob::tracking::SMatrixSym55(), 1);
      tracks.push_back( tempTrack );
    }    
  }
  
  return tracks;
}

void protoana::ProtoDUNEBeamlineUtils::reconfigure(fhicl::ParameterSet const& p){
  fBeamEventTag = p.get<art::InputTag>("BeamEventTag");
  fUseCERNCalibSelection = p.get<bool>("UseCERNCalibSelection");

  //Tracking parameters
  fBeamX = p.get<double>("BeamX");
  fBeamY = p.get<double>("BeamY");
  fBeamZ = p.get<double>("BeamZ");

  fRotateMonitorXZ = p.get<double>("RotateMonitorXZ"); 
  fRotateMonitorYZ = p.get<double>("RotateMonitorYZ"); 
  fRotateMonitorYX = p.get<double>("RotateMonitorYX"); 

  fFirstTrackingProfZ  = p.get<double>("FirstTrackingProfZ");
  fSecondTrackingProfZ = p.get<double>("SecondTrackingProfZ");
  fNP04FrontZ          = p.get<double>("NP04FrontZ");
  ////////////////////////


  //Momentum parameters
  fBeamBend = p.get<double>("BeamBend");
  L1 = p.get<double>("L1");
  L2 = p.get<double>("L2");
  L3 = p.get<double>("L3");

  fBProf1Shift = p.get< double >("BProf1Shift"); 
  fBProf2Shift = p.get< double >("BProf2Shift"); 
  fBProf3Shift = p.get< double >("BProf3Shift"); 
  ///////////////////////
  
  // Justin's stuff from BeamlineUtils
  fMomentumScaleFactor = p.get<float>("MomentumScaleFactor"); 
  fMomentumOffset = p.get<float>("MomentumOffset"); // GeV/c
  fTOFScaleFactor = p.get<float>("TOFScaleFactor");
  fTOFOffset = p.get<float>("TOFOffset"); // ns
}

double protoana::ProtoDUNEBeamlineUtils::GetPosition( short theFiber ){
  //Fibers are 1mm in width
  return 1.*(96 - theFiber) - .5;
}

TVector3 protoana::ProtoDUNEBeamlineUtils::ConvertMonitorCoordinates(double x, double y, double z, double zOffset){

  if( !rotated ) BeamMonitorBasisVectors();

  double off = fNP04FrontZ - zOffset;

  TVector3 old(x,y,z);

  double newX = x*MonitorBasisX.X() + y*MonitorBasisY.X() + off*fabs(MonitorBasisZ.X());
  double newY = x*MonitorBasisX.Y() + y*MonitorBasisY.Y() + off*fabs(MonitorBasisZ.Y());
  double newZ = x*MonitorBasisX.Z() + y*MonitorBasisY.Z() - off*fabs(MonitorBasisZ.Z());

  newX += fBeamX*10.;
  newY += fBeamY*10.;
  newZ += fBeamZ*10.;

  TVector3 result(newX/10., newY/10., newZ/10.);
  return result;
}

void protoana::ProtoDUNEBeamlineUtils::BeamMonitorBasisVectors(){
  MonitorBasisX = TVector3( 1., 0., 0. );
  MonitorBasisY = TVector3( 0., 1., 0. );
  MonitorBasisZ = TVector3( 0., 0., 1. );

  RotateMonitorVector(MonitorBasisX);
  RotateMonitorVector(MonitorBasisY);
  RotateMonitorVector(MonitorBasisZ);
  

  rotated = true;
}

void protoana::ProtoDUNEBeamlineUtils::RotateMonitorVector(TVector3 &vec){
  vec.RotateX( fRotateMonitorYZ * TMath::Pi()/180. );
  vec.RotateY( fRotateMonitorXZ * TMath::Pi()/180. );
}


TVector3 protoana::ProtoDUNEBeamlineUtils::ProjectToTPC(TVector3 firstPoint, TVector3 secondPoint){
  TVector3 dR = (secondPoint - firstPoint);

  double deltaZ = -1.*secondPoint.Z();
  double deltaX = deltaZ * (dR.X() / dR.Z());
  double deltaY = deltaZ * (dR.Y() / dR.Z());

  TVector3 lastPoint = secondPoint + TVector3(deltaX, deltaY, deltaZ);
  return lastPoint;
}

std::vector< double > protoana::ProtoDUNEBeamlineUtils::MomentumSpec( art::Event const & evt ){
 
  GetCurrent( evt );
 
  std::cout << "Got current: " << Current << std::endl;

  std::vector< double > theMomenta;

  if( Current < .000000001 ){
    std::cout << "Warning! Low magnet current. Momentum spectrometry invalid. Returning empty momenta vector." << std::endl;
    return theMomenta;
  }



  double LB = mag_P1*fabs(Current);
  double deltaI = fabs(Current) - mag_P4;

  if(deltaI>0) LB+= mag_P3*deltaI*deltaI;

  //Get the active fibers from the upstream tracking XBPF
  std::vector<short> BProf1Fibers = ActiveFibers.at(BProf1); 
  std::vector<short> BProf2Fibers = ActiveFibers.at(BProf2); 
  std::vector<short> BProf3Fibers = ActiveFibers.at(BProf3);


  std::cout << BProf1 << " has " << BProf1Fibers.size() << " active fibers" << std::endl;
  for(size_t i = 0; i < BProf1Fibers.size(); ++i){
    std::cout << BProf1Fibers[i] << " ";
  }
  std::cout << std::endl;

  std::cout << BProf2 << " has " << BProf2Fibers.size() << " active fibers" << std::endl;
  for(size_t i = 0; i < BProf2Fibers.size(); ++i){
    std::cout << BProf2Fibers[i] << " ";
  }
  std::cout << std::endl;

  std::cout << BProf3 << " has " << BProf3Fibers.size() << " active fibers" << std::endl;
  for(size_t i = 0; i < BProf3Fibers.size(); ++i){
    std::cout << BProf3Fibers[i] << " ";
  }
  std::cout << std::endl;


  if( (BProf1Fibers.size() < 1) || (BProf2Fibers.size() < 1) || (BProf3Fibers.size() < 1) ){
    std::cout << "Warning, at least one empty Beam Profiler. Not checking momentum" << std::endl;
    return theMomenta;
  }


  std::cout << "Getting all trio-wise hits" << std::endl;
  std::cout << "N1,N2,N3 " << BProf1Fibers.size()
            << " "         << BProf2Fibers.size() 
            << " "         << BProf3Fibers.size() << std::endl;

  for(size_t i1 = 0; i1 < BProf1Fibers.size(); ++i1){

    double x1,x2,x3;

    x1 = -1.*GetPosition(BProf1Fibers[i1])/1.E3;

    if (i1 < BProf1Fibers.size() - 1){
      if (BProf1Fibers[i1] == (BProf1Fibers[i1 + 1] - 1)){
        //Add .5 mm
        x1 += .0005;
      }
    }

    for(size_t i2 = 0; i2 < BProf2Fibers.size(); ++i2){

      x2 = -1.*GetPosition(BProf2Fibers[i2])/1.E3;

      if (i2 < BProf2Fibers.size() - 1){
        if (BProf2Fibers[i2] == (BProf2Fibers[i2 + 1] - 1)){
          //Add .5 mm
          x2 += .0005;
        }
      }

      for(size_t i3 = 0; i3 < BProf3Fibers.size(); ++i3){

        std::cout << "\t" << i1 << " " << i2 << " " << i3 << std::endl;

        x3 = -1.*GetPosition(BProf3Fibers[i3])/1.E3;

        if (i3 < BProf3Fibers.size() - 1){
          if (BProf3Fibers[i3] == (BProf3Fibers[i3 + 1] - 1)){
            //Add .5 mm
            x3 += .0005;
          }
        }

        //Calibrate the positions of the fibers in the
        //monitors
        x1 -= fBProf1Shift*1.e-3;
        x2 -= fBProf2Shift*1.e-3;
        x3 -= fBProf3Shift*1.e-3;

        double cosTheta = MomentumCosTheta(x1,x2,x3);        
        double momentum = 299792458*LB/(1.E9 * acos(cosTheta));

        theMomenta.push_back(momentum);

        if (i3 < BProf3Fibers.size() - 1){
          if (BProf3Fibers[i3] == (BProf3Fibers[i3 + 1] - 1)){
            //Skip the next
            ++i3;
          }
        }
        
      }

      if (i2 < BProf2Fibers.size() - 1){
        if (BProf2Fibers[i2] == (BProf2Fibers[i2 + 1] - 1)){
          //Skip the next
          ++i2;
        }
      }
    }

    if (i1 < BProf1Fibers.size() - 1){
      if (BProf1Fibers[i1] == (BProf1Fibers[i1 + 1] - 1)){
        //Skip the next
        ++i1;
      }
    }
  }

  return theMomenta;
}

double protoana::ProtoDUNEBeamlineUtils::MomentumCosTheta( double X1, double X2, double X3 ){
  double a =  (X2*L3 - X3*L2)*cos(fBeamBend)/(L3-L2);

 
  double numTerm = (a - X1)*( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) ) + L1*( L3 - L2 );

  double denomTerm1, denomTerm2, denom;
  denomTerm1 = sqrt( L1*L1 + (a - X1)*(a - X1) );
  denomTerm2 = sqrt( TMath::Power( ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) ),2)
                   + TMath::Power( ( (L3 - L2) ),2) );
  denom = denomTerm1 * denomTerm2;

  double cosTheta = numTerm/denom;  
  return cosTheta;
}

protoana::PossibleParticleCands protoana::ProtoDUNEBeamlineUtils::GetPIDCandidates( beam::ProtoDUNEBeamEvent const & beamevt, double nominal_momentum ){
  return GetPIDCandidates_CERNCalib(beamevt,nominal_momentum);
}

protoana::PossibleParticleCands protoana::ProtoDUNEBeamlineUtils::GetPIDCandidates( art::Event const & evt, double nominal_momentum ){
  auto beamHand = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  for(size_t iBeamEvent=0; iBeamEvent < beamHand->size(); iBeamEvent++)
  {
    return GetPIDCandidates(beamHand->at(iBeamEvent),nominal_momentum);
  }
  return protoana::PossibleParticleCands();
}


std::vector< int > protoana::ProtoDUNEBeamlineUtils::GetPID( beam::ProtoDUNEBeamEvent const & beamevt, double nominal_momentum ){
  const auto& thePIDCands = GetPIDCandidates(beamevt, nominal_momentum);
  std::vector< int > thePIDs = thePIDCands.getPDGCodes();
  return thePIDs;
}

std::vector<int> protoana::ProtoDUNEBeamlineUtils::GetPID( art::Event const & evt, double nominal_momentum ){
  auto beamHand = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  for(size_t iBeamEvent=0; iBeamEvent < beamHand->size(); iBeamEvent++)
  {
    return GetPID(beamHand->at(iBeamEvent),nominal_momentum);
  }
  return std::vector<int>();
}

protoana::PossibleParticleCands protoana::ProtoDUNEBeamlineUtils::GetPIDCandidates_CERNCalib( beam::ProtoDUNEBeamEvent const & beamevt, double nominal_momentum ){
 
  PossibleParticleCands candidates;
  
  //Check if momentum is in valid set
  std::vector< double > valid_momenta = {1., 2., 3., 6., 7.};
  if( std::find(valid_momenta.begin(), valid_momenta.end(), nominal_momentum) == valid_momenta.end() ){
    std::cout << "Reference momentum " << nominal_momentum << " not valid" << std::endl;
    return candidates;
  }

  //Get the high/low pressure Cerenkov info
  //int high_pressure_status, low_pressure_status; 
  int high_pressure_status = beamevt.GetCKov0Status();
  int low_pressure_status = beamevt.GetCKov1Status();
  
  //std::cout << "Pressures: " << beamevt.GetCKov0Pressure() << " " << beamevt.GetCKov1Pressure() << std::endl;
  //if( beamevt.GetCKov0Pressure() < beamevt.GetCKov1Pressure() ){
  //  high_pressure_status = beamevt.GetCKov1Status();
  //  low_pressure_status = beamevt.GetCKov0Status();
  //}
  //else{
  //  high_pressure_status = beamevt.GetCKov0Status();
  //  low_pressure_status = beamevt.GetCKov1Status();
  //}
  

  

  if( nominal_momentum == 1. ){
    if( beamevt.GetTOFChan() == -1 ){
      std::cout << "TOF invalid" << std::endl;
      return candidates;
    }
    if( low_pressure_status == -1 ){
      std::cout << "High pressure status invalid" << std::endl;
      return candidates;
    }

    const double & tof = beamevt.GetTOF();
    if ( 
        ((fUseCERNCalibSelection && tof < 105.) 
            || (!fUseCERNCalibSelection && tof < 170.))
        && low_pressure_status == 1 
       ) {
      candidates.electron = true;
    }
    else if ( 
        ((fUseCERNCalibSelection && tof < 110.) 
            || (!fUseCERNCalibSelection && tof < 170.))
        && low_pressure_status == 0 ){
      candidates.muon = true;
      candidates.pion = true;
    }
    else if ( 
        ((fUseCERNCalibSelection && tof > 110. && tof < 160.) 
            || (!fUseCERNCalibSelection && tof > 170.))
        && low_pressure_status == 0 ) {
      candidates.proton = true;
    }
  }
  else if( nominal_momentum == 2. ){
    if( beamevt.GetTOFChan() == -1 ){
      std::cout << "TOF invalid" << std::endl;
      return candidates;
    }
    if( low_pressure_status == -1 ){
      std::cout << "High pressure Cerenkov status invalid" << std::endl;
      return candidates;
    }

    const double & tof = beamevt.GetTOF();
    if ( 
        ((fUseCERNCalibSelection && tof < 105.) 
            || (!fUseCERNCalibSelection && tof < 160.))
        && low_pressure_status == 1 
       ) {
      candidates.electron = true;
    }
    else if ( 
        ((fUseCERNCalibSelection && tof < 103.) 
            || (!fUseCERNCalibSelection && tof < 160.))
        && low_pressure_status == 0 ){
      candidates.muon = true;
      candidates.pion = true;
    }
    else if ( 
        ((fUseCERNCalibSelection && tof > 103. && tof < 160.) 
            || (!fUseCERNCalibSelection && tof > 160.))
        && low_pressure_status == 0 ) {
      candidates.proton = true;
    }
  }
  else if( nominal_momentum == 3. ){
    if( high_pressure_status == -1 || low_pressure_status == -1 ){
      std::cout << "At least one Cerenkov status invalid " << std::endl;
      std::cout << "High: " << high_pressure_status << " Low: " << low_pressure_status << std::endl;
      return candidates;
    }
    else if ( low_pressure_status == 1 && high_pressure_status == 1 ) 
      candidates.electron = true;

    else if ( low_pressure_status == 0 && high_pressure_status == 1 ){
      candidates.muon = true;
      candidates.pion = true;
    }

    else{ // low, high = 0, 0
      candidates.proton = true;
      candidates.kaon = true; 
    }
  }
  else if( nominal_momentum == 6. || nominal_momentum == 7. ){
    if( high_pressure_status == -1 || low_pressure_status == -1 ){
      std::cout << "At least one Cerenkov status invalid " << std::endl;
      std::cout << "High: " << high_pressure_status << " Low: " << low_pressure_status << std::endl;
      return candidates;
    }
    else if ( low_pressure_status == 1 && high_pressure_status == 1 ){
      candidates.electron = true;
      candidates.muon = true;
      candidates.pion = true;
    }
    else if ( low_pressure_status == 0 && high_pressure_status == 1 ) 
      candidates.kaon = true; 

    else  // low, high = 0, 0
      candidates.proton = true;
  }

  return candidates;
 
}

std::vector< int > protoana::ProtoDUNEBeamlineUtils::GetPID_CERNCalib( beam::ProtoDUNEBeamEvent const & beamevt, double nominal_momentum ){
  const auto& thePIDCands = GetPIDCandidates_CERNCalib(beamevt, nominal_momentum);
  std::vector< int > thePIDs = thePIDCands.getPDGCodes();
  return thePIDs;
}

double protoana::ProtoDUNEBeamlineUtils::ComputeTOF( int pdg, double momentum ){
  if( momentum  < .0000001 ){ 
    std::cout << "Low Momentum" << std::endl;
    return -1.;
  }

  if( particle_mass.find( pdg ) == particle_mass.end() ){
    std::cout << "PDG " << pdg << " not found" <<  std::endl;
    return -1.;
  }

  double m = particle_mass[pdg];
  double tof = fTOFDist / c;
  tof = tof / sqrt( 1. - m*m / (momentum*momentum + m*m) );

  return tof * 1.e9;
}

double protoana::ProtoDUNEBeamlineUtils::ComputeMomentum( int pdg, double tof ){
  if( tof < .0000001 ){ 
    std::cout << "Low TOF" << std::endl;
    return -1.;
  }

  if( particle_mass.find( pdg ) == particle_mass.end() ){
    std::cout << "PDG " << pdg << " not found" <<  std::endl;
    return -1.;
  }

  double m = particle_mass[pdg];
  return m * sqrt(1. / ( 1. - fTOFDist*fTOFDist / ( 1.e-18*tof*tof*c*c ) )  - 1. );

}

// ----------------------------------------------------------------------------
bool protoana::ProtoDUNEBeamlineUtils::IsGoodBeamlineTrigger(art::Event const & evt) const{
  return IsGoodBeamlineTrigger( GetBeamEvent(evt) );
}

bool protoana::ProtoDUNEBeamlineUtils::IsGoodBeamlineTrigger(beam::ProtoDUNEBeamEvent const & beamEvent) const{
  return (beamEvent.GetTimingTrigger() == 12 && beamEvent.CheckIsMatched());
}


bool protoana::ProtoDUNEBeamlineUtils::HasPerfectBeamMomentum(art::Event const & evt) const{
  return HasPerfectBeamMomentum( GetBeamEvent(evt) );
}

bool protoana::ProtoDUNEBeamlineUtils::HasPerfectBeamMomentum(beam::ProtoDUNEBeamEvent const & beamEvent) const{

  //Get Active fibers
  std::vector< short > fibers_p1 = beamEvent.GetActiveFibers( "XBPF022697" );
  std::vector< short > fibers_p2 = beamEvent.GetActiveFibers( "XBPF022701" );
  std::vector< short > fibers_p3 = beamEvent.GetActiveFibers( "XBPF022702" );

  //Get any glitches and remove from the active
  std::array< short, 192 > glitch_mask = beamEvent.GetFBM( "XBPF022697" ).glitch_mask;
  int nP1 = 0;
  for( size_t i = 0; i < fibers_p1.size(); ++i ){
    if( !glitch_mask[ fibers_p1[i] ] )
      ++nP1;   
  }

  glitch_mask = beamEvent.GetFBM( "XBPF022701" ).glitch_mask;
  int nP2 = 0;
  for( size_t i = 0; i < fibers_p2.size(); ++i ){
    if( !glitch_mask[ fibers_p2[i] ] )
      ++nP2;   
  }

  glitch_mask = beamEvent.GetFBM( "XBPF022702" ).glitch_mask;
  int nP3 = 0;
  for( size_t i = 0; i < fibers_p3.size(); ++i ){
    if( !glitch_mask[ fibers_p3[i] ] )
      ++nP3;   
  }

  return ( nP1 == 1 && nP2 == 1 && nP3 == 1 );
}


std::vector<double> protoana::ProtoDUNEBeamlineUtils::GetBeamlineMass(art::Event const & evt) const{
  std::vector<double> result;
  
  const auto & massSquareds = GetBeamlineMassSquared(evt);
  for (const auto & massSquared: massSquareds)
  {
    const auto & mass = std::sqrt(massSquared);
    result.push_back(mass);
  }
  return result;
}

std::vector<double> protoana::ProtoDUNEBeamlineUtils::GetBeamlineMassSquared(art::Event const & evt) const{
  std::vector<double> tofs;
  std::vector<double> momenta;
  std::vector<double> result;

  auto const& beamEvents = evt.getProduct<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  for(auto const& beamEvent : beamEvents)
  {
    tofs.push_back(beamEvent.GetTOF());

    const std::vector<double> & beamMomenta = beamEvent.GetRecoBeamMomenta();
    for(size_t iMom=0; iMom < beamMomenta.size(); iMom++)
    {
      momenta.push_back(beamEvent.GetRecoBeamMomentum(iMom));
    }
  }
  for(const auto & tof : tofs)
  {
    for(const auto & momentum : momenta)
    {
        const double massSquared = std::pow(momentum*fMomentumScaleFactor-fMomentumOffset,2) 
                    * ((tof*fTOFScaleFactor-fTOFOffset)/(fTOFDist/c*1e9) - 1.);
        result.push_back(massSquared);
    }
  }
  return result;
}

const std::tuple<double,double,int,int> protoana::ProtoDUNEBeamlineUtils::GetBeamlineVars(art::Event const & evt) const {

  double momentum = -99999.;
  double tof = -99999.;
  int ckov0 = -99999.;
  int ckov1 = -99999.;

  auto const& beamEvents = evt.getProduct<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  for(auto const& beamEvent : beamEvents)
  {
    tof = beamEvent.GetTOF();
    ckov0 = beamEvent.GetCKov0Status();
    ckov1 = beamEvent.GetCKov1Status();

    const std::vector<double> & beamMomenta = beamEvent.GetRecoBeamMomenta();
    if (beamMomenta.size() > 0)
    {
        momentum = beamEvent.GetRecoBeamMomentum(0);
    }
  }
  return std::make_tuple(momentum, tof, ckov0, ckov1);
}

const std::tuple<double,double,int,int,int,double,double,int,int,bool> protoana::ProtoDUNEBeamlineUtils::GetBeamlineVarsAndStatus(art::Event const & evt) const {

  double momentum = -99999.;
  double tof = -99999.;
  int tofChannel = -99999.;
  int ckov0 = -99999.;
  int ckov1 = -99999.;
  double ckov0Pressure = -99999.;
  double ckov1Pressure = -99999.;
  int timingTrigger = -99999.;
  int BITrigger = -99999.;
  bool areBIAndTimingMatched = false;

  auto const& beamEvents = evt.getProduct<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
  for(auto const& beamEvent : beamEvents) {
    tof = beamEvent.GetTOF();
    tofChannel = beamEvent.GetTOFChan();
    ckov0 = beamEvent.GetCKov0Status();
    ckov1 = beamEvent.GetCKov1Status();
    ckov0Pressure = beamEvent.GetCKov0Pressure();
    ckov1Pressure = beamEvent.GetCKov1Pressure();
    timingTrigger = beamEvent.GetTimingTrigger();
    BITrigger = beamEvent.GetBITrigger();
    areBIAndTimingMatched = beamEvent.CheckIsMatched();

    const std::vector<double> & beamMomenta = beamEvent.GetRecoBeamMomenta();
    if (beamMomenta.empty()) {
    break;
  }

    momentum = beamEvent.GetRecoBeamMomentum(0);
}

  return std::make_tuple(momentum, tof, tofChannel, ckov0, ckov1, ckov0Pressure, ckov1Pressure, timingTrigger, BITrigger, areBIAndTimingMatched);
}
