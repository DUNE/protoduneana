////////////////////////////////////////////////////////////////////////
// Class:       PDVDbeamSelection
// Plugin Type: analyzer (Unknown Unknown)
// File:        PDVDbeamSelection_module.cc
//
// Generated at Fri Jan 23 11:56:46 2026 by Yoann Kermaidic using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/ArtDataHelper/MVAReader.h"

#include "dunecore/DuneObj/ProtoDUNEBeamSpill.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

#include "TTree.h"

namespace pdvd {
  class PDVDbeamSelection;
}


class pdvd::PDVDbeamSelection : public art::EDAnalyzer {
public:
  explicit PDVDbeamSelection(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDVDbeamSelection(PDVDbeamSelection const&) = delete;
  PDVDbeamSelection(PDVDbeamSelection&&) = delete;
  PDVDbeamSelection& operator=(PDVDbeamSelection const&) = delete;
  PDVDbeamSelection& operator=(PDVDbeamSelection&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  double MomentumCosTheta(double,double,double);
  float  GetRecoBeamMomentum(art::Ptr<beam::ProtoDUNEBeamEvent>);
  double GetPosition(std::string, int);

private:

  protoana::ProtoDUNETrackUtils trackUtil;
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  int         fLogLevel;
  std::string fBeamEventLabel;
  std::string fHitLabel;
  std::string fSliceLabel;
  std::string fPFParticleLabel;
  std::string fSpacePointLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPIDLabel;
  bool        fUseMVA;
  bool        fExcludeMuons;

  //////////
  // Needed for beam momentum reconstructed
  // based on BeamEvent_module

  double L1;
  double L2;
  double L3;

  //Hardware Parameters for magnetic field stuff
  double mag_P1;
  double mag_P3;
  double mag_P4;

  std::vector< std::string > fDevices;

  float fBProf1Shift;
  float fBProf2Shift;
  float fBProf3Shift;

  std::string firstBPROF1;
  std::string secondBPROF1;
  std::string BPROF2;
  std::string BPROF3;

  float fBeamBend;

  double fFiberDimension;

  std::map<std::string, std::string > fDeviceTypes;


  ///////////

  TTree* fTree;

  unsigned int fEventID;
  float fBeamTof;
  short fBeamCKov0;
  short fBeamCKov1;
  double fBeamMomentum;
  float fMagnetCurrent;

  int fNPFParticles;
  int fNPrimaryChildren;
  int fPrimaryPdg;
  int fPrimaryKey;

  std::vector<int>  fChildBeamPFPKey;
  std::vector<int>  fChildBeamPFPPdg;
  std::vector<int>  fChildBeamPFPChildren;

  float fPrimBeamCNNtrack;
  float fPrimBeamCNNshower;
  float fPrimBeamCNNmichel;

  float fPrimBeamPFPTrackScore;

  int   fPrimBeamTrackPid;
  float fPrimBeamTrackVertexX;
  float fPrimBeamTrackVertexY;
  float fPrimBeamTrackVertexZ;
  float fPrimBeamTrackDirX;
  float fPrimBeamTrackDirY;
  float fPrimBeamTrackDirZ;
  float fPrimBeamTrackLength;

  float fPrimBeamShowerVertexX;
  float fPrimBeamShowerVertexY;
  float fPrimBeamShowerVertexZ;
  float fPrimBeamShowerDirX;
  float fPrimBeamShowerDirY;
  float fPrimBeamShowerDirZ;
  float fPrimBeamShowerLength;
  float fPrimBeamShowerEnergy;
  float fPrimBeamShowerdEdx;
  float fPrimBeamShowerOpeningAngle;

  std::vector<float> fPrimBeamX;
  std::vector<float> fPrimBeamY;
  std::vector<float> fPrimBeamZ;

  std::vector<int>   fChildBeamTrackPid;

  std::vector<float> fChildBeamCNNtrack;
  std::vector<float> fChildBeamCNNshower;
  std::vector<float> fChildBeamCNNmichel;

  std::vector<float> fChildBeamPFPTrackScore;
  std::vector<float> fChildBeamTrackVertexX;
  std::vector<float> fChildBeamTrackVertexY;
  std::vector<float> fChildBeamTrackVertexZ;
  std::vector<float> fChildBeamTrackDirX;
  std::vector<float> fChildBeamTrackDirY;
  std::vector<float> fChildBeamTrackDirZ;
  std::vector<float> fChildBeamTrackLength;

  std::vector<float> fChildBeamShowerVertexX;
  std::vector<float> fChildBeamShowerVertexY;
  std::vector<float> fChildBeamShowerVertexZ;
  std::vector<float> fChildBeamShowerDirX;
  std::vector<float> fChildBeamShowerDirY;
  std::vector<float> fChildBeamShowerDirZ;
  std::vector<float> fChildBeamShowerLength;
  std::vector<float> fChildBeamShowerEnergy;
  std::vector<float> fChildBeamShowerdEdx;
  std::vector<float> fChildBeamShowerOpeningAngle;

  std::vector<std::vector<float>> fChildBeamX;
  std::vector<std::vector<float>> fChildBeamY;
  std::vector<std::vector<float>> fChildBeamZ;

};


pdvd::PDVDbeamSelection::PDVDbeamSelection(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fLogLevel(       p.get<int>("LogLevel")),
  fBeamEventLabel( p.get<std::string>("BeamEventLabel")), 
  fHitLabel(       p.get<std::string>("HitLabel")),
  fSliceLabel(     p.get<std::string>("SliceLabel")), 
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")), 
  fSpacePointLabel(p.get<std::string>("SpacePointLabel")), 
  fTrackLabel(     p.get<std::string>("TrackLabel")), 
  fShowerLabel(    p.get<std::string>("ShowerLabel")),
  fPIDLabel(       p.get<std::string>("PIDLabel")),
  fUseMVA(         p.get<bool>("UseMVA")),
  fExcludeMuons(   p.get<bool>("ExcludeMuons")),
  L1(p.get<double>("L1")), 
  L2(p.get<double>("L2")), 
  L3(p.get<double>("L3")),
  mag_P1(p.get<double>("mag_P1")), 
  mag_P3(p.get<double>("mag_P3")), 
  mag_P4(p.get<double>("mag_P4")),
  fDevices(p.get<std::vector<std::string>>("Devices")), 
  fBProf1Shift(p.get<double>("BProf1Shift")), 
  fBProf2Shift(p.get<double>("BProf2Shift")), 
  fBProf3Shift(p.get<double>("BProf3Shift")),
  firstBPROF1(p.get<std::string>("FirstBPROF1")),
  secondBPROF1(p.get<std::string>("SecondBPROF1")),
  BPROF2(p.get<std::string>("BPROF2")),
  BPROF3(p.get<std::string>("BPROF3")),
  fBeamBend(p.get<double>("BeamBend")),
  fFiberDimension(p.get<double>("FiberDimension", 1.))
{

  std::vector< std::pair<std::string, std::string> >  tempTypes = p.get<std::vector< std::pair<std::string, std::string> >>("DeviceTypes");
  fDeviceTypes  = std::map<std::string, std::string>(tempTypes.begin(), tempTypes.end() );

}

double pdvd::PDVDbeamSelection::MomentumCosTheta(double X1, double X2, double X3){
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

double pdvd::PDVDbeamSelection::GetPosition(std::string deviceName, int fiberIdx){
  if(fiberIdx > 192){ return -1.;}
  double size = fFiberDimension; //[deviceName];
  
  //Define 0th fiber as farthest positive. Last fiber is farthest negative. Center is between 96 and 97 
  double pos = size*(96 - fiberIdx) - size/2.;
  return pos;
}

float pdvd::PDVDbeamSelection::GetRecoBeamMomentum(art::Ptr<beam::ProtoDUNEBeamEvent> beaminfo){

  double momentum_full = 0.0;

  double LB = mag_P1*fabs(fMagnetCurrent);
  double deltaI = fabs(fMagnetCurrent) - mag_P4;
  if(deltaI>0) LB+= mag_P3*deltaI*deltaI;

  //Get the active fibers from the upstream tracking XBPF
  std::string firstBPROF1Type    = fDeviceTypes[firstBPROF1]; 
  std::string secondBPROF1Type   = fDeviceTypes[secondBPROF1]; 
  std::vector<short> BPROF1Fibers;
  std::string BPROF1Name;

  if (firstBPROF1Type == "horiz" && secondBPROF1Type == "vert"){
    for(size_t iF = 0; iF < 192; ++iF){
      if(beaminfo->GetFBM(firstBPROF1).fibers[iF] && !beaminfo->GetFBM(firstBPROF1).glitch_mask[iF]) BPROF1Fibers.push_back(iF); 
    }
    BPROF1Name = firstBPROF1;
  }
  else if(secondBPROF1Type == "horiz" && firstBPROF1Type == "vert"){
    for(size_t iF = 0; iF < 192; ++iF){
	if(beaminfo->GetFBM(secondBPROF1).fibers[iF] && !beaminfo->GetFBM(secondBPROF1).glitch_mask[iF]) BPROF1Fibers.push_back(iF);
    }
    BPROF1Name = secondBPROF1;
  }
  else{
    return 0.0;
  }

  //////////////////////////////////////////////

  if( (BPROF1Fibers.size() < 1) ){
    return 0.0;
  }
  //We have the active Fibers, now go through them.
  //Skip the second of any adjacents 
  
  //BPROF2////
  //
  
  std::vector<short> BPROF2Fibers;
  for(size_t iF = 0; iF < 192; ++iF){
    if(beaminfo->GetFBM(BPROF2).fibers[iF] && !beaminfo->GetFBM(BPROF2).glitch_mask[iF]) BPROF2Fibers.push_back(iF);
  }

  if( (BPROF2Fibers.size() < 1) ){
    return 0.0;
  }
  
  ////////////

  //BPROF3////
  //

  std::vector<short> BPROF3Fibers;
  for(size_t iF = 0; iF < 192; ++iF){
    if(beaminfo->GetFBM(BPROF3).fibers[iF] && !beaminfo->GetFBM(BPROF3).glitch_mask[iF]) BPROF3Fibers.push_back(iF);
  }

  if( (BPROF3Fibers.size() < 1) ){
    return 0.0;
  }

  for(size_t i1 = 0; i1 < BPROF1Fibers.size(); ++i1){
    double x1,x2,x3;

    x1 = -1.*GetPosition(BPROF1Name, BPROF1Fibers[i1])/1.E3;
    if (i1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[i1] == (BPROF1Fibers[i1 + 1] - 1)){
        //Add .5 mm
        x1 += .0005;
      }
    }

    for(size_t i2 = 0; i2 < BPROF2Fibers.size(); ++i2){
      x2 = -1.*GetPosition(BPROF2, BPROF2Fibers[i2])/1.E3;
      if (i2 < BPROF2Fibers.size() - 1){
        if (BPROF2Fibers[i2] == (BPROF2Fibers[i2 + 1] - 1)){
          //Add .5 mm
          x2 += .0005;
        }
      }

      for(size_t i3 = 0; i3 < BPROF3Fibers.size(); ++i3){
        
        x3 = -1.*GetPosition(BPROF3, BPROF3Fibers[i3])/1.E3;
        if (i3 < BPROF3Fibers.size() - 1){
          if (BPROF3Fibers[i3] == (BPROF3Fibers[i3 + 1] - 1)){
            //Add .5 mm
            x3 += .0005;
          }
        }

        //Calibrate the positions
        //-1.*( FiberPos ) -> -1.*( FiberPos + ShiftDist )
        // = -1.*FiberPos - ShiftDist
        x1 = x1 - fBProf1Shift*1.e-3; 
        x2 = x2 - fBProf2Shift*1.e-3; 
        x3 = x3 - fBProf3Shift*1.e-3; 

        double cosTheta_full = MomentumCosTheta(x1,x2,x3);        
        momentum_full = 299792458*LB/(1.E9 * acos(cosTheta_full));

        return momentum_full;
      }
    }
  }
  return 0.0;
}

void pdvd::PDVDbeamSelection::analyze(art::Event const& e)
{
  fEventID = e.id().event();
  fBeamTof = 0.0;
  fBeamMomentum = 0.0;

  // Get the beam instrumentation related objects
  art::ValidHandle<std::vector<beam::ProtoDUNEBeamEvent>> beamHandle = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventLabel);
  std::vector< art::Ptr<beam::ProtoDUNEBeamEvent> > beaminfo;

  if(!beamHandle.isValid()) return;

  art::fill_ptr_vector(beaminfo, beamHandle);

  if ( beaminfo.size() == 0 ) { std::cout << "WARNING: Beam event vector is empty." << std::endl; return;}
  if ( beaminfo.size() >  1 ) { std::cout << "WARNING: Beam event vector has size " << beaminfo.size() << std::endl;}

  // for beam particle selection strategy 
  // see p. 18 https://inspirehep.net/files/134bfcc0c74ca45b565e7e9514434c24
  int beamTrigFlag =  beaminfo[0]->GetTimingTrigger();
  if(fLogLevel>0) std::cout << "Beam event trigger is " << beamTrigFlag << std::endl;

  // beam inst. Time Of Flight
  fBeamTof = beaminfo[0]->GetTOF();
  if(fLogLevel>0) std::cout << "Beam event TOF is " << fBeamTof << std::endl;

  // beam inst. Cerenkov trigger
  fBeamCKov0 = beaminfo[0]->GetCKov0Status();
  fBeamCKov1 = beaminfo[0]->GetCKov1Status();
  if(fLogLevel>0) std::cout << "Beam event CKovs are " << fBeamCKov0 << " / " << fBeamCKov1 << std::endl;

  fMagnetCurrent = beaminfo[0]->GetMagnetCurrent();

  // beam inst. reconstructed momentum
  // see original implementation 
  // here: https://github.com/DUNE/duneprototypes/blob/develop/duneprototypes/Protodune/singlephase/BeamReco/BeamEvent_module.cc#L2490
  fBeamMomentum = GetRecoBeamMomentum(beaminfo[0]);

  if(fLogLevel>0) std::cout << "Beam event momentum is " << fBeamMomentum << std::endl;

  // Get the pandora related objects
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fSliceLabel, sliceHandle);

  if(!sliceHandle.isValid()) return;

  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  e.getByLabel(fPFParticleLabel, pfpHandle);

  if(!pfpHandle.isValid())   return;

  // Associate PFPs to tracks, showers and space points
  art::FindManyP<recob::Track>      pfpTrackAssoc(               pfpHandle, e, fTrackLabel);
  art::FindManyP<recob::Shower>     pfpShowerAssoc(              pfpHandle, e, fShowerLabel);
  art::FindManyP<anab::ParticleID>  pfpPIDAssoc(                 pfpHandle, e, fPIDLabel);
  art::FindManyP<recob::SpacePoint> pfpSPAssoc(                  pfpHandle, e, fSpacePointLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetaAssoc(pfpHandle, e, fPFParticleLabel);

  std::vector<art::Ptr<recob::Slice>> sliceVector;

  art::fill_ptr_vector(sliceVector, sliceHandle);
  art::FindManyP<recob::PFParticle> slicePFPAssoc(sliceHandle, e, fPFParticleLabel);

  // Get CNN hits scores
  anab::MVAReader<recob::Hit,4> * hitResults = 0x0;
  if (fUseMVA) hitResults = new anab::MVAReader<recob::Hit, 4>(e, "emtrkmichelid:emtrkmichel" );

  for(const art::Ptr<recob::Slice> &slice: sliceVector){
    std::vector<art::Ptr<recob::PFParticle>> slicePFPs(slicePFPAssoc.at(slice.key()));
 
    for(const art::Ptr<recob::PFParticle> &slicePFP : slicePFPs){
      //const bool isPrimary(slicePFP->IsPrimary());

      fPrimaryPdg        = -999;
      fPrimaryKey        = -999;
      fNPFParticles      = -999;
      fNPrimaryChildren  = -999;

      fPrimBeamCNNtrack  = -999;
      fPrimBeamCNNshower = -999;
      fPrimBeamCNNmichel = -999;

      fPrimBeamPFPTrackScore  = -999;
      fPrimBeamTrackPid       = -999;
      fPrimBeamTrackVertexX   = -999;
      fPrimBeamTrackVertexY   = -999;
      fPrimBeamTrackVertexZ   = -999;
      fPrimBeamTrackDirX      = -999;
      fPrimBeamTrackDirY      = -999;
      fPrimBeamTrackDirZ      = -999;
      fPrimBeamTrackLength    = -999;

      fPrimBeamShowerVertexX       = -999;
      fPrimBeamShowerVertexY       = -999;
      fPrimBeamShowerVertexZ       = -999;
      fPrimBeamShowerDirX          = -999;
      fPrimBeamShowerDirY          = -999;
      fPrimBeamShowerDirZ          = -999;
      fPrimBeamShowerLength        = -999;
      fPrimBeamShowerEnergy        = -999;
      fPrimBeamShowerdEdx          = -999;
      fPrimBeamShowerOpeningAngle  = -999;

      fPrimBeamX.clear();
      fPrimBeamY.clear();
      fPrimBeamZ.clear();

      fChildBeamPFPKey.clear();
      fChildBeamPFPPdg.clear();
      fChildBeamPFPChildren.clear();
     
      fChildBeamCNNtrack.clear();
      fChildBeamCNNshower.clear();
      fChildBeamCNNmichel.clear();

      fChildBeamPFPTrackScore.clear();

      fChildBeamTrackPid.clear();
      fChildBeamTrackVertexX.clear();
      fChildBeamTrackVertexY.clear();
      fChildBeamTrackVertexZ.clear();
      fChildBeamTrackDirX.clear();
      fChildBeamTrackDirY.clear();
      fChildBeamTrackDirZ.clear();
      fChildBeamTrackLength.clear();
      
      fChildBeamShowerVertexX.clear();
      fChildBeamShowerVertexY.clear();
      fChildBeamShowerVertexZ.clear();
      fChildBeamShowerDirX.clear();
      fChildBeamShowerDirY.clear();
      fChildBeamShowerDirZ.clear();
      fChildBeamShowerLength.clear();
      fChildBeamShowerEnergy.clear();
      fChildBeamShowerdEdx.clear();
      fChildBeamShowerOpeningAngle.clear();

      fChildBeamX.clear();
      fChildBeamY.clear();
      fChildBeamZ.clear();

      // Only consider primary particles to idenfy the incoming beam particle
      // only look at IsBeamParticle check for now
      // if(!isPrimary)                                           continue;

      // Sort out beam vs non-beam PFPs
      if(!pfpUtil.IsBeamParticle(*slicePFP, e, fPFParticleLabel)) continue;

      //////////////////////
      //
      // Starting with primary particles
      //

      // Skip primaires idendified as muons to save file size 
      if(fExcludeMuons && std::abs(slicePFP->PdgCode()) == 13) continue;
	
      if(fLogLevel>1) std::cout << "  found primary particle at slice " << slice.key() << " with PDG code: " << slicePFP->PdgCode() << std::endl;

      fPrimaryKey       = slice.key();
      fPrimaryPdg       = slicePFP->PdgCode(); 
      fNPFParticles     = slicePFPs.size();
      fNPrimaryChildren = slicePFP->NumDaughters();

      if( fUseMVA ){
	// Retrieve the primary particle CNN scores from collection plane hits
        const std::vector< art::Ptr< recob::Hit > > prim_hits
          = pfpUtil.GetPFParticleHitsFromPlane_Ptrs(*slicePFP, e, fPFParticleLabel, 2);

        fPrimBeamCNNtrack  = 0;
        fPrimBeamCNNshower = 0;
        fPrimBeamCNNmichel = 0;
       
	// Loop over all PFP hits from collection plane to compute the average score
	// - could add a charge weight to discard near noise hits
        for( size_t h = 0; h < prim_hits.size(); ++h ){
          std::array<float,4> cnn_out = hitResults->getOutput( prim_hits[h] );
          fPrimBeamCNNtrack  += cnn_out[ hitResults->getIndex("track") ]; 
          fPrimBeamCNNshower += cnn_out[ hitResults->getIndex("em") ]; 
          fPrimBeamCNNmichel += cnn_out[ hitResults->getIndex("michel") ]; 
        }
	fPrimBeamCNNtrack  /= prim_hits.size();
	fPrimBeamCNNshower /= prim_hits.size();
	fPrimBeamCNNmichel /= prim_hits.size();
      }

      // Retrieve primary space points
      std::vector<art::Ptr<recob::SpacePoint>> prim_spacepoints = pfpSPAssoc.at(slicePFP.key());
      
      for (const auto& sp : prim_spacepoints){
	// remove badly reconstructed space points      
	if(sp->XYZ()[2] < -200) continue;

	fPrimBeamX.push_back(sp->XYZ()[0]);
	fPrimBeamY.push_back(sp->XYZ()[1]);
	fPrimBeamZ.push_back(sp->XYZ()[2]);
      } // end of primary space points loop

      // Retrieve the track score for PFP metadata
      if (pfpMetaAssoc.isValid()){
	if (pfpMetaAssoc.at(slicePFP->Self()).size() > 0){
	  const larpandoraobj::PFParticleMetadata metaData = *((pfpMetaAssoc.at(slicePFP->Self())).at(0));
	  const std::map<std::string, float> &propMap = metaData.GetPropertiesMap();
	  if(propMap.count("TrackScore") > 0) fPrimBeamPFPTrackScore = propMap.at("TrackScore");
	  else fPrimBeamPFPTrackScore = -1;
	}
	else fPrimBeamPFPTrackScore = -1;
      }
      else fPrimBeamPFPTrackScore = -1;

      // Retrieve tracks from PFPs list
      std::vector<art::Ptr<recob::Track>> prim_tracks = pfpTrackAssoc.at(slicePFP.key());
      // Retrieve shower from PFPs list
      std::vector<art::Ptr<recob::Shower>> prim_showers = pfpShowerAssoc.at(slicePFP.key());

      if(prim_tracks.size() == 1){
	art::Ptr<recob::Track> track = prim_tracks.at(0);
	// distinguish track/shower in tuples
	// fChildBeamTrackFlag = 1);

        // Retrieve PFP PID from chi2 particle ID module	
	// see https://github.com/LArSoft/lardata/blob/develop/lardata/ArtDataHelper/Dumpers/DumpParticleIDs_module.cc#L110
	std::vector<anab::ParticleID> pids = trackUtil.GetRecoTrackPID(*track, e, fTrackLabel, fPIDLabel);

	// Loop over algorithms - should be one (Chi2)
	// for (const anab::sParticleIDAlgScores pid : pids) {
	for (const auto &pid : pids) {
	  auto scores = pid.ParticleIDAlgScores();
	  
	  for (const anab::sParticleIDAlgScores& score : scores) {
	    if(score.fPlaneMask != 2) continue;
	    fPrimBeamTrackPid = score.fAssumedPdg;
	  }
	}

	// store vertex infomation 
	fPrimBeamTrackVertexX = track->Vertex().X();
	fPrimBeamTrackVertexY = track->Vertex().Y();
	fPrimBeamTrackVertexZ = track->Vertex().Z();
	// Store direction information
	fPrimBeamTrackDirX = track->StartDirection().X();
	fPrimBeamTrackDirY = track->StartDirection().Y();
	fPrimBeamTrackDirZ = track->StartDirection().Z();
	// store length infomation 
	fPrimBeamTrackLength = track->Length();
      } // end of track condition - UseAllParticle flag set to true to record tracks and shower

      if(prim_showers.size() == 1){
	art::Ptr<recob::Shower> shower = prim_showers.at(0);
	// store vertex infomation 
	fPrimBeamShowerVertexX = shower->ShowerStart().X();
	fPrimBeamShowerVertexY = shower->ShowerStart().Y();
	fPrimBeamShowerVertexZ = shower->ShowerStart().Z();
	// Store direction information
	fPrimBeamShowerDirX = shower->Direction().X();
	fPrimBeamShowerDirY = shower->Direction().Y();
	fPrimBeamShowerDirZ = shower->Direction().Z();
	// store length infomation 
	fPrimBeamShowerLength = shower->Length();
	// store energy infomation 
	fPrimBeamShowerEnergy = shower->Energy()[2];
	// store energy loss infomation 
	fPrimBeamShowerdEdx = shower->dEdx()[2];
	// opening angle information
	fPrimBeamShowerOpeningAngle = shower->OpenAngle();
      } // end of shower condition

      //////////////////////
      //
      // Moving to daughter particles
      //

      std::vector<art::Ptr<recob::PFParticle>> beamPFPs(slicePFPAssoc.at(fPrimaryKey));

      // Loop over the PFPs of the daughter beam particle
      for(const art::Ptr<recob::PFParticle> &beamPFP : beamPFPs){

	fChildBeamPFPKey.push_back(beamPFP.key());
	fChildBeamPFPPdg.push_back(beamPFP->PdgCode());
	fChildBeamPFPChildren.push_back(beamPFP->NumDaughters());

	if( fUseMVA ){
	  // Retrieve the primary particle CNN scores from collection plane hits
	  const std::vector< art::Ptr< recob::Hit > > daughters_hits
	    = pfpUtil.GetPFParticleHitsFromPlane_Ptrs(*beamPFP, e, fPFParticleLabel, 2);

	  float tmpCNNtrack  = 0;
	  float tmpCNNshower = 0;
	  float tmpCNNmichel = 0;

   	  // Loop over all PFP hits from collection plane to compute the average score
          // - could add a charge weight to discard near noise hits
	  for( size_t h = 0; h < daughters_hits.size(); ++h ){
	    std::array<float,4> cnn_out = hitResults->getOutput( daughters_hits[h] );
	    tmpCNNtrack  += cnn_out[ hitResults->getIndex("track") ];
	    tmpCNNshower += cnn_out[ hitResults->getIndex("em") ];
	    tmpCNNmichel += cnn_out[ hitResults->getIndex("michel") ];
	  }
	  tmpCNNtrack  /= daughters_hits.size();
	  tmpCNNshower /= daughters_hits.size();
	  tmpCNNmichel /= daughters_hits.size();

	  fChildBeamCNNtrack.push_back(tmpCNNtrack);
	  fChildBeamCNNshower.push_back(tmpCNNshower);
	  fChildBeamCNNmichel.push_back(tmpCNNmichel);
	}

	// space points distribution
	std::vector<art::Ptr<recob::SpacePoint>> spacepoints = pfpSPAssoc.at(beamPFP.key());
	
	std::vector<float> tmpX, tmpY, tmpZ;
	for (const auto& sp : spacepoints){
	  // remove badly reconstructed space points	  
	  if(sp->XYZ()[2] < 0) continue;

	  tmpX.push_back(sp->XYZ()[0]);
	  tmpY.push_back(sp->XYZ()[1]);
	  tmpZ.push_back(sp->XYZ()[2]);
	} // end of daughters space points loop
	
	fChildBeamX.push_back(tmpX);
	fChildBeamY.push_back(tmpY);
	fChildBeamZ.push_back(tmpZ);

	// Retrieve the track score for PFP metadata
	if (pfpMetaAssoc.isValid()){
          if (pfpMetaAssoc.at(beamPFP->Self()).size() > 0){
            const larpandoraobj::PFParticleMetadata metaData = *((pfpMetaAssoc.at(beamPFP->Self())).at(0));
            const std::map<std::string, float> &propMap = metaData.GetPropertiesMap();
            if(propMap.count("TrackScore") > 0) fChildBeamPFPTrackScore.push_back(propMap.at("TrackScore"));
	    else fChildBeamPFPTrackScore.push_back(-1);
          }
	  else fChildBeamPFPTrackScore.push_back(-1);
        }
	else fChildBeamPFPTrackScore.push_back(-1);

	// Retrieve tracks from PFPs list
	std::vector<art::Ptr<recob::Track>> tracks = pfpTrackAssoc.at(beamPFP.key());
	// Retrieve shower from PFPs list
	std::vector<art::Ptr<recob::Shower>> showers = pfpShowerAssoc.at(beamPFP.key());

	if(tracks.size() == 1){
	  art::Ptr<recob::Track> track = tracks.at(0);
	  // distinguish track/shower in tuples

          std::vector<anab::ParticleID> pids = trackUtil.GetRecoTrackPID(*track, e, fTrackLabel, fPIDLabel);

	  int pdgScore = -999;

          // Loop over algorithms - should be one (Chi2)
          for (const auto &pid : pids) {
            auto scores = pid.ParticleIDAlgScores();

            for (const anab::sParticleIDAlgScores& score : scores) {
              if(score.fPlaneMask != 2) continue;
              pdgScore = score.fAssumedPdg;
            }
          }

	  fChildBeamTrackPid.push_back(pdgScore);

	  // store vertex infomation 
	  fChildBeamTrackVertexX.push_back(track->Vertex().X());
	  fChildBeamTrackVertexY.push_back(track->Vertex().Y());
	  fChildBeamTrackVertexZ.push_back(track->Vertex().Z());
          // Store direction information
	  fChildBeamTrackDirX.push_back(track->StartDirection().X());
	  fChildBeamTrackDirY.push_back(track->StartDirection().Y());
	  fChildBeamTrackDirZ.push_back(track->StartDirection().Z());
          // store length infomation 
          fChildBeamTrackLength.push_back(track->Length());
	} // end of track condition - UseAllParticle flag set to true to record tracks and shower
	
        if(showers.size() == 1){
	  art::Ptr<recob::Shower> shower = showers.at(0);
          // store vertex infomation 
	  fChildBeamShowerVertexX.push_back(shower->ShowerStart().X());
	  fChildBeamShowerVertexY.push_back(shower->ShowerStart().Y());
	  fChildBeamShowerVertexZ.push_back(shower->ShowerStart().Z());
          // Store direction information
	  fChildBeamShowerDirX.push_back(shower->Direction().X());
	  fChildBeamShowerDirY.push_back(shower->Direction().Y());
	  fChildBeamShowerDirZ.push_back(shower->Direction().Z());
          // store length infomation 
	  fChildBeamShowerLength.push_back(shower->Length());
          // store energy infomation 
	  fChildBeamShowerEnergy.push_back(shower->Energy()[2]);
          // store energy loss infomation 
	  fChildBeamShowerdEdx.push_back(shower->dEdx()[2]);
          // opening angle information
          fChildBeamShowerOpeningAngle.push_back(shower->OpenAngle());
	} // end of shower condition
      } // end of beamPFP loop
    
      // Filling only primary particles
      fTree->Fill();

    } // end of PFP loop
  } //end of slice loop
}

void pdvd::PDVDbeamSelection::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree","Pandora PFParticles");

  // event
  fTree->Branch("eventID",                  &fEventID);
  // beam instrumentation
  fTree->Branch("BeamTOF",                  &fBeamTof);
  fTree->Branch("BeamCKov0",                &fBeamCKov0);
  fTree->Branch("BeamCKov1",                &fBeamCKov1);
  fTree->Branch("BeamMomentum",             &fBeamMomentum);
  // pandora PFPs - primary 
  fTree->Branch("PrimaryKey",               &fPrimaryKey);
  fTree->Branch("PrimaryPdg",               &fPrimaryPdg);
  fTree->Branch("NPFParticles",             &fNPFParticles);
  fTree->Branch("NPrimaryChildren",         &fNPrimaryChildren);
  fTree->Branch("PrimBeamPFPCNNtrack",      &fPrimBeamCNNtrack);
  fTree->Branch("PrimBeamPFPCNNshower",     &fPrimBeamCNNshower);
  fTree->Branch("PrimBeamPFPCNNmichel",     &fPrimBeamCNNmichel);
  fTree->Branch("PrimBeamPFPTrackPid",      &fPrimBeamTrackPid);
  fTree->Branch("PrimBeamPFPTrackScore",    &fPrimBeamPFPTrackScore);
  fTree->Branch("PrimBeamPFPTrackVertexX",  &fPrimBeamTrackVertexX);
  fTree->Branch("PrimBeamPFPTrackVertexY",  &fPrimBeamTrackVertexY);
  fTree->Branch("PrimBeamPFPTrackVertexZ",  &fPrimBeamTrackVertexZ);
  fTree->Branch("PrimBeamPFPTrackDirX",     &fPrimBeamTrackDirX);
  fTree->Branch("PrimBeamPFPTrackDirY",     &fPrimBeamTrackDirY);
  fTree->Branch("PrimBeamPFPTrackDirZ",     &fPrimBeamTrackDirZ);
  fTree->Branch("PrimBeamPFPTrackLength",   &fPrimBeamTrackLength);
  fTree->Branch("PrimBeamPFPShowerVertexX", &fPrimBeamShowerVertexX);
  fTree->Branch("PrimBeamPFPShowerVertexY", &fPrimBeamShowerVertexY);
  fTree->Branch("PrimBeamPFPShowerVertexZ", &fPrimBeamShowerVertexZ);
  fTree->Branch("PrimBeamPFPShowerDirX",    &fPrimBeamShowerDirX);
  fTree->Branch("PrimBeamPFPShowerDirY",    &fPrimBeamShowerDirY);
  fTree->Branch("PrimBeamPFPShowerDirZ",    &fPrimBeamShowerDirZ);
  fTree->Branch("PrimBeamPFPShowerdEdx",    &fPrimBeamShowerdEdx);
  fTree->Branch("PrimBeamPFPShowerLength",  &fPrimBeamShowerLength);
  fTree->Branch("PrimBeamPFPShowerEnergy",  &fPrimBeamShowerEnergy);
  fTree->Branch("PrimBeamPFPX",             &fPrimBeamX);
  fTree->Branch("PrimBeamPFPY",             &fPrimBeamY);
  fTree->Branch("PrimBeamPFPZ",             &fPrimBeamZ);
  // pandora PFPs - children 
  fTree->Branch("ChildBeamPFPKey",          &fChildBeamPFPKey);
  fTree->Branch("ChildBeamPFPPdg",          &fChildBeamPFPPdg);
  fTree->Branch("ChildBeamPFPChildren",     &fChildBeamPFPChildren);
  fTree->Branch("ChildBeamPFPTrackScore",   &fChildBeamPFPTrackScore);
  fTree->Branch("ChildBeamPFPCNNtrack",     &fChildBeamCNNtrack);
  fTree->Branch("ChildBeamPFPCNNshower",    &fChildBeamCNNshower);
  fTree->Branch("ChildBeamPFPCNNmichel",    &fChildBeamCNNmichel);
  fTree->Branch("ChildBeamPFPTrackPid",     &fChildBeamTrackPid);
  fTree->Branch("ChildBeamPFPTrackVertexX", &fChildBeamTrackVertexX);
  fTree->Branch("ChildBeamPFPTrackVertexY", &fChildBeamTrackVertexY);
  fTree->Branch("ChildBeamPFPTrackVertexZ", &fChildBeamTrackVertexZ);
  fTree->Branch("ChildBeamPFPTrackDirX",    &fChildBeamTrackDirX);
  fTree->Branch("ChildBeamPFPTrackDirY",    &fChildBeamTrackDirY);
  fTree->Branch("ChildBeamPFPTrackDirZ",    &fChildBeamTrackDirZ);
  fTree->Branch("ChildBeamPFPTrackLength",  &fChildBeamTrackLength);
  fTree->Branch("ChildBeamPFPShowerVertexX",&fChildBeamShowerVertexX);
  fTree->Branch("ChildBeamPFPShowerVertexY",&fChildBeamShowerVertexY);
  fTree->Branch("ChildBeamPFPShowerVertexZ",&fChildBeamShowerVertexZ);
  fTree->Branch("ChildBeamPFPShowerDirX",   &fChildBeamShowerDirX);
  fTree->Branch("ChildBeamPFPShowerDirY",   &fChildBeamShowerDirY);
  fTree->Branch("ChildBeamPFPShowerDirZ",   &fChildBeamShowerDirZ);
  fTree->Branch("ChildBeamPFPShowerdEdx",   &fChildBeamShowerdEdx);
  fTree->Branch("ChildBeamPFPShowerLength", &fChildBeamShowerLength);
  fTree->Branch("ChildBeamPFPShowerEnergy", &fChildBeamShowerEnergy);
  fTree->Branch("ChildBeamPFPX",            &fChildBeamX);
  fTree->Branch("ChildBeamPFPY",            &fChildBeamY);
  fTree->Branch("ChildBeamPFPZ",            &fChildBeamZ);
}

void pdvd::PDVDbeamSelection::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(pdvd::PDVDbeamSelection)
