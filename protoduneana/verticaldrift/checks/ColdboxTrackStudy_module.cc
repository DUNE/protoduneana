////////////////////////////////////////////////////////////////////////
// Class:       ColdboxTrackStudy
// Plugin Type: analyzer (Unknown Unknown)
// File:        ColdboxTrackStudy_module.cc
//
// Generated at Tue Feb 20 07:41:24 2024 by Leila Haegel using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

// Default
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

// LarSoft
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/ServicePack.h" 
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// ROOT
#include "Math/ProbFunc.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TStyle.h"

// C++
using std::vector;
using std::string;

namespace pdvdana {
  class ColdboxTrackStudy;
}

class pdvdana::ColdboxTrackStudy : public art::EDAnalyzer {
public:
  explicit ColdboxTrackStudy(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  
  // Plugins should not be copied or assigned.
  ColdboxTrackStudy(ColdboxTrackStudy const&) = delete;
  ColdboxTrackStudy(ColdboxTrackStudy&&) = delete;
  ColdboxTrackStudy& operator=(ColdboxTrackStudy const&) = delete;
  ColdboxTrackStudy& operator=(ColdboxTrackStudy&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  // verbosity
  int      fLogLevel;

  // ROOT stored objects
  TTree*                 fEventTree;
  TTree*                 fTrackTree;

  // event information
  unsigned int fRun;
  int          fEventId;
  
  // track information
  string            fTrackModuleLabel;
  unsigned int      fNTracks;
  int               fTrackId;
  int               fThroughGoingTrack;
  float             fTrackLength;
  float             fTrackStartX;
  float             fTrackStartY;
  float             fTrackStartZ;
  float             fTrackVertexX;
  float             fTrackVertexY;
  float             fTrackVertexZ;
  float             fTrackEndX;
  float             fTrackEndY;
  float             fTrackEndZ;
  float             fTrackPhi;
  float             fTrackTheta;
  void              GetTrackAngles(const recob::Track& track, float &theta, float &phi) ;
    
  // track counter
  int nb_tracks = 0;
  int nb_track_throughgoing = 0 ;

  // calorimetry information
  string                 fCaloModuleLabel;
  vector<float>          fCaloRange;
  vector<float> fCaloHitX_plane;
  vector<float> fCaloHitX_plane0;
  vector<float> fCaloHitX_plane1;
  vector<float> fCaloHitX_plane2;
  vector<float> fCaloHitXcorr_plane;
  vector<float> fCaloHitXcorr_plane0;
  vector<float> fCaloHitXcorr_plane1;
  vector<float> fCaloHitXcorr_plane2;
  vector<vector<vector<float>>> fCaloHitXcorr_alltrks;
  vector<vector<float>> fCaloHitXcorr_plane_alltrks;
  vector<vector<float> > fCaloHitX;
  vector<vector<float> > fCaloHitXcorr;
  vector<float> fCaloHitY_plane;
  vector<float> fCaloHitY_plane0;
  vector<float> fCaloHitY_plane1;
  vector<float> fCaloHitY_plane2;
  vector<vector<float> > fCaloHitY;
  vector<float> fCaloHitZ_plane;
  vector<float> fCaloHitZ_plane0;
  vector<float> fCaloHitZ_plane1;
  vector<float> fCaloHitZ_plane2;
  vector<vector<float> > fCaloHitZ;
  vector<float> fCaloHitResRange_plane;
  vector<float> fCaloHitResRange_plane0;
  vector<float> fCaloHitResRange_plane1;
  vector<float> fCaloHitResRange_plane2;
  vector<vector<float> > fCaloHitResRange;
  vector<float> fCaloHitDqdx_plane;
  vector<float> fCaloHitDqdx_plane0;
  vector<float> fCaloHitDqdx_plane1;
  vector<float> fCaloHitDqdx_plane2;
  vector<vector<vector<float>>> fCaloHitDqdx_alltrks;
  vector<vector<float>> fCaloHitDqdx_plane_alltrks;
  vector<vector<float> > fCaloHitDqdx;
  vector<float> fCaloHitDedx_plane;
  vector<float> fCaloHitDedx_plane0;
  vector<float> fCaloHitDedx_plane1;
  vector<float> fCaloHitDedx_plane2;
  vector<vector<float> > fCaloHitDedx;
  vector<float> fCaloHitTrackPitch_plane;
  vector<float> fCaloHitTrackPitch_plane0;
  vector<float> fCaloHitTrackPitch_plane1;
  vector<float> fCaloHitTrackPitch_plane2;
  vector<vector<float> > fCaloHitTrackPitch;
  vector<int> fCaloHitTpIndices_plane;
  vector<int> fCaloHitTpIndices_plane0;
  vector<int> fCaloHitTpIndices_plane1;
  vector<int> fCaloHitTpIndices_plane2;
  vector<vector<int> >   fCaloHitTpIndices;

  // track hit information
  unsigned int           fNTrackHits;
  vector<int>            fTrackHitPlane;
  vector<int>            fTrackHitChannel;
  vector<int>            fTrackHitSignalType;
  vector<int>            fTrackHitStartTick;
  vector<int>            fTrackHitEndTick;
  vector<float>          fTrackHitPeakTime;
  vector<float>          fTrackHitPeakTimeSigma;
  vector<float>          fTrackHitIntegral;
  vector<float>          fTrackHitIntegralSigma;
  vector<float>          fTrackHitSummedADC;
  vector<float>          fTrackHitPeakAmplitude;
  vector<float>          fTrackHitPeakAmplitudeSigma;
  vector<float>          fTrackHitRMS;
  vector<float>          fTrackHitMultiplicity;
  vector<float>          fTrackHitGoodnessOfFit;
  vector<float>          fTrackHitLocalIndex;

  // space point information
  string                 fSpointModuleLabel;
  //unsigned int           fNSpoints;

  // hit information
  string                 fHitModuleLabel;
  unsigned int           fNHits;
  vector<int>            fHitPlane;
  vector<int>            fHitChannel;
  vector<int>            fHitSignalType;
  vector<int>            fHitStartTick;
  vector<int>            fHitEndTick;
  vector<float>          fHitPeakTime;
  vector<float>          fHitPeakTimeSigma;
  vector<float>          fHitIntegral;
  vector<float>          fHitIntegralSigma;
  vector<float>          fHitSummedADC;
  vector<float>          fHitPeakAmplitude;
  vector<float>          fHitPeakAmplitudeSigma;
  vector<float>          fHitRMS;
  vector<float>          fHitMultiplicity;
  vector<float>          fHitGoodnessOfFit;
  vector<float>          fHitLocalIndex;

  // wire information
  string                 fWireModuleLabel;
  unsigned int           fNWires;
  vector<int>            fWirePlane;
  vector<int>            fWireChannel;
  vector<int>            fWireNSignal;  
  vector<float>          fWireSignalROI;
  vector<vector<float>>  fWireSignalsROI;
  vector<vector<float>>  fWireSignals;
  
  // raw digit information
  string              fRawDigitModuleLabel;
  vector<int>         fRawDigitChannel;
  vector<int>         fRawDigitADCs;
  vector<vector<int>> fRawDigitsADCs;
  vector<float>       fRawDigitPedestal;
  
  // colbox-VD geometry
  const geo::Geometry* fCbvdGeom ;
  const geo::WireReadoutGeom* fWireReadoutGeom;
  unsigned int         fNPlanes;
  float                fCbvdXmin; 
  float                fCbvdXmax; 
  float                fCbvdYmin; 
  float                fCbvdYmax; 
  float                fCbvdZmin; 
  float                fCbvdZmax; 
  float                fCbvdDriftSpeed;
  float                fFiducialVolumeCutYZ;

  // options to save info
  float             fDqMinCut;
  float             fDqMaxCut;
  bool              fSaveAllChannelsInfo;
  bool              fSaveTrackHitsInfo;
  bool              fThroughGoingDownTrackOnly;
  bool CheckIfTrackIsThroughGoing(float trackStartVec[3], float trackEndVec[3]);    
  bool CheckCaloInfoValue(float caloValue, float range_min, float range_max, std::string caloValueName) ;
  bool CheckIfHitIsInFiducialVolume(float hitPosX, float hitPosY, float hitPosZ);    

  // constants 
  float  deg_to_rad = 3.14159/180.; 
  double conv_wcls_elec = 200. ;
  double elec_charge_fc = 0.000160217663 ;

  // dQ/dX fit
  int fFitDqdxNbins ;
  double fFitDqdxXmin ;  
  double fFitDqdxXmax ;
  vector<TF1*> fit_dqdx;
  vector<TF1*> fit_dqdx_corr;
  TF1* FitDqdx(TH1F* hist1d_dqdx, int i_plane, bool save) ;
    
  // dQ/dX vs time
  int fFitDqdxTimeNXbins ;
  int fFitDqdxTimeNYbins ;
  double fFitDqdxTimeXmin ;
  double fFitDqdxTimeXmax ;
  double fFitDqdxTimeYmin ;
  double fFitDqdxTimeYmax ;
  vector<double> e_lifetime ;
  vector<double> e_lifetime_error ;
  vector<TF1*> fit_dqdx_time;
  void FitDqdxTime(TH2F* hist2d_dQdX_time, TF1* fit_dqdx_time, int i_plane) ;

  // histograms and plots
  vector<TH1F*>    hist1d_dQdX;
  vector<TH1F*>    hist1d_dQdX_corr;
  vector<TH2F*>    hist2d_dQdX_dQ;
  vector<TH2F*>    hist2d_dQdX_phi;
  vector<TH2F*>    hist2d_dQdX_theta;
  vector<TH2F*>    hist2d_dQdX_time;
  vector<TH2F*>    hist2d_dQdX_yz;
  vector<TH2F*>    hist2d_dQdX_yx;
  vector<TH2F*>    hist2d_dQdX_zx;
  vector<TH2F*>    hist2d_y_z;

  vector<TCanvas*> canv_dQdX;
  vector<TCanvas*> canv_dQdX_fit;
  vector<TCanvas*> canv_dQdXcorr_fit;
  vector<TCanvas*> canv_dQdX_time;
  vector<TCanvas*> canv_dQdX_dQ;
  vector<TCanvas*> canv_dQdX_phi;
  vector<TCanvas*> canv_dQdX_theta;
  vector<TCanvas*> canv_y_z;
  vector<TCanvas*> canv_dQdX_yz;
  vector<TCanvas*> canv_dQdX_yx;
  vector<TCanvas*> canv_dQdX_zx;

  void MakePrettyTH1F(TH1F* hist);
  void Save1Hist1d(TCanvas* canv, TH1F* hist1, TString filename);
  void Save1Hist2d(TCanvas* canv, TH2F* hist, TString filename);
  void Save3Hist1d(TCanvas* canv, TH1F* hist1, TH1F* hist2, TH1F* hist3, TString filename);  
  void SaveDqdxFit(TCanvas* canv, TH1F* hist1, TF1* fct, TString filename);
  void SaveDqdxTimeFit(TCanvas* canv, TH2F* hist, TGraphErrors* graph, TF1* fct, TString filename);

};


double LandauGaussFct(double *x, double *par) {
   //See reference: https://root.cern/doc/master/langaus_8C.html
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
   // Numeric constants                                                    
   double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
   double mpshift  = -0.22278298;       // Landau maximum location
 
   // Control constants
   double np = 100.0;      // number of convolution steps
   double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

   // Variables
   double xx;
   double mpc;
   double fland;
   double sum = 0.0;
   double xlow,xupp;
   double step;
   double i;

   // MP shift correction
   mpc = par[1] - mpshift * par[0];

   // Range of convolution integral
   xlow = x[0] - sc * par[3];
   xupp = x[0] + sc * par[3];

   step = (xupp-xlow) / np;

   // Convolution integral of Landau and Gaussian by sum
   for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);

      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
   }
   
   double baseline = par[4]*(1-ROOT::Math::normal_cdf(double(x[0]),double(par[3]),double(par[1])));

   return (par[2] * step * sum * invsq2pi / par[3] + baseline);
}



/////////////////////////////////////////////////////////////////////////////
// Read the configuration options from the FHICL file
pdvdana::ColdboxTrackStudy::ColdboxTrackStudy(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fLogLevel(          	           p.get< int >("LogLevel")),
    fTrackModuleLabel(               p.get< std::string >("TrackModuleLabel")),
    fCaloModuleLabel(                p.get< std::string >("CaloModuleLabel")),
    fSpointModuleLabel(              p.get< std::string >("SpointModuleLabel")),
    fHitModuleLabel(                 p.get< std::string >("HitModuleLabel")),
    fWireModuleLabel(                p.get< std::string >("WireModuleLabel")),
    fRawDigitModuleLabel(            p.get< std::string >("RawDigitModuleLabel")),
    fFiducialVolumeCutYZ(            p.get< double >("FiducialVolumeCutYZ")),
    fDqMinCut(                       p.get< double >("DqMinCut")),
    fDqMaxCut(                       p.get< double >("DqMaxCut")),
    fSaveAllChannelsInfo(            p.get< bool >("SaveAllChannelsInfo")),
    fSaveTrackHitsInfo(              p.get< bool >("SaveTrackHitsInfo")),
    fThroughGoingDownTrackOnly(      p.get< bool >("ThroughGoingDownTrackOnly")),
    fFitDqdxNbins(                   p.get< double >("FitDqdxNbins")),
    fFitDqdxXmin(                    p.get< double >("FitDqdxXmin")),
    fFitDqdxXmax(                    p.get< double >("FitDqdxXmax")),
    fFitDqdxTimeNXbins(              p.get< double >("FitDqdxTimeNXbins")),
    fFitDqdxTimeNYbins(              p.get< double >("FitDqdxTimeNYbins")),
    fFitDqdxTimeXmin(                p.get< double >("FitDqdxTimeXmin")),
    fFitDqdxTimeXmax(                p.get< double >("FitDqdxTimeXmax")),
    fFitDqdxTimeYmin(                p.get< double >("FitDqdxTimeYmin")),
    fFitDqdxTimeYmax(                p.get< double >("FitDqdxTimeYmax"))
{
  // Calling detector geometry product 
  fCbvdGeom = &*art::ServiceHandle<geo::Geometry>();
  fWireReadoutGeom = &art::ServiceHandle<geo::WireReadout>()->Get();
}

bool pdvdana::ColdboxTrackStudy::CheckIfTrackIsThroughGoing(float trackStartVec[3], float trackEndVec[3])
{
  /////////////////////////////////////////////////////////////////////////////
  // Check if track is going downward through the coldbox

  bool trackIsThroughGoing = true ;

  // Select tracks starting & ending inside the detector in Y and Z
  if(trackStartVec[1] < fCbvdYmin || trackStartVec[1] > fCbvdYmax)   trackIsThroughGoing = false ;
  if(trackEndVec[1]   < fCbvdYmin || trackEndVec[1]   > fCbvdYmax)   trackIsThroughGoing = false ;
  if(trackStartVec[2] < fCbvdZmin || trackStartVec[2] > fCbvdZmax)   trackIsThroughGoing = false ;
  if(trackEndVec[2]   < fCbvdZmin || trackEndVec[2]   > fCbvdZmax)   trackIsThroughGoing = false ;
  /*  if(trackIsThroughGoing == false) {
    std::cout << "track Vec[1]: " << trackStartVec[1] << " - " << trackEndVec[1] << std::endl 
	      << " | fCbvdY edges: " << fCbvdYmin << " - " << fCbvdYmax << std::endl 
	      << "track Vec[2]: " << trackStartVec[2] << " - " << trackEndVec[2] << std::endl 
	      << " | fCbvdZ edges: " << fCbvdZmin << " - " << fCbvdZmax << std::endl 
	      << std::endl ;} */

  // Select tracks with x-length is at least 95% the height of the coldbox (and not more than 100%)
  // to be sure that the track cross the anode and the cathode
  if( fabs(trackEndVec[0] - trackStartVec[0]) < 0.95 * (fCbvdXmax-fCbvdXmin) )  trackIsThroughGoing = false ;
  if( fabs(trackEndVec[0] - trackStartVec[0]) > (fCbvdXmax-fCbvdXmin) )         trackIsThroughGoing = false ;
  /*  if(trackIsThroughGoing == false) {
    std::cout << "fabs(trackEndVec[0] - trackStartVec[0]) " << fabs(trackEndVec[0] - trackStartVec[0]) << std::endl 
	      << " | (fCbvdXmax-fCbvdXmin) ( + *0.95)" << (fCbvdXmax-fCbvdXmin) << " " << 0.95 * (fCbvdXmax-fCbvdXmin) 
	      << std::endl ; } */

  // Select tracks going downwards
  if( trackEndVec[0] > trackStartVec[0] ) trackIsThroughGoing = false ; 
  /* if(trackIsThroughGoing == false) {
    std::cout << " trackEndVec[0] > trackStartVec[0] " << trackEndVec[0] << "  " << trackStartVec[0]
    << std::endl ;  } */

  return trackIsThroughGoing ;
}

bool pdvdana::ColdboxTrackStudy::CheckCaloInfoValue(float caloValue, float range_min, float range_max, std::string caloValueName)
{
  bool caloValueOk = true ;

  if( (caloValue < range_min) || (caloValue > range_max) )
    {
      caloValueOk = false ;
      std::cout << "!!! " << caloValueName << " = " << caloValue
		<< " => outside range [" << range_min << " - " << range_max << "]" << std::endl ;
    }
  return caloValueOk;
} 

bool pdvdana::ColdboxTrackStudy::CheckIfHitIsInFiducialVolume(float hitPosX, float hitPosY, float hitPosZ)
{
  bool isInFiducialVolume = true ;
  if ( (hitPosX < 0) || (hitPosX > fCbvdXmax - fCbvdXmin) ) isInFiducialVolume = false ;
  if ( (hitPosY < fCbvdYmin + fFiducialVolumeCutYZ) || (hitPosY > fCbvdYmax - fFiducialVolumeCutYZ) ) isInFiducialVolume = false ;
  if ( (hitPosZ < fCbvdZmin + fFiducialVolumeCutYZ) || (hitPosZ > fCbvdZmax - fFiducialVolumeCutYZ) ) isInFiducialVolume = false ;
  return isInFiducialVolume;
}  

void pdvdana::ColdboxTrackStudy::GetTrackAngles(const recob::Track& track, float &theta, float &phi){
  // Theta is 90(180) deg. for horizontal(vertical down-going) tracks
  // Phi is 0(90) deg. for tracks along the beam(collection strips) direction

  float x_len = track.End().X()-track.Start().X();
  float y_len = track.End().Y()-track.Start().Y();
  float z_len = track.End().Z()-track.Start().Z();

  phi = 0;
  if(z_len != 0) phi = abs(atan(y_len/z_len));

  if(     y_len >= 0 && z_len > 0) phi =               phi  / deg_to_rad;
  else if(y_len >= 0 && z_len < 0) phi = ( TMath::Pi()-phi) / deg_to_rad;
  else if(y_len < 0  && z_len < 0) phi = (-TMath::Pi()+phi) / deg_to_rad;
  else if(y_len < 0  && z_len > 0) phi =              -phi  / deg_to_rad;

  theta = (TMath::Pi()/2. - atan((x_len)/sqrt(pow(y_len,2)+pow(z_len,2))))/ deg_to_rad;
}

void pdvdana::ColdboxTrackStudy::analyze(art::Event const& evt)
{
  /////////////////////////////////////////////////////////////////////////////
  // Get event data
  // (1 event = 1 time 4-ms time window)
  
  fEventId = evt.id().event();
  fRun = evt.run();
  if( fLogLevel >= 2 ) std::cout << "---------- Run " << fRun << " - Start analysing event " << fEventId << " ..." << std::endl;
  
  /////////////////////////////////////////////////////////////////////////////
  // Save event raw digits information if desired
  // (raw ADC counts measured on each channel)
  
  if(fSaveAllChannelsInfo) {
    auto const& rawDigitHandle = *evt.getValidHandle<vector<raw::RawDigit>>(fRawDigitModuleLabel);
    for (unsigned i_chan = 0; i_chan < rawDigitHandle.size() ; ++ i_chan) {
      fRawDigitChannel.push_back(rawDigitHandle[i_chan].Channel()) ;
      fRawDigitPedestal.push_back(rawDigitHandle[i_chan].GetPedestal()) ;
      fRawDigitADCs.clear();
      for (size_t i_tick = 0 ; i_tick < rawDigitHandle[i_chan].NADC() ; i_tick ++)
	fRawDigitADCs.push_back(rawDigitHandle[i_chan].ADC(i_tick)) ;
      fRawDigitsADCs.push_back(fRawDigitADCs) ;
    } // end of loop over raw digits
    if( fLogLevel >= 3 ) std::cout << "--- Number of rawDigit handles found:  " << rawDigitHandle.size() << std::endl;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  // Save event wire information if desired
  // (deconvoluted signals fitted with Gaussian per event)

  if(fSaveAllChannelsInfo) {
    auto const& wireHandle = *evt.getValidHandle<vector<recob::Wire>>(fWireModuleLabel);  
    fNWires = wireHandle.size() ;
    if( fLogLevel >= 3 ) std::cout << "--- Number of wire handles found:  " << fNWires  << std::endl;
    
    for (unsigned i_wire = 0; i_wire < wireHandle.size() ; ++ i_wire) {
      fWirePlane.push_back(wireHandle[i_wire].View()) ;
      fWireChannel.push_back(wireHandle[i_wire].Channel());
      fWireNSignal.push_back(wireHandle[i_wire].NSignal());
      fWireSignals.push_back(wireHandle[i_wire].Signal()) ;
      
      fWireSignalROI.clear() ;
      for (long unsigned int i=0; i< wireHandle[i_wire].SignalROI().size(); i++) 
	       fWireSignalROI.push_back(wireHandle[i_wire].SignalROI()[i]) ;
      fWireSignalsROI.push_back(fWireSignalROI);
    }
  }

  //////////////////////////////////////////////////////
  // Save event hit information

  // Get all the hits in this event
  auto const hitHandle = evt.getValidHandle<vector<recob::Hit>>(fHitModuleLabel);  
  fNHits = hitHandle->size() ;
  
  // Loop over hits to save the information
  for (unsigned i_hit = 0; i_hit < fNHits ; ++ i_hit) {
    const recob::Hit& hit = hitHandle->at(i_hit);
    fHitPlane.push_back(hit.View()) ;
    fHitChannel.push_back(hit.Channel()) ;
    fHitSignalType.push_back(hit.SignalType()) ;
    fHitStartTick.push_back(hit.StartTick()) ;
    fHitEndTick.push_back(hit.EndTick()) ;
    fHitPeakTime.push_back(hit.PeakTime()) ;
    fHitPeakTimeSigma.push_back(hit.SigmaPeakTime()) ;
    fHitIntegral.push_back(hit.Integral()) ;
    fHitIntegralSigma.push_back(hit.SigmaIntegral()) ;
    fHitSummedADC.push_back(hit.ROISummedADC()) ;
    fHitPeakAmplitude.push_back(hit.PeakAmplitude()) ;
    fHitPeakAmplitudeSigma.push_back(hit.SigmaPeakAmplitude()) ;
    fHitRMS.push_back(hit.RMS()) ;
    fHitMultiplicity.push_back(hit.Multiplicity()) ;
    fHitGoodnessOfFit.push_back(hit.GoodnessOfFit()) ;
    fHitLocalIndex.push_back(hit.LocalIndex()) ;
  } // end of loop over event hits

  // Access the space points associated to hits in this event (to be continued...)
  art::FindManyP<recob::SpacePoint> spacepointsHits(hitHandle, evt, fSpointModuleLabel);
  
  //////////////////////////////////////////////////////
  // Storing the information in a TTree (1 entry = 1 event)
  fEventTree -> Fill() ;    

  // Get all the tracks in this event 
  auto trackHandle   = evt.getValidHandle<vector<recob::Track>>(fTrackModuleLabel);
  fNTracks = trackHandle->size() ;
  if( fLogLevel >= 2 ) std::cout << "--- Number of tracks found:  " << fNTracks  << std::endl;
  
  // Access the hits associated to tracks in this event
  art::FindManyP<recob::Hit, recob::TrackHitMeta> tracksHits(trackHandle, evt, fTrackModuleLabel);
  int n_track_hits_total = 0 ;
  
  // Access the calorimetry information associated to tracks in this event
  art::FindManyP<anab::Calorimetry> tracksCalos(trackHandle, evt, fCaloModuleLabel) ; 
  
  //////////////////////////////////////////////////////
  // Loop over tracks
  
  for (unsigned i_trk = 0; i_trk < fNTracks ; ++ i_trk) {

    //////////////////////////////////////////////////////
    // Get track information

    // Track ID
    const recob::Track& track = trackHandle->at(i_trk);
    fTrackId = track.ID() ;
    if( fLogLevel >= 2 ) std::cout << "--- Start analyzing track " << fTrackId ; 
    nb_tracks ++ ;
    
    // Track geometrical information 
    fTrackLength = track.Length();
    fTrackStartX = track.Start().X() ;
    fTrackStartY = track.Start().Y() ;
    fTrackStartZ = track.Start().Z() ;
    fTrackVertexX = track.Vertex().X() ;
    fTrackVertexY = track.Vertex().Y() ;
    fTrackVertexZ = track.Vertex().Z() ;
    fTrackEndX = track.End().X() ;
    fTrackEndY = track.End().Y() ;
    fTrackEndZ = track.End().Z() ;

    // Check if track is compatible with cosmic track (through-going, downward) 
    // if the fhicl file requires a cosmic track and this one is not, pass
    float trackStartVec[3] = {fTrackStartX, fTrackStartY, fTrackStartZ} ;
    float trackEndVec[3] = {fTrackEndX, fTrackEndY, fTrackEndZ} ;
    bool trackIsThroughGoing = CheckIfTrackIsThroughGoing(trackStartVec, trackEndVec) ;
    fThroughGoingTrack  = int(trackIsThroughGoing) ;
 
    if(trackIsThroughGoing == true) {
      if( fLogLevel >= 3 ) std::cout << ": through-going, downward track found " << fThroughGoingTrack << std::endl ;
      nb_track_throughgoing ++ ;
    }
    else {
      if( fThroughGoingDownTrackOnly == true ) {
        	std::cout << std::endl ;
        	continue ; }
    }

    // Get tranck theta, phi angles
    GetTrackAngles(track, fTrackTheta, fTrackPhi) ;
    if (fTrackTheta < 100) std::cout << "/tmp theta " << fTrackTheta << std::endl;
    //std::cout << "/tmp theta, phi " << fTrackTheta << "  " <<  fTrackPhi << std::endl;

    //////////////////////////////////////////////////////
    // Get hits associated to track
  
    // This object get the hits associated to the track
    const std::vector<art::Ptr<recob::Hit>>& trackHits = tracksHits.at(i_trk);
    fNTrackHits = trackHits.size() ;
    if( fLogLevel >= 3 ) std::cout << "- " << fNTrackHits << " hits are associated to this track" << std::endl;  
    n_track_hits_total += fNTrackHits ;
    
    // Clear the vector of hit information
    fTrackHitPlane.clear() ;
    fTrackHitChannel.clear() ;
    fTrackHitSignalType.clear() ;
    fTrackHitStartTick.clear() ;
    fTrackHitEndTick.clear() ;
    fTrackHitPeakTime.clear() ;
    fTrackHitPeakTimeSigma.clear() ;
    fTrackHitIntegral.clear() ;
    fTrackHitIntegralSigma.clear() ;
    fTrackHitSummedADC.clear() ;
    fTrackHitPeakAmplitude.clear() ;
    fTrackHitPeakAmplitudeSigma.clear() ;
    fTrackHitRMS.clear() ;
    fTrackHitMultiplicity.clear() ;
    fTrackHitGoodnessOfFit.clear() ;
    fTrackHitLocalIndex.clear() ;
     
    // Loop over hits associated to tracks
    for (unsigned i_hit = 0; i_hit < fNTrackHits ; ++ i_hit) {
      
      // Save track hit information
      fTrackHitPlane.push_back(trackHits[i_hit] -> View()) ;
      fTrackHitChannel.push_back(trackHits[i_hit] -> Channel()) ;
      fTrackHitSignalType.push_back(trackHits[i_hit] -> SignalType()) ;
      fTrackHitStartTick.push_back(trackHits[i_hit] -> StartTick()) ;
      fTrackHitEndTick.push_back(trackHits[i_hit] -> EndTick()) ;
      fTrackHitPeakTime.push_back(trackHits[i_hit] -> PeakTime()) ;
      fTrackHitPeakTimeSigma.push_back(trackHits[i_hit] -> SigmaPeakTime()) ;
      fTrackHitIntegral.push_back(trackHits[i_hit] -> Integral()) ;
      fTrackHitIntegralSigma.push_back(trackHits[i_hit] -> SigmaIntegral()) ;
      fTrackHitSummedADC.push_back(trackHits[i_hit] -> ROISummedADC()) ;
      fTrackHitPeakAmplitude.push_back(trackHits[i_hit] -> PeakAmplitude()) ;
      fTrackHitPeakAmplitudeSigma.push_back(trackHits[i_hit] -> SigmaPeakAmplitude()) ;
      fTrackHitRMS.push_back(trackHits[i_hit] -> RMS()) ;
      fTrackHitMultiplicity.push_back(trackHits[i_hit] -> Multiplicity()) ;
      fTrackHitGoodnessOfFit.push_back(trackHits[i_hit] -> GoodnessOfFit()) ;
      fTrackHitLocalIndex.push_back(trackHits[i_hit] -> LocalIndex()) ;

      // Print hit information to screen
      if( fLogLevel >= 5 )       std::cout << "- Hit " << i_hit << ": " << " plane = " << trackHits[i_hit] -> View() << " ; signal type = " << trackHits[i_hit] -> SignalType() << " ; channel = " << trackHits[i_hit] -> Channel() << " ; wire ID = " << trackHits[i_hit] -> WireID() << " ; start tick = " << trackHits[i_hit] -> StartTick() << " ; end tick = " << trackHits[i_hit] -> EndTick() << " ; peak time = " << trackHits[i_hit] -> PeakTime() << " +/- " << trackHits[i_hit] -> SigmaPeakTime() << " tick units "<< " ; integral = " << trackHits[i_hit] -> Integral() << " +/- " << trackHits[i_hit] -> SigmaIntegral() << " tick x ADC units"<< " ; summed ADC = " << trackHits[i_hit] -> ROISummedADC() << " ; peak amplitude = " << trackHits[i_hit] -> PeakAmplitude() << " +/- " << trackHits[i_hit] -> SigmaPeakAmplitude() << " ADC units"<< " ; multiplicity = " << trackHits[i_hit] -> Multiplicity() << " ; goodness of fit = " << trackHits[i_hit] -> GoodnessOfFit() << " ; local index = " << trackHits[i_hit] -> LocalIndex() << std::endl ;
            

    } // end of loop over track hits
    
    //////////////////////////////////////////////////////
    // Get calorimetry hits associated to track
    
    std::vector<art::Ptr<anab::Calorimetry>> trackCaloHits = tracksCalos.at(i_trk);

    // Clear the track vectors
    fCaloHitX.clear() ;
    fCaloHitXcorr.clear() ;
    fCaloHitY.clear() ;
    fCaloHitZ.clear() ;
    fCaloHitDqdx.clear() ;
    fCaloHitDedx.clear() ;
    fCaloHitResRange.clear() ;
    fCaloHitTpIndices.clear() ;
    fCaloHitTrackPitch.clear() ;
    fCaloHitXcorr_plane_alltrks.clear() ;
    fCaloHitDqdx_plane_alltrks.clear() ;
       
    // Loop over planes associated to calo hits
    // Note: planes are filled in reversed order (plane 0 = collection plane / plane 1 = induction plane 2 / plane 2 = induction plane 1)
    for (unsigned i_plane = 0; i_plane < trackCaloHits.size() ; ++ i_plane) { 

      // Get the range
      if(fLogLevel > 3) CheckCaloInfoValue((trackCaloHits[i_plane]-> Range()), 0, 450, "trackCaloHits.Range") ;
      fCaloRange.push_back(trackCaloHits[i_plane] -> Range()) ;

      // Clear the plane vectors
      fCaloHitDqdx_plane.clear() ;
      fCaloHitDedx_plane.clear() ;
      fCaloHitX_plane.clear() ;
      fCaloHitXcorr_plane.clear() ;
      fCaloHitY_plane.clear() ;
      fCaloHitZ_plane.clear() ;
      fCaloHitResRange_plane.clear() ;
      fCaloHitTrackPitch_plane.clear() ;
      fCaloHitTpIndices_plane.clear() ;
      
      // Loop over calorimetric hit clusters
      for (unsigned i_pt = 0 ; i_pt < trackCaloHits[i_plane] -> dQdx().size() ; ++ i_pt) {
		
      	// Get the spatial location and check that it is inside the fiducial volume
      	fCaloHitX_plane.push_back((trackCaloHits[i_plane]-> XYZ())[i_pt].X());
        //fCaloHitXcorr_plane.push_back((trackCaloHits[i_plane]-> XYZ())[i_pt].X() - fTrackStartX);
        fCaloHitXcorr_plane.push_back(fTrackStartX - (trackCaloHits[i_plane]-> XYZ())[i_pt].X());
      	fCaloHitY_plane.push_back((trackCaloHits[i_plane]-> XYZ())[i_pt].Y());
      	fCaloHitZ_plane.push_back((trackCaloHits[i_plane]-> XYZ())[i_pt].Z());

        if(CheckIfHitIsInFiducialVolume( fTrackStartX - (trackCaloHits[i_plane]-> XYZ())[i_pt].X(), (trackCaloHits[i_plane]-> XYZ())[i_pt].Y(), (trackCaloHits[i_plane]-> XYZ())[i_pt].Z()) == false) continue ;

        // Get the residual range
        if ( CheckCaloInfoValue((trackCaloHits[i_plane]-> ResidualRange())[i_pt], 0, 500, "trackCaloHits.ResRange") == false) continue ;
        fCaloHitResRange_plane.push_back((trackCaloHits[i_plane]-> ResidualRange())[i_pt]) ;

        // Get the track pitch at hit location
        fCaloHitTrackPitch_plane.push_back((trackCaloHits[i_plane]->TrkPitchVec())[i_pt]) ;
        fCaloHitTrackPitch_plane.push_back((trackCaloHits[i_plane]->TpIndices())[i_pt]) ;

      	// Get the dQ/dX
      	if ( CheckCaloInfoValue( (trackCaloHits[i_plane]-> dQdx())[i_pt], 0, 10000, "trackCaloHits.dQdx") == false) continue ;
        double dqdx = (trackCaloHits[i_plane]-> dQdx())[i_pt] * conv_wcls_elec * elec_charge_fc ;
        double dq = dqdx * (trackCaloHits[i_plane]->TrkPitchVec())[i_pt] ;
        
        hist2d_dQdX_dQ[i_plane] -> Fill(dq, dqdx) ;
        
        if (dq < fDqMinCut || dq > fDqMaxCut) continue ; //{
         // std::cout << "/tmp dq out of bounds | dqdx: " << dqdx << " | pitch: " << (trackCaloHits[i_plane]->TrkPitchVec())[i_pt] << " | dq: " << dq << std::endl ;
         // continue ;}

      	fCaloHitDqdx_plane.push_back(dqdx) ;
      	hist1d_dQdX[i_plane] -> Fill(dqdx) ;
        /*
        std::cout << "/tmp plane  " << i_plane << " | hit " << i_pt 
                  << " | x " << fTrackStartX - (trackCaloHits[i_plane]-> XYZ())[i_pt].X()
                  << " | dqdx " << (trackCaloHits[i_plane]-> dQdx())[i_pt] * conv_wcls_elec * elec_charge_fc 
                  << " | pitch " << (trackCaloHits[i_plane]->TrkPitchVec())[i_pt]
                  << " | dqdx*pitch " << (trackCaloHits[i_plane]-> dQdx())[i_pt] * conv_wcls_elec * elec_charge_fc * (trackCaloHits[i_plane]->TrkPitchVec())[i_pt]
                  << std::endl ; */

        // Get the dQ/dX vs time
        double driftTime = (fTrackStartX - (trackCaloHits[i_plane]-> XYZ())[i_pt].X()) / fCbvdDriftSpeed ;
      	double driftDqdx = (trackCaloHits[i_plane] -> dQdx())[i_pt] * conv_wcls_elec * elec_charge_fc ;
      	hist2d_dQdX_time[i_plane] -> Fill(driftTime, driftDqdx) ;
        hist2d_dQdX_phi[i_plane] -> Fill(fTrackPhi, dqdx) ;
        hist2d_dQdX_theta[i_plane] -> Fill(fTrackTheta, dqdx) ;
        hist2d_dQdX_yz[i_plane]->Fill((trackCaloHits[i_plane]-> XYZ())[i_pt].Y(), (trackCaloHits[i_plane]-> XYZ())[i_pt].Z(), dqdx);
        hist2d_dQdX_yx[i_plane]->Fill((trackCaloHits[i_plane]-> XYZ())[i_pt].Y(), (trackCaloHits[i_plane]-> XYZ())[i_pt].X(), dqdx);
        hist2d_dQdX_zx[i_plane]->Fill((trackCaloHits[i_plane]-> XYZ())[i_pt].Z(), (trackCaloHits[i_plane]-> XYZ())[i_pt].X(), dqdx);
        hist2d_y_z[i_plane] -> Fill((trackCaloHits[i_plane]-> XYZ())[i_pt].Y(), (trackCaloHits[i_plane]-> XYZ())[i_pt].Z()) ;

      	// Get the dE/dX
      	if ( CheckCaloInfoValue((trackCaloHits[i_plane]-> dEdx())[i_pt], 0, 10000, "trackCaloHits.dEdx") == false) continue ;
     	fCaloHitDedx_plane.push_back((trackCaloHits[i_plane]-> dEdx())[i_pt]) ;
      }

      // Fill the track vector with calorimetric information
        fCaloHitX.push_back(fCaloHitX_plane) ;
        fCaloHitXcorr.push_back(fCaloHitXcorr_plane) ;
        fCaloHitY.push_back(fCaloHitY_plane) ;
        fCaloHitZ.push_back(fCaloHitZ_plane) ;
        fCaloHitDqdx.push_back(fCaloHitDqdx_plane) ;
        fCaloHitDedx.push_back(fCaloHitDedx_plane) ;
        fCaloHitResRange.push_back(fCaloHitResRange_plane) ;
        fCaloHitTpIndices.push_back(fCaloHitTpIndices_plane) ;
        fCaloHitTrackPitch.push_back(fCaloHitTrackPitch_plane) ;

        fCaloHitXcorr_plane_alltrks.push_back(fCaloHitXcorr_plane) ;
        fCaloHitDqdx_plane_alltrks.push_back(fCaloHitDqdx_plane) ;

    } // end of loop over track calo planes
       
    //////////////////////////////////////////////////////
    // Storing the information in a TTree (1 entry = 1 track)
    
    // Ceck if we save only through-going tracks
    if(fThroughGoingDownTrackOnly == true && trackIsThroughGoing == false)   continue ;

    // Separate the calorimetric information by plane
    fCaloHitX_plane0 = fCaloHitX[2] ;
    fCaloHitX_plane1 = fCaloHitX[1] ;
    fCaloHitX_plane2 = fCaloHitX[0] ;
    fCaloHitXcorr_plane0 = fCaloHitXcorr[2] ;
    fCaloHitXcorr_plane1 = fCaloHitXcorr[1] ;
    fCaloHitXcorr_plane2 = fCaloHitXcorr[0] ;

    fCaloHitY_plane0 = fCaloHitY[2] ;
    fCaloHitY_plane1 = fCaloHitY[1] ;
    fCaloHitY_plane2 = fCaloHitY[0] ;
    fCaloHitZ_plane0 = fCaloHitZ[2] ;
    fCaloHitZ_plane1 = fCaloHitZ[1] ;
    fCaloHitZ_plane2 = fCaloHitZ[0] ;
    fCaloHitDqdx_plane0 = fCaloHitDqdx[2] ;
    fCaloHitDqdx_plane1 = fCaloHitDqdx[1] ;
    fCaloHitDqdx_plane2 = fCaloHitDqdx[0] ;
    fCaloHitDedx_plane0 = fCaloHitDedx[2] ;
    fCaloHitDedx_plane1 = fCaloHitDedx[1] ;
    fCaloHitDedx_plane2 = fCaloHitDedx[0] ;

    fCaloHitResRange_plane0 = fCaloHitResRange[2] ;
    fCaloHitResRange_plane1 = fCaloHitResRange[1] ;
    fCaloHitResRange_plane2 = fCaloHitResRange[0] ;
    fCaloHitTrackPitch_plane0 = fCaloHitTrackPitch[2] ;
    fCaloHitTrackPitch_plane1 = fCaloHitTrackPitch[1] ;
    fCaloHitTrackPitch_plane2 = fCaloHitTrackPitch[0] ;
    fCaloHitTpIndices_plane0 = fCaloHitTpIndices[2] ;
    fCaloHitTpIndices_plane1 = fCaloHitTpIndices[1] ;
    fCaloHitTpIndices_plane2 = fCaloHitTpIndices[0] ;
    fTrackTree -> Fill() ;    

    // Save those vector for electron lifetime correction
    fCaloHitXcorr_alltrks.push_back(fCaloHitXcorr_plane_alltrks) ; 
    fCaloHitDqdx_alltrks.push_back(fCaloHitDqdx_plane_alltrks) ; 

  } // end of loop over tracks

  // Print the total number of hits, and how many were associated to tracks
  if( fLogLevel >= 3 && fThroughGoingDownTrackOnly == false ) std::cout << "------ " << hitHandle->size() << " hits were found in this event, of which " << n_track_hits_total << " were associated to tracks. " << std::endl << std::endl ;

}


void pdvdana::ColdboxTrackStudy::beginJob()
{
  //////////////////////////////////////////////////////
  // Read FHiCL file options
  
  std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl
	    << "<> Starting ColdboxTrackStudy with options: " << std::endl 
	    << "   > save all channels info : " << fSaveAllChannelsInfo << std::endl
	    << "   > select only through-going downward tracks: " << fThroughGoingDownTrackOnly << std::endl 
	    << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl ;

  /////////////////////////////////////////////////////////////////////////////
  // Find the coldbox edges

  // Get the coldbox TPC edges (1 TPC = 1 CRU = 1/2 CRP)
  unsigned int NCbvdTpc = fCbvdGeom->NTPC() ;
  
  std::vector <float> cbvdTpcXmin(NCbvdTpc) ;
  std::vector <float> cbvdTpcXmax(NCbvdTpc) ;
  std::vector <float> cbvdTpcYmin(NCbvdTpc) ;
  std::vector <float> cbvdTpcYmax(NCbvdTpc) ;
  std::vector <float> cbvdTpcZmin(NCbvdTpc) ;
  std::vector <float> cbvdTpcZmax(NCbvdTpc) ;

  for(unsigned i_tpc=0 ; i_tpc < NCbvdTpc ; i_tpc ++){
    geo::TPCID cbvdTpcId{0, i_tpc};
    
    cbvdTpcXmin[i_tpc] = fCbvdGeom->TPC(cbvdTpcId).BoundingBox().MinX();
    cbvdTpcXmax[i_tpc] = fCbvdGeom->TPC(cbvdTpcId).BoundingBox().MaxX();
    cbvdTpcYmin[i_tpc] = fCbvdGeom->TPC(cbvdTpcId).BoundingBox().MinY();
    cbvdTpcYmax[i_tpc] = fCbvdGeom->TPC(cbvdTpcId).BoundingBox().MaxY();
    cbvdTpcZmin[i_tpc] = fCbvdGeom->TPC(cbvdTpcId).BoundingBox().MinZ();
    cbvdTpcZmax[i_tpc] = fCbvdGeom->TPC(cbvdTpcId).BoundingBox().MaxZ();
  }

  // Get the coldbox edges
  fCbvdXmin = *std::min_element(cbvdTpcXmin.begin(), cbvdTpcXmin.end());
  fCbvdXmax = *std::max_element(cbvdTpcXmax.begin(), cbvdTpcXmax.end());
  fCbvdYmin = *std::min_element(cbvdTpcYmin.begin(), cbvdTpcYmin.end());
  fCbvdYmax = *std::max_element(cbvdTpcYmax.begin(), cbvdTpcYmax.end());
  fCbvdZmin = *std::min_element(cbvdTpcZmin.begin(), cbvdTpcZmin.end());
  fCbvdZmax = *std::max_element(cbvdTpcZmax.begin(), cbvdTpcZmax.end());
  
  if( fLogLevel >= 3 ) std::cout << "<> Coldbox edges: " << std::endl 
				 << " > X: " << fCbvdXmin << " - " << fCbvdXmax << std::endl 
				 << " > Y: " << fCbvdYmin << " - " << fCbvdYmax << std::endl 
				 << " > Z: " << fCbvdZmin << " - " << fCbvdZmax << std::endl 
				 << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl ;				
  
  // Get the drift velocity - now hard coded, but need to fix the seg fault
  //  auto const cbvdPpties = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  // fCbvdDriftSpeed = cbvdPpties.DriftVelocity();
  auto const cbvdPpties   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  fCbvdDriftSpeed = cbvdPpties.DriftVelocity();
  std::cout << "fCbvdDriftSpeed " <<  fCbvdDriftSpeed << std::endl ;
//  fCbvdDriftSpeed = 0.143849 ;

  if( fLogLevel >= 3 ) std::cout << "<> Coldbox drift speed: " << fCbvdDriftSpeed << " cm/us" << std::endl;
 
  //////////////////////////////////////////////////////
  // Open ROOT output file and set tree & branches

  // Attach stored objects to output file
  art::ServiceHandle<art::TFileService> tfs;

  // Create output trees
  fEventTree = tfs->make<TTree>("eventTree","Store information about events" );
  fTrackTree = tfs->make<TTree>("trackTree","Store information about tracks" );

  // Get the number of anode planes to resize some vectors
  fNPlanes = fWireReadoutGeom->Nplanes();
  if( fNPlanes != 3) std::cout << "!!! Warning! Number of planes found is: " << fNPlanes << std::endl ;

  //////////////////////////////////////////////////////
  // Add branches to the event tree

  fEventTree->Branch("EventId",          &fEventId,      "EventId/i"   );
  fTrackTree->Branch("NTracks",          &fNTracks,      "NTracks/i"   );
  if(fSaveAllChannelsInfo) {
    //--- raw digits
    fEventTree->Branch("RawDigitChannel",  &fRawDigitChannel) ;
    fEventTree->Branch("RawDigitPedestal", &fRawDigitPedestal) ;
    fEventTree->Branch("RawDigitsADCs",    &fRawDigitsADCs) ;
    //--- wires 
    fEventTree->Branch("NWires",          &fNWires) ;
    fEventTree->Branch("WirePlane",       &fWirePlane) ;
    fEventTree->Branch("WireChannel",     &fWireChannel) ;
    fEventTree->Branch("WireNSignal",     &fWireNSignal);
    fEventTree->Branch("WireSignalsROI",  &fWireSignalsROI);
    fEventTree->Branch("WireSignals",     &fWireSignals);
    //--- hits 
    fEventTree->Branch("NHits",            &fNHits,        "NHits/i"   );
    fEventTree->Branch("HitPlane",         &fHitPlane);
    fEventTree->Branch("HitChannel",       &fHitChannel);
    fEventTree->Branch("HitSignalType",    &fHitSignalType);
    fEventTree->Branch("HitStartTick",     &fHitStartTick);
    fEventTree->Branch("HitEndTick",       &fHitEndTick);
    fEventTree->Branch("HitPeakTime",      &fHitPeakTime);
    fEventTree->Branch("HitPeakTimeSigma", &fHitPeakTimeSigma);
    fEventTree->Branch("HitIntegral",      &fHitIntegral);
    fEventTree->Branch("HitIntegralSigma", &fHitIntegralSigma);
    fEventTree->Branch("HitSummedADC",     &fHitSummedADC);
    fEventTree->Branch("HitPeakAmplitude", &fHitPeakAmplitude);
    fEventTree->Branch("HitPeakAmplitudeSigma", &fHitPeakAmplitudeSigma);
    fEventTree->Branch("HitRMS",           &fHitRMS);
    fEventTree->Branch("HitMultiplicity",  &fHitMultiplicity);
    fEventTree->Branch("HitGoodnessOfFit", &fHitGoodnessOfFit);
    fEventTree->Branch("HitLocalIndex",    &fHitLocalIndex);
  }

  //////////////////////////////////////////////////////
  // Add branches to the track tree

  fTrackTree->Branch("EventId",      &fEventId,      "EventId/i"   );
  //--- track information 
  fTrackTree->Branch("NTracks",      &fNTracks,      "NTracks/i"   );
  fTrackTree->Branch("TrackId",      &fTrackId,      "TrackId/i"   );
  fTrackTree->Branch("TrackThroughGoing", &fThroughGoingTrack, "ThroughGoingTrack/i"   );
  fTrackTree->Branch("TrackStartX",  &fTrackStartX,  "TrackStartX/F" ); 
  fTrackTree->Branch("TrackStartY",  &fTrackStartY,  "TrackStartY/F" ); 
  fTrackTree->Branch("TrackStartZ",  &fTrackStartZ,  "TrackStartZ/F" ); 
  fTrackTree->Branch("TrackVertexX", &fTrackVertexX, "TrackVertexX/F" ); 
  fTrackTree->Branch("TrackVertexY", &fTrackVertexY, "TrackVertexY/F" ); 
  fTrackTree->Branch("TrackVertexZ", &fTrackVertexZ, "TrackVertexZ/F" ); 
  fTrackTree->Branch("TrackEndX",    &fTrackEndX,    "TrackEndX/F" ); 
  fTrackTree->Branch("TrackEndY",    &fTrackEndY,    "TrackEndY/F" ); 
  fTrackTree->Branch("TrackEndZ",    &fTrackEndZ,    "TrackEndZ/F" ); 
  fTrackTree->Branch("TrackPhi",     &fTrackPhi,     "TrackPhi/F" ); 
  fTrackTree->Branch("TrackTheta",   &fTrackTheta,   "TrackTheta/F" ); 
    //--- tracks hits
  if(fSaveTrackHitsInfo) {
    fTrackTree->Branch("NTrackHits",   &fNTrackHits,   "NTrackHits/i"   );
    fTrackTree->Branch("TrackHitPlane",         &fTrackHitPlane);
    fTrackTree->Branch("TrackHitChannel",       &fTrackHitChannel);
    fTrackTree->Branch("TrackHitSignalType",    &fTrackHitSignalType);
    fTrackTree->Branch("TrackHitStartTick",     &fTrackHitStartTick);
    fTrackTree->Branch("TrackHitEndTick",       &fTrackHitEndTick);
    fTrackTree->Branch("TrackHitPeakTime",      &fTrackHitPeakTime);
    fTrackTree->Branch("TrackHitPeakTimeSigma", &fTrackHitPeakTimeSigma);
    fTrackTree->Branch("TrackHitIntegral",      &fTrackHitIntegral);
    fTrackTree->Branch("TrackHitIntegralSigma", &fTrackHitIntegralSigma);
    fTrackTree->Branch("TrackHitSummedADC",     &fTrackHitSummedADC);
    fTrackTree->Branch("TrackHitPeakAmplitude", &fTrackHitPeakAmplitude);
    fTrackTree->Branch("TrackHitPeakAmplitudeSigma", &fTrackHitPeakAmplitudeSigma);
    fTrackTree->Branch("TrackHitRMS",           &fTrackHitRMS);
    fTrackTree->Branch("TrackHitMultiplicity",  &fTrackHitMultiplicity);
    fTrackTree->Branch("TrackHitGoodnessOfFit", &fTrackHitGoodnessOfFit);
    fTrackTree->Branch("TrackHitLocalIndex",    &fTrackHitLocalIndex);
    //--- calorimetric info
    fTrackTree->Branch("TrackCaloRange",        &fCaloRange) ;
    fTrackTree->Branch("TrackCaloHitX_plane0",  &fCaloHitX_plane0) ;
    fTrackTree->Branch("TrackCaloHitX_plane1",  &fCaloHitX_plane1) ;
    fTrackTree->Branch("TrackCaloHitX_plane2",  &fCaloHitX_plane2) ;
    fTrackTree->Branch("TrackCaloHitXcorr_plane0",  &fCaloHitXcorr_plane0) ;
    fTrackTree->Branch("TrackCaloHitXcorr_plane1",  &fCaloHitXcorr_plane1) ;
    fTrackTree->Branch("TrackCaloHitXcorr_plane2",  &fCaloHitXcorr_plane2) ;
    fTrackTree->Branch("TrackCaloHitY_plane0",  &fCaloHitY_plane0) ;
    fTrackTree->Branch("TrackCaloHitY_plane1",  &fCaloHitY_plane1) ;
    fTrackTree->Branch("TrackCaloHitY_plane2",  &fCaloHitY_plane2) ;
    fTrackTree->Branch("TrackCaloHitZ_plane0",  &fCaloHitZ_plane0) ;
    fTrackTree->Branch("TrackCaloHitZ_plane1",  &fCaloHitZ_plane1) ;
    fTrackTree->Branch("TrackCaloHitZ_plane2",  &fCaloHitZ_plane2) ;
    fTrackTree->Branch("TrackCaloHitDqdx_plane0",  &fCaloHitDqdx_plane0) ;
    fTrackTree->Branch("TrackCaloHitDqdx_plane1",  &fCaloHitDqdx_plane1) ;
    fTrackTree->Branch("TrackCaloHitDqdx_plane2",  &fCaloHitDqdx_plane2) ;
    fTrackTree->Branch("TrackCaloHitDedx_plane0",  &fCaloHitDedx_plane0) ;
    fTrackTree->Branch("TrackCaloHitDedx_plane1",  &fCaloHitDedx_plane1) ;
    fTrackTree->Branch("TrackCaloHitDedx_plane2",  &fCaloHitDedx_plane2) ;
    fTrackTree->Branch("TrackCaloHitResRange_plane0",  &fCaloHitResRange_plane0) ;
    fTrackTree->Branch("TrackCaloHitResRange_plane1",  &fCaloHitResRange_plane1) ;
    fTrackTree->Branch("TrackCaloHitResRange_plane2",  &fCaloHitResRange_plane2) ;
    fTrackTree->Branch("TrackCaloHitTrkPitch_plane0",  &fCaloHitTrackPitch_plane0) ;
    fTrackTree->Branch("TrackCaloHitTrkPitch_plane1",  &fCaloHitTrackPitch_plane1) ;
    fTrackTree->Branch("TrackCaloHitTrkPitch_plane2",  &fCaloHitTrackPitch_plane2) ;
    fTrackTree->Branch("TrackCaloHitTpIndices_plane0",  &fCaloHitTpIndices_plane0) ;
    fTrackTree->Branch("TrackCaloHitTpIndices_plane1",  &fCaloHitTpIndices_plane1) ;
    fTrackTree->Branch("TrackCaloHitTpIndices_plane2",  &fCaloHitTpIndices_plane2) ;
  }

  //--- create histograms

  hist1d_dQdX.resize(fNPlanes);
  hist1d_dQdX_corr.resize(fNPlanes);
  hist2d_dQdX_time.resize(fNPlanes); 
  hist2d_dQdX_dQ.resize(fNPlanes); 
  hist2d_dQdX_phi.resize(fNPlanes); 
  hist2d_dQdX_theta.resize(fNPlanes); 
  hist2d_y_z.resize(fNPlanes); 
  hist2d_dQdX_yz.resize(fNPlanes); 
  hist2d_dQdX_yx.resize(fNPlanes); 
  hist2d_dQdX_zx.resize(fNPlanes); 
  
  for (unsigned i_plane = 0; i_plane < fNPlanes; i_plane++) {
    hist1d_dQdX[i_plane] = tfs->make<TH1F>(Form("hist1d_dQdX_plane%d", -i_plane+2), ";dQ/dX [fC/cm];Counts",  fFitDqdxNbins, 0, 20);
    hist1d_dQdX_corr[i_plane] = tfs->make<TH1F>(Form("hist1d_dQdX_corr_plane%d", -i_plane+2), ";dQ/dX [fC/cm];Counts",  fFitDqdxNbins, 0, 20);
    hist2d_dQdX_time[i_plane] = tfs->make<TH2F>(Form("hist2d_dQdX_time_plane%d", -i_plane+2), ";drift time [#mu s]; dQ/dX [fC/cm]; Counts", fFitDqdxTimeNXbins, 0, 160, fFitDqdxTimeNYbins, 0, 20);
    hist2d_dQdX_dQ[i_plane] = tfs->make<TH2F>(Form("hist2d_dQdX_dQ_plane%d", -i_plane+2), ";dQ [fC]; dQ/dX [fC/cm]; Counts", 200, 0., 20., 200, 0., 20.);
    hist2d_dQdX_phi[i_plane] = tfs->make<TH2F>(Form("hist2d_dQdX_phi_plane%d", -i_plane+2), ";track #phi [#circ]; dQ/dX [fC/cm]; Counts", 180, 0., 180., 200, 0., 20.);
    hist2d_dQdX_theta[i_plane] = tfs->make<TH2F>(Form("hist2d_dQdX_theta_plane%d", -i_plane+2), ";track #theta [#circ]; dQ/dX [fC/cm]; Counts", 90, 90., 180., 200, 0., 20.);
    hist2d_dQdX_yz[i_plane] = tfs->make<TH2F>(Form("hist2d_dQdX_yz_plane%d", -i_plane+2), ";y [cm]; z [cm]; dQ/dX [fC/cm]", 340, -170, 170., 300, 0., 300.);
    hist2d_dQdX_yx[i_plane] = tfs->make<TH2F>(Form("hist2d_dQdX_yx_plane%d", -i_plane+2), ";y [cm]; x [cm]; dQ/dX [fC/cm]", 340, -170, 170., 24, -12., 12.);
    hist2d_dQdX_zx[i_plane] = tfs->make<TH2F>(Form("hist2d_dQdX_zx_plane%d", -i_plane+2), ";z [cm]; x [cm]; dQ/dX [fC/cm]", 300, 0., 300., 24, -12., 12.);
    hist2d_y_z[i_plane] = tfs->make<TH2F>(Form("hist2d_y_z_plane%d", -i_plane+2), ";y [cm]; z [cm]; Counts", 340, -170, 170., 300, 0., 300.);
  }
}

void pdvdana::ColdboxTrackStudy::endJob()
{
  // Objects are saved in the output ROOT file
  art::ServiceHandle<art::TFileService> tfs;

  fit_dqdx.resize(fNPlanes) ;
  fit_dqdx_corr.resize(fNPlanes) ;
  fit_dqdx_time.resize(fNPlanes) ;

  canv_dQdX.resize(fNPlanes);
  canv_dQdX_fit.resize(fNPlanes);
  canv_dQdXcorr_fit.resize(fNPlanes);
  canv_dQdX_time.resize(fNPlanes);
  canv_dQdX_dQ.resize(fNPlanes);
  canv_dQdX_phi.resize(fNPlanes);
  canv_dQdX_theta.resize(fNPlanes);
  canv_dQdX_yz.resize(fNPlanes);
  canv_dQdX_yx.resize(fNPlanes);
  canv_dQdX_zx.resize(fNPlanes);
  canv_y_z.resize(fNPlanes);

  // Analyse the calorimetry deposit plane by plane
  for (unsigned i_plane = 0; i_plane < fNPlanes; i_plane++){

    // Calorimetric planes are saved in inversed order 
    int calo_plane = - i_plane + 2 ;

    //////////////////////////////////////////////////////
    // Save histograms

    // Save the dQ/dX 2D distributions
    canv_dQdX_dQ[i_plane]   = tfs->make<TCanvas>() ; 
    Save1Hist2d(canv_dQdX_dQ[i_plane], hist2d_dQdX_dQ[i_plane], Form("canv_dQdX_dQ_plane%d", calo_plane)) ;
    
    canv_dQdX_phi[i_plane]   = tfs->make<TCanvas>() ; 
    Save1Hist2d(canv_dQdX_phi[i_plane], hist2d_dQdX_phi[i_plane], Form("canv_dQdX_phi_plane%d", calo_plane)) ;
    
    canv_dQdX_theta[i_plane]   = tfs->make<TCanvas>() ; 
    Save1Hist2d(canv_dQdX_theta[i_plane], hist2d_dQdX_theta[i_plane], Form("canv_dQdX_theta_plane%d", calo_plane)) ;

    canv_y_z[i_plane]   = tfs->make<TCanvas>() ; 
    Save1Hist2d(canv_y_z[i_plane], hist2d_y_z[i_plane], Form("canv_y_z_plane%d", calo_plane)) ;

    canv_dQdX_yz[i_plane]   = tfs->make<TCanvas>() ; 
    Save1Hist2d(canv_dQdX_yz[i_plane], hist2d_dQdX_yz[i_plane], Form("canv_dQdX_yz_plane%d", calo_plane)) ;

    canv_dQdX_yx[i_plane]   = tfs->make<TCanvas>() ; 
    Save1Hist2d(canv_dQdX_yx[i_plane], hist2d_dQdX_yx[i_plane], Form("canv_dQdX_yx_plane%d", calo_plane)) ;

    canv_dQdX_zx[i_plane]   = tfs->make<TCanvas>() ; 
    Save1Hist2d(canv_dQdX_zx[i_plane], hist2d_dQdX_zx[i_plane], Form("canv_dQdX_zx_plane%d", calo_plane)) ;

    // Save the dQ/dX 1D histogram pre-fit
    canv_dQdX[i_plane]   = tfs->make<TCanvas>() ; 
    Save1Hist1d(canv_dQdX[i_plane], hist1d_dQdX[i_plane], Form("canv_dQdX_plane%d", calo_plane)) ;

    //////////////////////////////////////////////////////
    // Fit the electron lifetime
    /*
    if( fLogLevel >= 1 ) std::cout << "------ Starting electron lifetime fit for plane " << calo_plane << std::endl ;
    FitDqdxTime(hist2d_dQdX_time[i_plane], fit_dqdx_time[i_plane], calo_plane) ;

    // Fit the dQ/dX distribution
    if( fLogLevel >= 1 ) std::cout << "------ Starting dQ/dX fit with Landau-Gauss function for plane " << calo_plane << std::endl ;
    fit_dqdx[i_plane] = FitDqdx(hist1d_dQdX[i_plane], calo_plane, true) ; 

    canv_dQdX_fit[i_plane] = tfs->make<TCanvas>() ; 
    SaveDqdxFit(canv_dQdX_fit[i_plane], hist1d_dQdX[i_plane], fit_dqdx[i_plane], Form("canv_dQdX_fit_plane%d", calo_plane)) ;
    
    if( fLogLevel >= 1 ){
      std::cout << "--- dQ/dX parameters have been fitted to: " << std::endl ;
      for (Int_t i_par = 0 ; i_par < fit_dqdx[i_plane] -> GetNpar() ; i_par ++)      std::cout << " - " << fit_dqdx[i_plane] -> GetParName(i_par) << ": " << fit_dqdx[i_plane] -> GetParameter(i_par) << std::endl ;
    }

    ////////////////////////////////////////////////////////////////
    // Apply electron lifetime correction to dQ/dX distribution

    // Fill dQ/dX histograms with corrected distribution if electron lifetime has been measured
    if (e_lifetime[i_plane] > 0)
    {
      std::cout << "check " << fCaloHitXcorr_alltrks[0][0][0] << std::endl ;
       for(unsigned int i_trk=0 ; i_trk<fCaloHitDqdx_alltrks.size() ; i_trk++) {
            for(unsigned int i_hit=0 ; i_hit<fCaloHitDqdx_alltrks[i_trk][i_plane].size() ; i_hit++) {
              double e_corr = TMath::Exp(fCaloHitXcorr_alltrks[i_trk][i_plane][i_hit] / (e_lifetime[i_plane] * fCbvdDriftSpeed)) ;
      
              hist1d_dQdX_corr[i_plane] -> Fill(fCaloHitDqdx_alltrks[i_trk][i_plane][i_hit] * e_corr) ;
            }
        }

      // Perform dQ/dX fit    
      if( fLogLevel >= 1 ) std::cout << "------ Starting corrected dQ/dX fit with Landau-Gauss function for plane " << calo_plane << std::endl ;
      fit_dqdx_corr[i_plane] = FitDqdx(hist1d_dQdX_corr[i_plane], calo_plane, true) ; 

      canv_dQdXcorr_fit[i_plane] = tfs->make<TCanvas>() ; 
      SaveDqdxFit(canv_dQdXcorr_fit[i_plane], hist1d_dQdX_corr[i_plane], fit_dqdx_corr[i_plane], Form("canv_dQdXcorr_fit_plane%d", calo_plane)) ;
      
      if( fLogLevel >= 1 ){
        std::cout << "--- dQ/dX parameters have been fitted to: " << std::endl ;
        for (Int_t i_par = 0 ; i_par < fit_dqdx_corr[i_plane] -> GetNpar() ; i_par ++)      std::cout << " - " << fit_dqdx_corr[i_plane] -> GetParName(i_par) << ": " << fit_dqdx_corr[i_plane] -> GetParameter(i_par) << std::endl ;
      }
      }
  */
  
  } // end of loop over planes
}

TF1* pdvdana::ColdboxTrackStudy::FitDqdx(TH1F* hist1d_dqdx, int plane, bool save) { 

  // Create the function
  art::ServiceHandle<art::TFileService> tfs;
  TF1* fit_dqdx = tfs->make<TF1>(Form("fit_dqdx_lg_plane%d",-plane+2), LandauGaussFct, fFitDqdxXmin, fFitDqdxXmax, 5);
  fit_dqdx -> SetParNames("Landau width","Landau MPV","Integral","Gaussian Width", "Step");

  // Set the parameters initial values
  fit_dqdx -> SetParameters(1.,
          hist1d_dqdx -> GetBinCenter(hist1d_dqdx -> GetMaximumBin()),
          hist1d_dqdx -> Integral(),
          0.5,
          hist1d_dqdx -> GetBinContent(hist1d_dqdx -> FindBin(5.)));
  
  fit_dqdx -> SetParLimits(0, 0.                             , 3.0);
  fit_dqdx -> SetParLimits(1, 0.9 * fit_dqdx->GetParameter(1), 1.5 * fit_dqdx->GetParameter(1));
  fit_dqdx -> SetParLimits(2, 0.                             , 5.0 * fit_dqdx->GetParameter(2));
  fit_dqdx -> SetParLimits(3, 0.                             , 3) ; 
  fit_dqdx -> SetParLimits(4, 0.3 * fit_dqdx->GetParameter(4), 3.0 * fit_dqdx->GetParameter(4));

  if( fLogLevel >= 3 ){ 
    std::cout << "--- Initial parameters are set to: " << std::endl ;
    for (Int_t i_par = 0 ; i_par < fit_dqdx -> GetNpar() ; i_par ++)      std::cout << " - " << fit_dqdx -> GetParName(i_par) << ": " << fit_dqdx -> GetParameter(i_par) << std::endl ;    
  }

  // Perform fit and save if needed
  hist1d_dqdx -> Fit(Form("fit_dqdx_lg_plane%d",-plane+2), "RQ");

  if(save == true) fit_dqdx -> Write() ;

  return fit_dqdx; 
}

void pdvdana::ColdboxTrackStudy::FitDqdxTime(TH2F* hist2d_dQdX_time, TF1* fit_dqdx_time, int i_plane) {

  // Fit of dQ/dX MPV vs time 
  if( fLogLevel >= 4 ) std::cout << "--- " << hist2d_dQdX_time -> GetNbinsX() << " dQ/dX distributions will be fitted " <<std::endl ;
  
  // Create graph to store dQ/dX MPV vs drift time
  art::ServiceHandle<art::TFileService> tfs;
  TGraphErrors* graph_dQdX_time = tfs->make<TGraphErrors>();
  graph_dQdX_time -> SetName(Form("graph_dQdX_time_plane%d", i_plane));

  // Loop over drift time bin to get dQ/dX MPV
  for(int i_bin=1 ; i_bin < hist2d_dQdX_time -> GetNbinsX(); i_bin ++){
    
    // Project dQ/dX 1d histogram for this drift time bin
    TH1F* hist1d_dqdx_time = (TH1F*)hist2d_dQdX_time -> ProjectionY("hist1d_dqdx_time", i_bin, i_bin+1);
    if (hist1d_dqdx_time -> GetEntries() < 10){
      if( fLogLevel >= 4 ) std::cout << " - Skipping drift time bin " << i_bin << " due to lack of statistics. " << std::endl ;
      continue ; 
    }
    
    // Fit dQ/dX 1d histogram to get MPV in this drift time bin
    TF1* fit_dqdx_timebin = FitDqdx(hist1d_dqdx_time, i_plane, false);
    double fit_error = fit_dqdx_timebin->GetParError(1) / fit_dqdx_timebin->GetParameter(1) ;
    if (fit_error < 0.0001 || fit_error > 0.1) {
      if( fLogLevel >= 4 )	  std::cout << " - Skippping drift time bin " << i_bin << " due to not trusting dQ/dX MPV: " << fit_dqdx_timebin->GetParameter(1) << " +/- " << fit_dqdx_timebin->GetParError(1)  << std::endl ;
      continue ;   
    }
    
    // Store fitted dQ/dX MPV in TGraph
    graph_dQdX_time -> SetPoint(graph_dQdX_time -> GetN(), hist2d_dQdX_time -> GetXaxis() -> GetBinCenter(i_bin), fit_dqdx_timebin->GetParameter(1) ) ;
    graph_dQdX_time -> SetPointError(graph_dQdX_time -> GetN() - 1, 0, fit_dqdx_timebin->GetParError(1)) ;

  } // end of loop over time bins

  // Fit dQ/dX MPV vs time with an exponential function
  fit_dqdx_time = tfs->make<TF1>(Form("fit_dqdx_time_plane%d",i_plane), "[0]*TMath::Exp(-x*[1])", fFitDqdxTimeXmin, fFitDqdxTimeXmax) ;
  fit_dqdx_time -> SetParNames("Scale", "Tau") ;
  
  // Check that there is enough statistics
  if (graph_dQdX_time -> GetN() > 2) {
    if( fLogLevel >= 2 ) std::cout << "--- Fitting e- lifetime with : " << graph_dQdX_time -> GetN() << " points." << std::endl ;              

    // Initialise the function
    fit_dqdx_time -> SetParameter(0, graph_dQdX_time->GetY()[5]);
    fit_dqdx_time -> SetParameter(1, 0.) ;
    fit_dqdx_time -> SetParLimits(0, fFitDqdxTimeYmin, fFitDqdxTimeYmax) ;
    fit_dqdx_time -> SetParLimits(1, 1e-5,1e-1) ;

    if( fLogLevel >= 5 ){ 
      std::cout << "--- Initial parameters are set to: " << std::endl ;
      for (Int_t i_par = 0 ; i_par < fit_dqdx_time -> GetNpar() ; i_par ++)      std::cout << " - " << fit_dqdx_time -> GetParName(i_par) << ": " << fit_dqdx_time -> GetParameter(i_par) << std::endl ;    
    }

    // Perform the fit
    graph_dQdX_time -> Fit(Form("fit_dqdx_time_plane%d",i_plane), "R") ;
    
    if( fLogLevel >= 3 ){ 
      std::cout << "--- Parameters are fitted to: " << std::endl ;
      for (Int_t i_par = 0 ; i_par < fit_dqdx_time -> GetNpar() ; i_par ++)      std::cout << " - " << fit_dqdx_time -> GetParName(i_par) << ": " << fit_dqdx_time -> GetParameter(i_par) << std::endl ;    
    }

    e_lifetime.push_back(1. / fit_dqdx_time -> GetParameter(1)) ;
    e_lifetime_error.push_back(fit_dqdx_time -> GetParError(1) / pow(fit_dqdx_time -> GetParameter(1),2)) ;
    
    // Print result and save plot
    if( fLogLevel >= 1 ) std::cout << "--- Electron lifetime is found to be: " << 1. / fit_dqdx_time -> GetParameter(1) << " +/- " << fit_dqdx_time -> GetParError(1) / pow(fit_dqdx_time -> GetParameter(1),2) << " microseconds." << std::endl ;  				     

    canv_dQdX_time[i_plane] = tfs->make<TCanvas>() ; 
    SaveDqdxTimeFit(canv_dQdX_time[i_plane], hist2d_dQdX_time, graph_dQdX_time, fit_dqdx_time, Form("canv_dQdX_time_fit_plane%d", i_plane)) ;
  }
  else {
    if( fLogLevel >= 1 ) std::cout << "--- Not enough dQ/dX points to fit electron lifetime. Skipping." << std::endl ;
  } 

  return ;
}

void pdvdana::ColdboxTrackStudy::MakePrettyTH1F(TH1F* hist)
{
  hist -> SetLineColor(kCyan+2) ;
  hist -> SetLineWidth(2) ;
  hist -> SetFillColor(kCyan-8);
  hist -> GetXaxis() -> SetLabelSize(.05);
  hist -> GetXaxis() -> SetTitleSize(.05);
  hist -> GetYaxis() -> SetLabelSize(.05);
  hist -> GetYaxis() -> SetTitleSize(.05);
 
  return ;
}

void pdvdana::ColdboxTrackStudy::Save1Hist1d(TCanvas* canv, TH1F* hist1, TString filename)
{
  art::ServiceHandle<art::TFileService> tfs;
  gStyle -> SetCanvasColor(0);
  gStyle -> SetOptStat(111111);
  
  canv -> SetLeftMargin(0.15);
  canv -> SetBottomMargin(0.15);
  
  canv -> cd() ;
  
  MakePrettyTH1F(hist1) ; 
  hist1 -> Draw();
  canv -> Write(filename) ;
  canv -> SaveAs(filename+".pdf") ;
  //canv -> SaveAs(filename+".png") ;
  //canv -> SaveAs(filename+".C") ;
  //canv -> SaveAs(filename+".root") ;
  return ;
}

void pdvdana::ColdboxTrackStudy::SaveDqdxFit(TCanvas* canv, TH1F* hist1, TF1* fct, TString filename)
{
  gStyle -> SetCanvasColor(0);
  gStyle -> SetOptStat(111111);
  
  canv -> SetLeftMargin(0.15);
  canv -> SetBottomMargin(0.15);
  
  canv -> cd() ;
  MakePrettyTH1F(hist1) ; 
  hist1 -> Draw();
  //fct -> Draw("SAME") ;
  canv -> Write(filename) ;
  canv -> SaveAs(filename+".pdf") ;
  //canv -> SaveAs(filename+".png") ;
  //canv -> SaveAs(filename+".root") ;
  //canv -> SaveAs(filename+".C") ;
  
  return ;
}


void pdvdana::ColdboxTrackStudy::SaveDqdxTimeFit(TCanvas* canv, TH2F* hist, TGraphErrors* graph, TF1* fct, TString filename)
{
  gStyle -> SetCanvasColor(0);
  gStyle -> SetOptStat(0);
  gStyle -> SetPalette(kViridis);

  canv -> SetLeftMargin(0.15);
  canv -> SetBottomMargin(0.15);
  
  canv -> cd() ;
  hist -> Draw("COLZ");

  graph -> SetMarkerColor(kRed);
  graph -> SetMarkerColor(kDot);
  graph -> Draw("PSAME");
  
  fct -> SetLineStyle(7); 
  fct -> SetLineColor(kRed);
  fct -> Draw("SAME") ;

  canv -> Write(filename) ;
  canv -> SaveAs(filename+".pdf") ;
  //canv -> SaveAs(filename+".png") ;
  //canv -> SaveAs(filename+".root") ;
  //canv -> SaveAs(filename+".C") ;
  
  return ;  
}

void pdvdana::ColdboxTrackStudy::Save1Hist2d(TCanvas* canv, TH2F* hist2d, TString filename)
{
  gStyle -> SetCanvasColor(0);
  gStyle -> SetOptStat(0);
  gStyle -> SetPalette(kViridis);
  
  canv -> SetLeftMargin(0.15);
  canv -> SetBottomMargin(0.15);
  
  canv -> cd() ;
  hist2d -> Draw("colz");
  canv -> Write(filename) ;
  canv -> SaveAs(filename+".pdf") ;
  //canv -> SaveAs(filename+".png") ;
  canv -> SaveAs(filename+".C") ;
  //canv -> SaveAs(filename+".root") ;
  
  return ;
}

void pdvdana::ColdboxTrackStudy::Save3Hist1d(TCanvas* canv, TH1F* hist1, TH1F* hist2, TH1F* hist3, TString filename)
{
    gStyle -> SetCanvasColor(0);
    gStyle -> SetOptStat(111111);

    canv -> Divide(1,3);

    canv -> cd(1) ;
    MakePrettyTH1F(hist1) ; 
    hist1 -> Draw();

    canv -> cd(2) ;
    MakePrettyTH1F(hist2) ; 
    hist2 -> Draw();

    canv -> cd(3) ;
    MakePrettyTH1F(hist3) ; 
    hist3 -> Draw();
    
    canv -> Write(filename) ;
    canv -> SaveAs(filename+".pdf") ;
    //canv -> SaveAs(filename+".png") ;
    //canv -> SaveAs(filename+".C") ;
    //canv -> SaveAs(filename+".root") ;

    return ;
}


DEFINE_ART_MODULE(pdvdana::ColdboxTrackStudy)
