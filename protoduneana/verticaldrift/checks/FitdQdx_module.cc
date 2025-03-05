////////////////////////////////////////////////////////////////////////
// Class:       FitdQdx
// Plugin Type: analyzer (Unknown Unknown)
// File:        FitdQdx_module.cc
//
// Generated at Fri Feb 10 16:52:50 2023 by Yoann Kermaidic using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include <limits>  // std::numeric_limits<>

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

//LArSoft
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/ServicePack.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// ROOT
#include "Math/ProbFunc.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPaveStats.h"

using std::vector;
using std::string;

namespace pdvdana {
  class FitdQdx;
}

namespace {
  constexpr unsigned int int_max_as_unsigned_int{std::numeric_limits<int>::max()};
}


class pdvdana::FitdQdx : public art::EDAnalyzer {
public:
  explicit FitdQdx(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FitdQdx(FitdQdx const&) = delete;
  FitdQdx(FitdQdx&&) = delete;
  FitdQdx& operator=(FitdQdx const&) = delete;
  FitdQdx& operator=(FitdQdx&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  float    degTOrad = 3.14159/180.; // rad / deg
  float    fEltofC  = 1./1.60E-4;      // e- / fC
  float    fADCtoEl = 5E-3;         // ADC x tick / e-

private:

  int      fLogLevel;
  string   fTrackModuleLabel;

  // Cut parameters
  float    fTrackLenMin;
  float    fTrackLenMax;
  float    fYZfid;
  float    fQMin;

  // Retrieve electron speed in LAr
  float    fDriftSpeed;

  int      iDtBins;
  float    fDtMin;
  float    fDtMax;
  float    fDtFitMin;
  float    fDtFitMax;

  int      iDQdxBins;
  float    fDQdxMin;
  float    fDQdxMax;
  float    fDQdxFitMin;
  float    fDQdxFitMax;

  // Space charge related flags
  bool     bFieldDistortion;
  bool     bTrackIsFieldDistortionCorrected;
  bool     bFieldDistortionCorrectionXSign;

  // Detector properties
  unsigned int fNtpcs;
  unsigned int fNplanes;
  float    fHeight;
  float    fXmin = 1e6;
  float    fXmax =-1e6;
  float    fYmin = 1e6;
  float    fYmax =-1e6;
  float    fZmin = 1e6;
  float    fZmax =-1e6;

  vector<vector<vector<int> > > yxw;
  vector<vector<vector<int> > > zxw;
  vector<vector<vector<int> > > yzw;

  // detector geometry
  const geo::Geometry* fGeom;
  geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

  // Stored objects
  TTree*                fTree;
  vector<TH1F*>         h_dqdx_cor;
  vector<TH2F*>         h2_dqdx_yx;
  vector<TH2F*>         h2_dqdx_zx;
  vector<TH2F*>         h2_dqdx_yz;
  vector<TH2F*>         h2_dqdx_dt;
  vector<TH2F*>         h2_dqdx_t;
  vector<TH2F*>         h2_dqdx_p;
  vector<TCanvas*>      c_dqdx;
  vector<TCanvas*>      c_dqdx_z;
  vector<TF1*>          f_dqdx_z;
  vector<TGraphErrors*> g_dqdx_z;

  // Track-wise information
  int   fEventNum;
  int   fTrackId;
  float fTrackLength;
  float fTrackTheta;
  float fTrackPhi;
  float fTrackStartT;
  float fTrackStartX;
  float fTrackStartY;
  float fTrackStartZ;
  float fTrackEndT;
  float fTrackEndX;
  float fTrackEndY;
  float fTrackEndZ;
  // Plane-wise information
  vector<float>          fTrackQ;    // Total track charge in each view
  vector<int>            fNsnippets;
  vector<int>            fNvalid;
  vector<int>            fNhits;
  // Hit-wise information
  vector<vector<int> >   fCut;
  vector<vector<float> > fPosX;
  vector<vector<float> > fPosY;
  vector<vector<float> > fPosZ;
  vector<vector<float> > fQ;
  vector<vector<float> > fDqdx;

  bool   SetPhiFlag(float fPhi);
  bool   HitsInYXVolume(float fX, float fY, float fZ);
  bool   HitsInZXVolume(float fX, float fY, float fZ);
  bool   HitsInYZVolume(float fX, float fY, float fZ);
  bool   HitsInFiducialVolume(float fX, float fY, float fZ);
  void   GetAngles(const recob::Track& track, bool isRevert, float &theta, float &phi);
  TF1*   LGfit(TH1F* hin);

  //GnocchiCalorimetry module functions
  bool HitIsValid(const art::Ptr<recob::Hit> hit,
                  const recob::TrackHitMeta* thm,
                  const recob::Track& track);
  std::vector<std::vector<unsigned>> OrganizeHitsSnippets(
      const std::vector<art::Ptr<recob::Hit>>& hits,
      const std::vector<const recob::TrackHitMeta*>& thms,
      const recob::Track& track,
      unsigned nplanes,
      vector<int> &nhits);
  float GetPitch(const recob::Track& track,
                 const art::Ptr<recob::Hit> hit,
                 const recob::TrackHitMeta* meta);

  geo::Point_t GetLocationAtWires(const recob::Track& track,
                                    const art::Ptr<recob::Hit> hit,
                                    const recob::TrackHitMeta* meta);
  geo::Point_t WireToTrajectoryPosition(const geo::Point_t& loc, const geo::TPCID& tpc);
  geo::Point_t TrajectoryToWirePosition(const geo::Point_t& loc, const geo::TPCID& tpc);
};


pdvdana::FitdQdx::FitdQdx(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fLogLevel(                       p.get< int >("LogLevel")),
  fTrackModuleLabel(               p.get< std::string  >("TrackModuleLabel")),
  fTrackLenMin(                    p.get< float >("TrackLenMin")),
  fTrackLenMax(                    p.get< float >("TrackLenMax")),
  fYZfid(                          p.get< float >("YZfid")),
  fQMin(                           p.get< float >("QMin")),
  iDtBins(                         p.get< int   >("DtBins") ),
  fDtMin(                          p.get< float >("DtMin") ),
  fDtMax(                          p.get< float >("DtMax") ),
  fDtFitMin(                       p.get< float >("DtFitMin") ),
  fDtFitMax(                       p.get< float >("DtFitMax") ),
  iDQdxBins(                       p.get< int   >("DQdxBins") ),
  fDQdxMin(                        p.get< float >("DQdxMin") ),
  fDQdxMax(                        p.get< float >("DQdxMax") ),
  fDQdxFitMin(                     p.get< float >("DQdxFitMin") ),
  fDQdxFitMax(                     p.get< float >("DQdxFitMax") ),
  bFieldDistortion(                p.get< bool  >("FieldDistortion")),
  bTrackIsFieldDistortionCorrected(p.get< bool  >("TrackIsFieldDistortionCorrected")),
  bFieldDistortionCorrectionXSign( p.get< bool  >("FieldDistortionCorrectionXSign"))
{
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
}

void pdvdana::FitdQdx::analyze(art::Event const& e)
{
    fEventNum = e.id().event();

    if( fLogLevel >= 2 ) std::cout << "Start analysing event " << fEventNum << " ..." << std::endl;

    // Get list of reconstructed tracks
    auto Tracks   = e.getValidHandle<vector<recob::Track>>(fTrackModuleLabel);
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmHits(Tracks, e, fTrackModuleLabel);

    if( fLogLevel >= 2 ) std::cout << "  #tracks:    " << Tracks->size()    << std::endl;

    // Loop over tracks
    for (unsigned trk = 0; trk < Tracks->size(); ++trk) {

      const recob::Track& track = Tracks->at(trk);

      // Initialize output vectors
      fTrackLength= -999;
      fTrackTheta = -999;
      fTrackPhi   = -999;
      fTrackStartT= 10000;
      fTrackStartX= -999;
      fTrackStartY= -999;
      fTrackStartZ= -999;
      fTrackEndT  = -999;
      fTrackEndX  = -999;
      fTrackEndY  = -999;
      fTrackEndZ  = -999;

      for (unsigned plane_i = 0; plane_i < fNplanes; plane_i++){
        fNsnippets[plane_i]  = 0;
        fNvalid[plane_i]     = 0;
        fTrackQ[plane_i]     = -999;
        fPosX[plane_i].resize(1,-999);
        fPosY[plane_i].resize(1,-999);
        fPosZ[plane_i].resize(1,-999);
        fDqdx[plane_i].resize(1,-999);
        fCut[plane_i].resize(1,-999);
        fQ[plane_i].resize(1,-999);
      }

      // Cut tracks that are not within specified thresholds
      if(track.Length() < fTrackLenMin || track.Length() > fTrackLenMax) continue;

      // Retrieve begining and end of the track
      fTrackStartX = track.Start().X();
      fTrackStartY = track.Start().Y();
      fTrackStartZ = track.Start().Z();
      fTrackEndX   = track.End().X();
      fTrackEndY   = track.End().Y();
      fTrackEndZ   = track.End().Z();

      // Remove tracks that start outside the detector in the horizontal plane
      if(fTrackStartY < fYmin || fTrackStartY > fYmax) continue;
      if(fTrackStartZ < fZmin || fTrackStartZ > fZmax) continue;

      // Remove tracks that end   outside the detector in the horizontal plane
      if(fTrackEndY < fYmin   || fTrackEndY > fYmax)   continue;
      if(fTrackEndZ < fZmin   || fTrackEndZ > fZmax)   continue;

      // Remove tracks that are not crossing the anode/cathode planes
      // -> unknown drift time
      if(fabs(fTrackStartX - fTrackEndX) < 0.95*fHeight || fabs(fTrackStartX - fTrackEndX) > 1.05*fHeight) continue;

      // Check if the track is reconstructed upwards or downwards
      // for now, assume all reconstructed muons should be going downwards
      bool isRevert = false;
      if(fTrackEndX - fTrackStartX > 0) isRevert = true;

      fTrackLength = track.Length();
      GetAngles(track,isRevert,fTrackTheta,fTrackPhi);

      if( fLogLevel >= 3 ) std::cout << "    track " << trk << " -> length:    " << track.Length() << " - theta: " << fTrackTheta << " - phi: " << fTrackPhi << std::endl;
      if( fLogLevel >= 4 ) std::cout << "    track first/last valid point: " <<  track.FirstValidPoint() << "/" << track.LastValidPoint() << std::endl;

      // if the track passed selection criteria, loop over hits to retrieve pulse width
      const std::vector<art::Ptr<recob::Hit>>& Hits       = fmHits.at(trk);
      const std::vector<const recob::TrackHitMeta*>& thms = fmHits.data(trk);

      // follow GnocchiCalorimatery module to sort hits by plane and snippets
      std::vector<std::vector<unsigned>> hit_indices = OrganizeHitsSnippets(Hits, thms, track, fNplanes, fNhits);


      // Loop over planes
      for (unsigned plane_i = 0; plane_i < fNplanes; plane_i++){

        if( fLogLevel >= 3 )
          std::cout << "    #hits in plane " << plane_i << ": " << hit_indices[plane_i].size() << std::endl;

        unsigned Nhits = hit_indices[plane_i].size();

        fNsnippets[plane_i] = Nhits;
        fTrackQ[plane_i]    = 0;
        fPosX[plane_i].clear();    fPosX[plane_i].resize(Nhits,-999);
        fPosY[plane_i].clear();    fPosY[plane_i].resize(Nhits,-999);
        fPosZ[plane_i].clear();    fPosZ[plane_i].resize(Nhits,-999);
        fDqdx[plane_i].clear();    fDqdx[plane_i].resize(Nhits,-999);
        fCut[plane_i].clear();     fCut[plane_i].resize(Nhits,-999);
        fQ[plane_i].clear();       fQ[plane_i].resize(Nhits,-999);

        // Loop over hits
        for (unsigned hit_i = 0; hit_i < Nhits; hit_i++) {
          if( fLogLevel >= 4 ) std::cout << "      hit " << hit_i << std::endl;

          unsigned hit_index = hit_indices[plane_i][hit_i];
          unsigned tms_index = thms[hit_index]->Index();

          if( fLogLevel >= 4 ) std::cout << "        -> indexes: " << hit_index << "/" << tms_index << "-" << tms_index << " -> time: " << Hits[hit_index]->PeakTime()  << " start/end: " << fTrackStartT << "/" << fTrackEndT << std::endl;

          //////////////////
          // CUT BLOCK
          //

          // Remove hits associated to the track that are not valid points
          if(!track.HasValidPoint(tms_index)){if( fLogLevel >= 4 ) std::cout << "        #1 nok" << std::endl;      continue;}
          if( fLogLevel >= 6 ) std::cout << "        #1 ok" << std::endl;

          // Remove first and last valid point that are most probably
          // starting/ending in the anode/cathode/field cage
          if(tms_index == track.FirstValidPoint()){if( fLogLevel >= 4 ) std::cout << "        #2 nok" << std::endl; continue;}
          if( fLogLevel >= 6 ) std::cout << "        #2 ok" << std::endl;
          if(tms_index == track.LastValidPoint()){ if( fLogLevel >= 4 ) std::cout << "        #3 nok" << std::endl; continue;}
          if( fLogLevel >= 6 ) std::cout << "        #3 ok" << std::endl;

          ////////////////

          float pitch  = GetPitch(track, Hits[hit_index], thms[hit_index]);
          float charge = Hits[hit_index]->ROISummedADC() / fADCtoEl / fEltofC;

          // Retrieve (x,y,z) point coordinates
          float fX = track.LocationAtPoint(tms_index).X();
          float fY = track.LocationAtPoint(tms_index).Y();
          float fZ = track.LocationAtPoint(tms_index).Z();

          // Shift the vertical hits coordinate assuming all muons are going downwards
          // due to prefered reconstructed direction along the beam axis (Z coordinate) in pandora
          if(isRevert) fX = fHeight + fTrackStartX - fX;
          else         fX = fTrackStartX - fX;

          if( fLogLevel >= 5 )
            std::cout << "        -> ( " << fX << "," << fY << "," << fZ << ") - pitch: " << pitch << " - charge:" << charge << " - dQdx: " << charge/pitch << std::endl;

          // Store quantities of interest into vector -> output file
          fNvalid[plane_i]++;
          fPosX[plane_i][hit_i]    = fX;
          fPosY[plane_i][hit_i]    = fY;
          fPosZ[plane_i][hit_i]    = fZ;
          fTrackQ[plane_i]        += charge;
          fQ[plane_i][hit_i]       = charge;
          fDqdx[plane_i][hit_i]    = charge/pitch;
          fCut[plane_i][hit_i]     = 2;

          if(fTrackStartT > Hits[hit_index]->PeakTime()) fTrackStartT = Hits[hit_index]->PeakTime();
          if(fTrackEndT   < Hits[hit_index]->PeakTime()) fTrackEndT   = Hits[hit_index]->PeakTime();

          if( fLogLevel >= 6 ) std::cout << "        vectors set... " << std::endl;

          // Remove hits from noise pedestal
          if(charge < fQMin){if( fLogLevel >= 4 ) std::cout << "        #4 nok" << std::endl;                       continue;}
          if( fLogLevel >= 6 ) std::cout << "        #4 ok" << std::endl;
          fCut[plane_i][hit_i]  = 1;

          // Remove hits that reconstructed outside the fiducial volume
          if(!HitsInFiducialVolume(fX,fY,fZ)){if( fLogLevel >= 4 ) std::cout << "        #5 nok" << std::endl;      continue;}
          if( fLogLevel >= 6 ) std::cout << "        #5 ok" << std::endl;

          // Fill 2D histogram later used for drift time-wise fitting of dQ/dx distributions
          h2_dqdx_dt[plane_i]->Fill(fX/fDriftSpeed,charge/pitch);

          // Not cut by any criteria (enters dQ/dx vs drift time fit)
          fCut[plane_i][hit_i]  = 0;

        } // end of hit loop
      } // end of plane loop

      fTree->Fill();
      fTrackId++;

    } // end of track loop
} // end of event loop

bool pdvdana::FitdQdx::SetPhiFlag(float phi){

  // along plane 0 strip
  if(phi > 145  && phi < 155)  return true;
  if(phi > -35  && phi < -25)  return true;
  // along plane 1 strip
  if(phi > -155 && phi < -145) return true;
  if(phi > 25   && phi < 35)   return true;
  // along plane 2 strip
  if(phi > -95  && phi < -85)  return true;
  if(phi > 85   && phi < 95 )  return true;

  return false;
}

bool pdvdana::FitdQdx::HitsInYXVolume(float fX, float fY, float fZ){

  if(fX <= 0            || fX >= fXmax-fXmin)  return false;
  if(fY <= fYmin        || fY >= fYmax)        return false;
  if(fZ <= fZmin+fYZfid || fZ >= fZmax-fYZfid) return false;

  return true;
}
bool pdvdana::FitdQdx::HitsInZXVolume(float fX, float fY, float fZ){

  if(fX <= 0            || fX >= fXmax-fXmin)  return false;
  if(fY <= fYmin+fYZfid || fY >= fYmax-fYZfid) return false;
  if(fZ <= fZmin        || fZ >= fZmax)        return false;

  return true;
}
bool pdvdana::FitdQdx::HitsInYZVolume(float fX, float fY, float fZ){

  if(fX <= 0     || fX >= fXmax-fXmin)  return false;
  if(fY <= fYmin || fY >= fYmax)        return false;
  if(fZ <= fZmin || fZ >= fZmax)        return false;

  return true;
}
bool pdvdana::FitdQdx::HitsInFiducialVolume(float fX, float fY, float fZ){

  if(fX <= 0            || fX >= fXmax-fXmin)  return false;
  if(fY <= fYmin+fYZfid || fY >= fYmax-fYZfid) return false;
  if(fZ <= fZmin+fYZfid || fZ >= fZmax-fYZfid) return false;

  return true;
}

void pdvdana::FitdQdx::beginJob()
{
  fNtpcs   = fGeom->NTPC();
  fNplanes = fWireReadoutGeom.Nplanes();

  if( fLogLevel >= 1 ){
    std::cout << "  #TPCs:       "  << fNtpcs   << std::endl;
    std::cout << "  #planes:     "  << fNplanes << std::endl;
    for(unsigned i=0;i<fNtpcs;i++){
      geo::TPCID tpcid{0, i};
      std::cout << "  TPC " << i << " center: ("<< fGeom->TPC(tpcid).GetCenter().X()      << "," << fGeom->TPC(tpcid).GetCenter().Y()      << ","<< fGeom->TPC(tpcid).GetCenter().Z() << ")"
                                 << " box:  ["  << fGeom->TPC(tpcid).BoundingBox().MinX() << "," << fGeom->TPC(tpcid).BoundingBox().MaxX() << "]"
                                        << "["  << fGeom->TPC(tpcid).BoundingBox().MinY() << "," << fGeom->TPC(tpcid).BoundingBox().MaxY() << "]"
                                        << "["  << fGeom->TPC(tpcid).BoundingBox().MinZ() << "," << fGeom->TPC(tpcid).BoundingBox().MaxZ() << "]" << std::endl;

      if(fXmin > fGeom->TPC(tpcid).BoundingBox().MinX()) fXmin = fGeom->TPC(tpcid).BoundingBox().MinX();
      if(fXmax < fGeom->TPC(tpcid).BoundingBox().MaxX()) fXmax = fGeom->TPC(tpcid).BoundingBox().MaxX();
      if(fYmin > fGeom->TPC(tpcid).BoundingBox().MinY()) fYmin = fGeom->TPC(tpcid).BoundingBox().MinY();
      if(fYmax < fGeom->TPC(tpcid).BoundingBox().MaxY()) fYmax = fGeom->TPC(tpcid).BoundingBox().MaxY();
      if(fZmin > fGeom->TPC(tpcid).BoundingBox().MinZ()) fZmin = fGeom->TPC(tpcid).BoundingBox().MinZ();
      if(fZmax < fGeom->TPC(tpcid).BoundingBox().MaxZ()) fZmax = fGeom->TPC(tpcid).BoundingBox().MaxZ();
    }
    for(unsigned vw=0;vw<fNplanes;vw++){
      std::cout << "  pitch in view: " << vw << " ";
      if(vw == 0)      std::cout << fWireReadoutGeom.Plane({0, 0, geo::kU}).WirePitch() << std::endl;
      else if(vw == 1) std::cout << fWireReadoutGeom.Plane({0, 0, geo::kV}).WirePitch() << std::endl;
      else if(vw == 2) std::cout << fWireReadoutGeom.Plane({0, 0, geo::kZ}).WirePitch() << std::endl;
    }
  }
  fHeight = fXmax-fXmin;

  std::cout << "  Geometry boundaries: [" << fXmin << "," << fXmax << "] [" << fYmin << "," << fYmax << "] ["<< fZmin << "," << fZmax << "]" << std::endl;

  int iXbins = 5*int(fXmax-fXmin)+1;
  int iYbins =   int(fYmax-fYmin)+1;
  int iZbins =   int(fZmax-fZmin)+1;

  std::cout << "  Geometry binning: (" << iXbins << "," << iYbins << "," << iZbins << ")" << std::endl;

  auto const detProp   = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  fDriftSpeed = detProp.DriftVelocity();

  std::cout << "  Detector drift speed: " << fDriftSpeed << " cm/us" << std::endl;

  // Set size of output vectors (# of planes)
  c_dqdx.resize(fNplanes);
  c_dqdx_z.resize(fNplanes);
  yxw.resize(fNplanes);
  zxw.resize(fNplanes);
  yzw.resize(fNplanes);
  h_dqdx_cor.resize(fNplanes);
  h2_dqdx_yx.resize(fNplanes);
  h2_dqdx_zx.resize(fNplanes);
  h2_dqdx_yz.resize(fNplanes);
  h2_dqdx_dt.resize(fNplanes);
  h2_dqdx_t.resize(fNplanes);
  h2_dqdx_p.resize(fNplanes);
  f_dqdx_z.resize(fNplanes);
  g_dqdx_z.resize(fNplanes);
  fTrackQ.resize(fNplanes);
  fNsnippets.resize(fNplanes);
  fNvalid.resize(fNplanes);
  fPosX.resize(fNplanes);
  fPosY.resize(fNplanes);
  fPosZ.resize(fNplanes);
  fDqdx.resize(fNplanes);
  fCut.resize(fNplanes);
  fQ.resize(fNplanes);

  // Attach TH2 to output file
  art::ServiceHandle<art::TFileService> tfs;
  for (unsigned plane_i = 0; plane_i < fNplanes; plane_i++){
    yxw[plane_i].resize(iYbins);
    zxw[plane_i].resize(iZbins);
    yzw[plane_i].resize(iYbins);
    for(int y_i = 0; y_i < iYbins; y_i++){yzw[plane_i][y_i].resize(iZbins,0); yxw[plane_i][y_i].resize(iXbins,0);}
    for(int z_i = 0; z_i < iZbins; z_i++) zxw[plane_i][z_i].resize(iXbins,0);

    h2_dqdx_yz[plane_i] = tfs->make<TH2F>(Form("h2_dqdx_yz_%d",plane_i),Form("Plane %d;Y (cm);Z (cm);<dQ/dx> (fC/cm)",plane_i),iYbins,fYmin,fYmax,iZbins,fZmin,fZmax);
    h2_dqdx_yx[plane_i] = tfs->make<TH2F>(Form("h2_dqdx_yx_%d",plane_i),Form("Plane %d;Y (cm);X (cm);<dQ/dx> (fC/cm)",plane_i),iYbins,fYmin,fYmax,iXbins,0,fXmax-fXmin);
    h2_dqdx_zx[plane_i] = tfs->make<TH2F>(Form("h2_dqdx_zx_%d",plane_i),Form("Plane %d;Z (cm);X (cm);<dQ/dx> (fC/cm)",plane_i),iZbins,fZmin,fZmax,iXbins,0,fXmax-fXmin);
    h2_dqdx_dt[plane_i] = tfs->make<TH2F>(Form("h2_dqdx_dt_%d",plane_i),Form("Plane %d;Drift time (#mus);dQ/dx (fC/cm)",plane_i),iDtBins,fDtMin,fDtMax,iDQdxBins,fDQdxMin,fDQdxMax);
    h2_dqdx_t[plane_i]  = tfs->make<TH2F>(Form("h2_dqdx_t_%d",plane_i), Form("Plane %d;Theta (deg);dQ/dx (fC/cm)",plane_i),360,180.,180.,iDQdxBins,fDQdxMin,fDQdxMax);
    h2_dqdx_p[plane_i]  = tfs->make<TH2F>(Form("h2_dqdx_p_%d",plane_i), Form("Plane %d;Phi (deg);dQ/dx (fC/cm)",plane_i),360,-180.,180.,iDQdxBins,fDQdxMin,fDQdxMax);
    h_dqdx_cor[plane_i] = tfs->make<TH1F>(Form("h_dqdx_cor_%d",plane_i),Form("Plane %d;dQ/dx (fC/cm); Counts (/%.1f fC/cm)",plane_i,(fDQdxMax-iDQdxBins)/iDQdxBins),iDQdxBins,fDQdxMin,fDQdxMax);
  }

  fTree = tfs->make<TTree>("fitdqdxTree","Store dqdx info" );
  fTree->Branch("EventNum",   &fEventNum,   "EventNum/i"   );
  fTree->Branch("TrackId",    &fTrackId,    "TrackId/i"    );
  fTree->Branch("TrackLen",   &fTrackLength,"TrackLen/F"   );
  fTree->Branch("TrackTheta", &fTrackTheta, "TrackTheta/F" );
  fTree->Branch("TrackPhi",   &fTrackPhi,   "TrackPhi/F"   );
  fTree->Branch("TrackStartT",&fTrackStartT,"TrackStartT/F");
  fTree->Branch("TrackStartX",&fTrackStartX,"TrackStartX/F");
  fTree->Branch("TrackStartY",&fTrackStartY,"TrackStartY/F");
  fTree->Branch("TrackStartZ",&fTrackStartZ,"TrackStartZ/F");
  fTree->Branch("TrackEndT",  &fTrackEndT,  "TrackEndT/F"  );
  fTree->Branch("TrackEndX",  &fTrackEndX,  "TrackEndX/F"  );
  fTree->Branch("TrackEndY",  &fTrackEndY,  "TrackEndY/F"  );
  fTree->Branch("TrackEndZ",  &fTrackEndZ,  "TrackEndZ/F"  );
  fTree->Branch("TrackQ",     &fTrackQ   );
  fTree->Branch("nSnippets",  &fNsnippets);
  fTree->Branch("nHits",      &fNhits    );
  fTree->Branch("nValid",     &fNvalid   );
  fTree->Branch("HitCut",     &fCut      );
  fTree->Branch("HitX",       &fPosX     );
  fTree->Branch("HitY",       &fPosY     );
  fTree->Branch("HitZ",       &fPosZ     );
  fTree->Branch("HitQ",       &fQ        );
  fTree->Branch("HitDqdx",    &fDqdx     );
}

double langaufun(double *x, double *par) {
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

TF1* pdvdana::FitdQdx::LGfit(TH1F* hin){

  TF1* fit = new TF1("fit",langaufun,5,20,5);
  fit->SetParameters(1.,
                     hin->GetBinCenter(hin->GetMaximumBin()),
                     hin->Integral(),
                     0.5,
                     hin->GetBinContent(hin->FindBin(5.)));

  fit->SetParLimits(0,0.                      ,3.0);
  fit->SetParLimits(1,0.9*fit->GetParameter(1),1.5*fit->GetParameter(1));
  fit->SetParLimits(2,0.                      ,5.0*fit->GetParameter(2));
  fit->SetParLimits(3,0.                      ,3.0);
  fit->SetParLimits(4,0.3*fit->GetParameter(4),3.0*fit->GetParameter(4));

  fit->SetParNames("Width","MP","Area","GSigma","Step");

  hin->Fit("fit","QR");

  return fit;
}

void pdvdana::FitdQdx::endJob()
{
  if( fLogLevel >= 1 ) std::cout << "Start fitting" << std::endl;

  // Setup TF1 and TGraph to be stored in the output file
  art::ServiceHandle<art::TFileService> tfs;
  for (unsigned plane_i = 0; plane_i < fNplanes; plane_i++){
    c_dqdx_z[plane_i] = tfs->make<TCanvas>(Form("c_dqdx_z_%d",plane_i),Form("dQdx versus dist to anode - plane %d",plane_i),700,500);
    c_dqdx[plane_i]   = tfs->make<TCanvas>(Form("c_dqdx_%d",plane_i),  Form("dQdx distribution - plane %d",plane_i),700,500);
    f_dqdx_z[plane_i] = tfs->make<TF1>(Form("f_dqdx_z_%d",plane_i),"[0]*TMath::Exp(-x*[1])",fDtFitMin,fDtFitMax);
    g_dqdx_z[plane_i] = tfs->make<TGraphErrors>();
    g_dqdx_z[plane_i]->SetName(Form("g_dqdx_z_%d",plane_i));
    g_dqdx_z[plane_i]->SetMarkerStyle(24);
  }

  vector<float> tau(fNplanes,0);
  vector<float> dtau(fNplanes,0);

  for (unsigned plane_i = 0; plane_i < fNplanes; plane_i++){
    // TF1 for fitting dQ/dx in each individual TH2 bin

    for(int bin=1;bin<h2_dqdx_dt[fNplanes-1]->GetNbinsX();bin++){
      if( fLogLevel >= 2 ) std::cout << "Fit window " << bin << "/" << h2_dqdx_dt[plane_i]->GetNbinsX() << std::endl;

      // Retrieve drift time-wise dQ/dx distribution from TH2
      TH1F* hFit = (TH1F*)h2_dqdx_dt[plane_i]->ProjectionY("",bin,bin+1);

      if(hFit->GetEntries() < 10){
        if( fLogLevel >= 1 ) std::cout << "  WARNING: Not enough entries in bin " << bin << ". Skipped." << std::endl;
        continue;
      }

      // Perform gaus convolved landau fit of each dQ/dx distribution
      TF1* f_dqdx = LGfit(hFit);

      // Remove points with unexpectedly too small or too large uncertainty from the final fit
      if(f_dqdx->GetParError(1)/f_dqdx->GetParameter(1) < 0.05/100. || f_dqdx->GetParError(1)/f_dqdx->GetParameter(1) > 5./100.) continue;

      // Store individual dQ/dx mean values into a TGraph
      g_dqdx_z[plane_i]->SetPoint(g_dqdx_z[plane_i]->GetN(),h2_dqdx_dt[fNplanes-1]->GetXaxis()->GetBinCenter(bin),f_dqdx->GetParameter(1));
      g_dqdx_z[plane_i]->SetPointError(g_dqdx_z[plane_i]->GetN()-1,0,f_dqdx->GetParError(1));

    } // end of bins loop

    // Fit the binned dQ/dx as a function of drift time

    if(g_dqdx_z[plane_i]->GetN()>2){
      f_dqdx_z[plane_i]->SetLineStyle(7); f_dqdx_z[plane_i]->SetLineColor(kRed);
      f_dqdx_z[plane_i]->SetParameter(0,g_dqdx_z[plane_i]->GetY()[5]);
      f_dqdx_z[plane_i]->SetParameter(1,0.01);
      f_dqdx_z[plane_i]->SetParLimits(0,fDQdxFitMin,fDQdxFitMax);
      f_dqdx_z[plane_i]->SetParLimits(1,1e-4,1e-2);
      g_dqdx_z[plane_i]->Fit(f_dqdx_z[plane_i],"R");

      tau[plane_i]  = 1./(f_dqdx_z[plane_i]->GetParameter(1));
      dtau[plane_i] = f_dqdx_z[plane_i]->GetParError(1)/pow(f_dqdx_z[plane_i]->GetParameter(1),2);

      if( fLogLevel >= 1 ) std::cout << "  Electron lifetime: " << tau[plane_i] << " +/- " << dtau[plane_i] << " us" << std::endl;
    }
    else if( fLogLevel >= 1 ) std::cout << "  WARNING: Not enough data point for e- lifetime fit" << std::endl;
  } // end of plane loop

  //////////////////////////////
  // e- lifetime corrected dQ/dx

  for(int trk = 0;trk < fTree->GetEntries(); trk++){
    fTree->GetEntry(trk);

    // Check if the track has been recontructed along one of the strips
    bool isTrkAlongStrip = SetPhiFlag(fTrackPhi);

    for(unsigned int plane_i = 0;plane_i<fNplanes;plane_i++){
        for(unsigned int hit_i = 0; hit_i < fCut[plane_i].size(); hit_i++){

          if(fCut[plane_i][hit_i] < 0) continue;

          // Correct bare dQ/dx for e- lifetime
          float dqdx = fDqdx[plane_i][hit_i] * TMath::Exp(fPosX[plane_i][hit_i]/(tau[plane_i]*fDriftSpeed));

          if(fCut[plane_i][hit_i] == 0){
            // Fill 1D dQ/dx distributions
            h_dqdx_cor[plane_i]->Fill(dqdx);
            // Fill angular dependance 2D histograms
            h2_dqdx_t[plane_i] ->Fill(fTrackTheta,dqdx);
            h2_dqdx_p[plane_i] ->Fill(fTrackPhi,  dqdx);
          }

          // Fill 2D histograms used for position dependent checks
          // -> remove tracks along one of the strips to avoid <dQ/dx> overshoot
          if(!isTrkAlongStrip && fCut[plane_i][hit_i] <= 1){
            // YX 2D map
            if(HitsInYXVolume(fPosX[plane_i][hit_i],fPosY[plane_i][hit_i],fPosZ[plane_i][hit_i])){
              yxw[plane_i][h2_dqdx_yx[plane_i]->GetXaxis()->FindBin(fPosY[plane_i][hit_i])-1][h2_dqdx_yx[plane_i]->GetYaxis()->FindBin(fPosX[plane_i][hit_i])-1]++;
              h2_dqdx_yx[plane_i]->Fill(fPosY[plane_i][hit_i],fPosX[plane_i][hit_i],dqdx);
            }
            // ZX 2D map
            if(HitsInZXVolume(fPosX[plane_i][hit_i],fPosY[plane_i][hit_i],fPosZ[plane_i][hit_i])){
              zxw[plane_i][h2_dqdx_zx[plane_i]->GetXaxis()->FindBin(fPosZ[plane_i][hit_i])-1][h2_dqdx_zx[plane_i]->GetYaxis()->FindBin(fPosX[plane_i][hit_i])-1]++;
              h2_dqdx_zx[plane_i]->Fill(fPosZ[plane_i][hit_i],fPosX[plane_i][hit_i],dqdx);
            }
            // YZ 2D map
            if(HitsInYZVolume(fPosX[plane_i][hit_i],fPosY[plane_i][hit_i],fPosZ[plane_i][hit_i])){
              yzw[plane_i][h2_dqdx_yz[plane_i]->GetXaxis()->FindBin(fPosY[plane_i][hit_i])-1][h2_dqdx_yz[plane_i]->GetYaxis()->FindBin(fPosZ[plane_i][hit_i])-1]++;
              h2_dqdx_yz[plane_i]->Fill(fPosY[plane_i][hit_i],fPosZ[plane_i][hit_i],dqdx);
            }
          }
        } // end of hits loop
    } // end of planes loop
  } // end of tracks loop

  //////////////////////////////

  vector<float> mDqdx(fNplanes,0);
  vector<float> dmDqdx(fNplanes,0);
  vector<float> wDqdx(fNplanes,0);
  vector<float> dwDqdx(fNplanes,0);

  // Perform gaus convolved landau + step fit of the global dQ/dx distribution
  for(unsigned int plane_i = 0;plane_i<fNplanes;plane_i++){
    TF1* f_dqdx = LGfit(h_dqdx_cor[plane_i]);

    mDqdx[plane_i]   = f_dqdx->GetParameter(1);
    dmDqdx[plane_i]  = f_dqdx->GetParError(1);
    wDqdx[plane_i]   = f_dqdx->GetParameter(0);
    dwDqdx[plane_i]  = f_dqdx->GetParError(0);

    if( fLogLevel >= 1 ){
      std::cout << "  Plane " << plane_i << ":" << std::endl;
      std::cout << "    Mean dQ/dx = " << f_dqdx->GetParameter(1) << " +/- " << f_dqdx->GetParError(1) << " fC/cm" << std::endl;
      std::cout << "    Width      = " << f_dqdx->GetParameter(0) << " +/- " << f_dqdx->GetParError(0) << " fC/cm" << std::endl;
    }
  }

  // Plot histo into TCanvas for nice automatic rendering
  vector<vector<TLatex*> > lpars_1D(fNplanes,vector<TLatex*>(2));
  vector<TPaveStats*>      st_1D(fNplanes);
  vector<TList*>           list_1D(fNplanes);

  vector<vector<TLatex*> > lpars(fNplanes,vector<TLatex*>(2));
  vector<TPaveStats*>      st(fNplanes);
  vector<TList*>           list(fNplanes);

  for(unsigned int plane_i = 0;plane_i<fNplanes;plane_i++){
    c_dqdx[plane_i]->cd();
    h_dqdx_cor[plane_i]->GetXaxis()->SetNdivisions(409);
    h_dqdx_cor[plane_i]->GetYaxis()->SetNdivisions(409);
    h_dqdx_cor[plane_i]->SetLineColor(kBlack);
    h_dqdx_cor[plane_i]->SetStats("e");

    c_dqdx[plane_i]->Update();
    h_dqdx_cor[plane_i]->Draw("hist");
    c_dqdx[plane_i]->Update();

    st_1D[plane_i] = (TPaveStats*)h_dqdx_cor[plane_i]->GetListOfFunctions()->FindObject("stats");
    st_1D[plane_i]->SetOptStat(000000100);
    st_1D[plane_i]->SetOptFit(1);
    gPad->Modified(); gPad->Update();
    st_1D[plane_i]->SetX1NDC(0.45);
    st_1D[plane_i]->SetX2NDC(0.89);
    st_1D[plane_i]->SetY1NDC(0.7);
    st_1D[plane_i]->SetY2NDC(0.89);

    st_1D[plane_i]->SetName(Form("mystats_1D_%d",plane_i));

    c_dqdx[plane_i]->Update();
    h_dqdx_cor[plane_i]->SetStats(0);
    gPad->Modified(); gPad->Update();
    c_dqdx[plane_i]->Update();
/*
    list_1D[plane_i] = st_1D[plane_i]->GetListOfLines();
    lpars_1D[plane_i][0] = new TLatex(0, 0, Form("<dQ/dx_{0}> = %.2f +/- %.2f (fC/cm)",mDqdx[plane_i],mDqdx[plane_i]));
    lpars_1D[plane_i][1] = new TLatex(0, 0, Form("Width       = %.2f +/- %.2f (fC/cm)",wDqdx[plane_i],wDqdx[plane_i]));
    for(unsigned int par=0;par<2;par++){
      lpars_1D[plane_i][par]->SetTextFont(42);
      lpars_1D[plane_i][par]->SetTextSize(0.04);
      list_1D[plane_i]->Add(lpars_1D[plane_i][par]);
    }
*/

    c_dqdx_z[plane_i]->cd();
    h2_dqdx_dt[plane_i]->SetMinimum(1); //h2_dqdx_dt[plane_i]->SetMaximum(600);
    h2_dqdx_dt[plane_i]->GetXaxis()->SetNdivisions(409);
    h2_dqdx_dt[plane_i]->GetYaxis()->SetNdivisions(409);
    h2_dqdx_dt[plane_i]->SetStats("e");

    c_dqdx_z[plane_i]->Update();
    h2_dqdx_dt[plane_i]->Draw("colz");
    c_dqdx_z[plane_i]->Update();

    g_dqdx_z[plane_i]->Draw("P same");

    st[plane_i] = (TPaveStats*)h2_dqdx_dt[plane_i]->GetListOfFunctions()->FindObject("stats");
    st[plane_i]->SetOptStat(000000100);
    st[plane_i]->SetOptFit(1);
    gPad->Modified(); gPad->Update();
    st[plane_i]->SetX1NDC(0.45);
    st[plane_i]->SetX2NDC(0.89);
    st[plane_i]->SetY1NDC(0.7);
    st[plane_i]->SetY2NDC(0.89);

    st[plane_i]->SetName(Form("mystats_%d",plane_i));

    c_dqdx_z[plane_i]->Update();
    h2_dqdx_dt[plane_i]->SetStats(0);
    gPad->Modified(); gPad->Update();
    c_dqdx_z[plane_i]->Update();

    list[plane_i] = st[plane_i]->GetListOfLines();
    lpars[plane_i][0] = new TLatex(0, 0, Form("dQ/dx_{0} = %.3f +/- %.3f (fC/cm)",f_dqdx_z[plane_i]->GetParameter(0),f_dqdx_z[plane_i]->GetParError(0)));
    lpars[plane_i][1] = new TLatex(0, 0, Form("#tau     = %.1f +/- %.1f (#mus)",tau[plane_i],dtau[plane_i]));
    for(unsigned int par=0;par<2;par++){
      lpars[plane_i][par]->SetTextFont(42);
      lpars[plane_i][par]->SetTextSize(0.04);
      list[plane_i]->Add(lpars[plane_i][par]);
    }
    gPad->Modified(); gPad->Update();
  }

  // Write objects to output files
  for(unsigned plane_i = 0;plane_i<fNplanes;plane_i++){
    // Bin-wise normalize dQ/dx 2D YX maps
    for(int y_i = 0; y_i < h2_dqdx_yx[plane_i]->GetNbinsX(); y_i++)
      for(int x_i = 0; x_i < h2_dqdx_yx[plane_i]->GetNbinsY(); x_i++){
        if(yxw[plane_i][y_i][x_i] == 0) continue;
        h2_dqdx_yx[plane_i]->SetBinContent(y_i+1,x_i+1,h2_dqdx_yx[plane_i]->GetBinContent(y_i+1,x_i+1)/yxw[plane_i][y_i][x_i]);
      }
    h2_dqdx_yx[plane_i]->SetStats(0);
    h2_dqdx_yx[plane_i]->SetMinimum(0);
    h2_dqdx_yx[plane_i]->SetMaximum(20);

    // Bin-wise normalize dQ/dx 2D YX maps
    for(int z_i = 0; z_i < h2_dqdx_zx[plane_i]->GetNbinsX(); z_i++)
      for(int x_i = 0; x_i < h2_dqdx_zx[plane_i]->GetNbinsY(); x_i++){
        if(zxw[plane_i][z_i][x_i] == 0) continue;
        h2_dqdx_zx[plane_i]->SetBinContent(z_i+1,x_i+1,h2_dqdx_zx[plane_i]->GetBinContent(z_i+1,x_i+1)/zxw[plane_i][z_i][x_i]);
      }
    h2_dqdx_zx[plane_i]->SetStats(0);
    h2_dqdx_zx[plane_i]->SetMinimum(0);
    h2_dqdx_zx[plane_i]->SetMaximum(20);

    // Bin-wise normalize dQ/dx 2D YX maps
    for(int y_i = 0; y_i < h2_dqdx_yz[plane_i]->GetNbinsX(); y_i++)
      for(int z_i = 0; z_i < h2_dqdx_yz[plane_i]->GetNbinsY(); z_i++){
        if(yzw[plane_i][y_i][z_i] == 0) continue;
        h2_dqdx_yz[plane_i]->SetBinContent(y_i+1,z_i+1,h2_dqdx_yz[plane_i]->GetBinContent(y_i+1,z_i+1)/yzw[plane_i][y_i][z_i]);
      }

    // Set histogram rendering
    h2_dqdx_yz[plane_i]->SetStats(0);
    h2_dqdx_yz[plane_i]->SetMinimum(0);
    h2_dqdx_yz[plane_i]->SetMaximum(20);

    // Write objects to output files
    c_dqdx_z[plane_i]->Write();
    f_dqdx_z[plane_i]->Write();
    g_dqdx_z[plane_i]->Write();
  }

  // End of the game
  std::cout << "Done " << std::endl;
  std::cout << " " << std::endl;
}

float pdvdana::FitdQdx::GetPitch(const recob::Track& track,
                                  const art::Ptr<recob::Hit> hit,
                                  const recob::TrackHitMeta* meta)
{
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  float angleToVert =
    fWireReadoutGeom.WireAngleToVertical(hit->View(), hit->WireID().asPlaneID()) - 0.5 * ::util::pi<>();

  geo::Vector_t dir;

  // "dir" should be the direction that the wires see. If the track already has the field
  // distortion corrections applied, then we need to de-apply them to get the direction as
  // seen by the wire planes
  if (sce->EnableCalSpatialSCE() && bFieldDistortion &&
      bTrackIsFieldDistortionCorrected) {
    geo::Point_t loc = track.LocationAtPoint(meta->Index());

    // compute the dir of the track trajectory
    geo::Vector_t track_dir = track.DirectionAtPoint(meta->Index());
    geo::Point_t loc_mdx = loc - track_dir * (fWireReadoutGeom.Plane({0, 0, hit->View()}).WirePitch() / 2.);
    geo::Point_t loc_pdx = loc + track_dir * (fWireReadoutGeom.Plane({0, 0, hit->View()}).WirePitch() / 2.);

    loc_mdx = TrajectoryToWirePosition(loc_mdx, hit->WireID());
    loc_pdx = TrajectoryToWirePosition(loc_pdx, hit->WireID());

    // Direction at wires
    dir = (loc_pdx - loc_mdx) / (loc_mdx - loc_pdx).r();
  }
  // If there is no space charge or the track is not yet corrected, then the dir
  // is the track is what we want
  else {
    dir = track.DirectionAtPoint(meta->Index());
  }

  float cosgamma = std::abs(std::sin(angleToVert) * dir.Y() + std::cos(angleToVert) * dir.Z());
  float pitch;
  if (cosgamma) { pitch = fWireReadoutGeom.Plane({0, 0, hit->View()}).WirePitch() / cosgamma; }
  else {
    pitch = 0.;
  }

  // now take the pitch computed on the wires and correct it back to the particle trajectory
  geo::Point_t loc_w = GetLocationAtWires(track, hit, meta);

  geo::Point_t locw_pdx_traj = WireToTrajectoryPosition(loc_w + pitch * dir, hit->WireID());
  geo::Point_t loc = WireToTrajectoryPosition(loc_w, hit->WireID());

  pitch = (locw_pdx_traj - loc).R();

  return pitch;
}

void pdvdana::FitdQdx::GetAngles(const recob::Track& track, bool isRevert, float &theta, float &phi){
  // Theta is 90(180) deg. for horizontal(vertical down-going) tracks
  // Phi is 0(90) deg. for tracks along the beam(collection strips) direction

  float x_len = track.End().X()-track.Start().X();
  float y_len = track.End().Y()-track.Start().Y();
  float z_len = track.End().Z()-track.Start().Z();

  if(isRevert){x_len *= -1.; z_len *= -1.;}

  phi = 0;
  if(z_len != 0) phi = abs(atan(y_len/z_len));

  if(     y_len >= 0 && z_len > 0) phi =               phi  / degTOrad;
  else if(y_len >= 0 && z_len < 0) phi = ( TMath::Pi()-phi) / degTOrad;
  else if(y_len < 0  && z_len < 0) phi = (-TMath::Pi()+phi) / degTOrad;
  else if(y_len < 0  && z_len > 0) phi =              -phi  / degTOrad;


  theta = (TMath::Pi()/2. - atan((x_len)/sqrt(pow(y_len,2)+pow(z_len,2))))/ degTOrad;
}

geo::Point_t pdvdana::FitdQdx::TrajectoryToWirePosition(const geo::Point_t& loc,
                                                        const geo::TPCID& tpc)
{
  // See Gnocchi calorimetry module for reference
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  geo::Point_t ret = loc;

  if (sce->EnableCalSpatialSCE() && bFieldDistortion) {
    // Returned X is the drift -- multiply by the drift direction to undo this
    int corr = fGeom->TPC(tpc).DriftDir().X();

    geo::Vector_t offset = sce->GetPosOffsets(ret);

    ret.SetX(ret.X() + corr * bFieldDistortionCorrectionXSign * offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;
}

geo::Point_t pdvdana::FitdQdx::GetLocationAtWires(const recob::Track& track,
                                                  const art::Ptr<recob::Hit> hit,
                                                  const recob::TrackHitMeta* meta)
{
  // See Gnocchi calorimetry module for reference
  geo::Point_t loc = track.LocationAtPoint(meta->Index());
  return bTrackIsFieldDistortionCorrected ? TrajectoryToWirePosition(loc, hit->WireID()) :
                                                     loc;
}

geo::Point_t pdvdana::FitdQdx::WireToTrajectoryPosition(const geo::Point_t& loc,
                                                                const geo::TPCID& tpc)
{
  // See Gnocchi calorimetry module for reference
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  geo::Point_t ret = loc;

  if (sce->EnableCalSpatialSCE() && bFieldDistortion) {
    geo::Vector_t offset = sce->GetCalPosOffsets(ret, tpc.TPC);

    ret.SetX(ret.X() + bFieldDistortionCorrectionXSign * offset.X());
    ret.SetY(ret.Y() + offset.Y());
    ret.SetZ(ret.Z() + offset.Z());
  }

  return ret;
}

std::vector<std::vector<unsigned>> pdvdana::FitdQdx::OrganizeHitsSnippets(
  // See Gnocchi calorimetry module for reference
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const std::vector<const recob::TrackHitMeta*>& thms,
  const recob::Track& track,
  unsigned nplanes,
  vector<int> &nhits)
{
  // In this case, we need to only accept one hit in each snippet
  // Snippets are counted by the Start, End, and Wire. If all these are the same for a hit, then they are on the same snippet.
  //
  // If there are multiple valid hits on the same snippet, we need a way to pick the best one.
  // (TODO: find a good way). The current method is to take the one with the highest charge integral.
  struct HitIdentifier {
    int startTick;
    int endTick;
    int wire;
    float integral;

    // construct
    explicit HitIdentifier(const recob::Hit& hit)
      : startTick(hit.StartTick())
      , endTick(hit.EndTick())
      , wire(hit.WireID().Wire)
      , integral(hit.Integral())
    {}

    // Defines whether two hits are on the same snippet
    inline bool operator==(const HitIdentifier& rhs) const
    {
      return startTick == rhs.startTick && endTick == rhs.endTick && wire == rhs.wire;
    }

    // Defines which hit to pick between two both on the same snippet
    inline bool operator>(const HitIdentifier& rhs) const { return integral > rhs.integral; }
  };

  nhits.resize(nplanes,0);

  std::vector<std::vector<unsigned>> ret(nplanes);
  std::vector<std::vector<HitIdentifier>> hit_idents(nplanes);
  for (unsigned i = 0; i < hits.size(); i++) {
    if (HitIsValid(hits[i], thms[i], track)) {
      HitIdentifier this_ident(*hits[i]);

      // check if we have found a hit on this snippet before
      bool found_snippet = false;
      auto const plane = hits[i]->WireID().Plane;
      for (unsigned j = 0; j < ret[plane].size(); j++) {
        if (this_ident == hit_idents[plane][j]) {
          found_snippet = true;
          nhits[plane]++;
          if (this_ident > hit_idents[plane][j]) {
            ret[plane][j] = i;
            hit_idents[plane][j] = this_ident;
          }
          break;
        }
      }
      if (!found_snippet) {
        ret[plane].push_back(i);
        hit_idents[plane].push_back(this_ident);
        nhits[plane]++;
      }
    }
  }
  return ret;
}

bool pdvdana::FitdQdx::HitIsValid(const art::Ptr<recob::Hit> hit,
                                  const recob::TrackHitMeta* thm,
                                  const recob::Track& track)
{
  // See Gnocchi calorimetry module for reference
  if (thm->Index() == int_max_as_unsigned_int) return false;
  if (!track.HasValidPoint(thm->Index())) return false;
  return true;
}





DEFINE_ART_MODULE(pdvdana::FitdQdx)
