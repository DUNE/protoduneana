////////////////////////////////////////////////////////////////////////
// Class:       GetAnodeCathodePulses 
// Plugin Type: analyzer (Unknown Unknown)
// File:        GetAnodeCathodePulses_module.cc
//
// Generated at Mon Jan 30 17:27:10 2023 by Yoann Kermaidic using cetskelgen
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
#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"

// ROOT
#include <TStyle.h>
#include "TTree.h"
#include "TH2.h"


using std::vector;
using std::string;

namespace {
  constexpr unsigned int int_max_as_unsigned_int{std::numeric_limits<int>::max()};
}

namespace pddpana {
  class GetAnodeCathodePulses;
}


class pddpana::GetAnodeCathodePulses : public art::EDAnalyzer {
public:
  explicit GetAnodeCathodePulses(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GetAnodeCathodePulses(GetAnodeCathodePulses const&) = delete;
  GetAnodeCathodePulses(GetAnodeCathodePulses&&) = delete;
  GetAnodeCathodePulses& operator=(GetAnodeCathodePulses const&) = delete;
  GetAnodeCathodePulses& operator=(GetAnodeCathodePulses&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  float    degTOrad = 3.14159/180.; // rad / deg
  float    e_charge = 1.60E-4;      // fC / e-
  float    calarea  = 14.275E-3;    // ADC x tick / e-

private:

  // Declare member data here.
  int      fLogLevel;
  string   fWireModuleLabel;
  string   fTrackModuleLabel;
  float    fFiducialCut;
  float    fTrackCurvature;
  float    fTrackDelta;
  float    fTrackDeltaStop;
  float    fTrackLenMin;
  float    fTrackLenMax;
  float    fTrackThetaMin;
  float    fTrackThetaMax;

  unsigned fTotalWires;
  unsigned fTotalTracks;
  unsigned int fNplanes;

  string fTrackTypes[2] = {"cathode","anode"};

  // summary tree
  TTree   *fTree;

  unsigned fEventNum;
  unsigned fTrackId;
  int      fTrackTyp;
  float    fLength;
  float    fTheta;
  float    fPhi;
  vector<int>   fNwf;
  vector<float> fRMS;
  vector<float> fDqdx;

  // cm - NP02 active anode volume
  float x_min = -300;
  float x_max = 300;
  float y_min = -100;
  float y_max = 0;
  float z_min = 300;
  float z_max = 400;
  float l_min = 0;
  float l_max = sqrt(pow(x_max-x_min,2)+pow(y_max-y_min,2)+pow(z_max-z_min,2));
 
  // cm - Reco uncertainty margin
  float eps   = 5;

  // cm - hit spacing threshold (holes in mis-reconstructed tracks)
  float hs_thrs = 50;

  // Store valid track start and end
  struct Edge{
    float x;
    float y;
  };

  // detector geometry
  geo::WireReadoutGeom const& fWireReadoutGeom = art::ServiceHandle<geo::WireReadout>()->Get();

  // Store good and bad tracks hits into TH2D
  vector<TH1F*> hAvgWire;

  bool HitIsValid(const art::Ptr<recob::Hit> hit,
                  const recob::TrackHitMeta* thm,
                  const recob::Track& track);
  std::vector<std::vector<unsigned>> OrganizeHitsSnippets(
      const std::vector<art::Ptr<recob::Hit>>& hits,
      const std::vector<const recob::TrackHitMeta*>& thms,
      const recob::Track& track,
      unsigned nplanes);
  void  GetAngles(const recob::Track& track, float &theta, float &phi);
  float GetWidth(TH1F* hist, float thrs);
  float GetPitch(const recob::Track& track,
                 const art::Ptr<recob::Hit> hit,
                 const recob::TrackHitMeta* meta);
};


pddpana::GetAnodeCathodePulses::GetAnodeCathodePulses(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
  fLogLevel( p.get< int >("LogLevel") ),
  fWireModuleLabel(  p.get< std::string  >("WireModuleLabel") ),
  fTrackModuleLabel( p.get< std::string  >("TrackModuleLabel") ),
  fFiducialCut(   p.get< float  >("FiducialCut") ),
  fTrackCurvature(p.get< float  >("TrackCurvature") ), // cm - hit to straight track spacing threshold
  fTrackDelta(    p.get< float  >("TrackDelta") ),
  fTrackDeltaStop(p.get< float  >("TrackDeltaStop") ),
  fTrackLenMin(   p.get< float  >("TrackLenMin") ),
  fTrackLenMax(   p.get< float  >("TrackLenMax") ),
  fTrackThetaMin( p.get< float  >("TrackThetaMin") ),
  fTrackThetaMax( p.get< float  >("TrackThetaMax") )
{}

void pddpana::GetAnodeCathodePulses::analyze(art::Event const& e)
{
  if( fLogLevel >= 1 ) std::cout << "Start analysing file ..." << std::endl;

  fNplanes = fWireReadoutGeom.Nplanes();
 
  bool point_set = false;
  Edge bot_edge;
  Edge top_edge;

  int hitout = 0;         // # of hits outside specs
  int isgood = 0;         // Flag good = 0 or bad = 1 tracks

  // get tracks
  auto Wires    = e.getValidHandle<vector<recob::Wire>>(fWireModuleLabel);
  fTotalWires  += Wires->size();  
  auto Tracks   = e.getValidHandle<vector<recob::Track>>(fTrackModuleLabel);
  fTotalTracks += Tracks->size();  
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmHits(Tracks, e, fTrackModuleLabel);
  
  if( fLogLevel >= 1 ){
    std::cout << "The event contains "<< Wires->size()  << " wires\n";
    std::cout << "The event contains "<< Tracks->size() << " tracks\n";
  }
  
  //
  fEventNum = e.id().event();

  // loop over tracks
  for (unsigned trktyp = 0; trktyp <  sizeof(fTrackTypes) / sizeof(string); ++trktyp) {
    int ntrks[3] = {0,0,0};   // Number of tracks taken into account in the analysis

    string fTrackType = fTrackTypes[trktyp];
    fTrackTyp = 2;
    if(fTrackType == "cathode")    fTrackTyp = 0;
    else if(fTrackType == "anode") fTrackTyp = 1;

    for (unsigned trk = 0; trk < Tracks->size(); ++trk) {
    const recob::Track& track = Tracks->at(trk);

    bot_edge.x = 0;
    bot_edge.y = 0;
    top_edge.x = 0;
    top_edge.y = 0;

    isgood = 0;
    hitout = 0;
    point_set = false;

    fLength = track.Length(0);
    fTheta = 0, fPhi = 0;
    GetAngles(track,fTheta,fPhi);

    if( fLogLevel >= 1 ){
      std::cout << "  #" << trk << " length: "<< fLength <<" cm / theta: " << fTheta << " deg \n";
    }

    // Filter out tracks that are shorter/longer than specified threshold
    if(fLength < fTrackLenMin   || fLength > fTrackLenMax)  {ntrks[2]++; continue;}
    // Filter out tracks that are too horizontal or too vertical to speed-up processing 
    if(fTheta  < fTrackThetaMin || fTheta  > fTrackThetaMax){ntrks[2]++; continue;}      
  
    bool is_out = false;

    if(fTrackType == "anode"){
      // Filter out tracks that are starting outside the top anode plane fiducial volume
      for(size_t hit=0; hit<track.NPoints(); ++hit){
        if(track.LocationAtPoint(hit).X() == -999) continue;
        if(track.LocationAtPoint(hit).Y() - y_min < fFiducialCut || track.LocationAtPoint(hit).Y() - y_min > 100 - fFiducialCut) is_out = true;
        if(track.LocationAtPoint(hit).Z() - z_min < fFiducialCut || track.LocationAtPoint(hit).Z() - z_min > 100 - fFiducialCut) is_out = true;
        break;
      }
    }
    else if(fTrackType == "cathode"){
      // Filter out tracks that are ending outside the bottom cathode plane fiducial volume
      for(size_t hit=track.NPoints()-1; hit>0; --hit){
        if(track.LocationAtPoint(hit).X() == -999) continue;
        if(track.LocationAtPoint(hit).Y() - y_min < fFiducialCut || track.LocationAtPoint(hit).Y() - y_min > 100 - fFiducialCut) is_out = true;
        if(track.LocationAtPoint(hit).Z() - z_min < fFiducialCut || track.LocationAtPoint(hit).Z() - z_min > 100 - fFiducialCut) is_out = true;
        break;
      }
    }

    if(is_out){ ntrks[2]++; continue;}

    if(fLength < l_min-eps || fLength  > l_max+eps){
      if(fLogLevel >= 1) std::cout << "  BAD: Length [" << l_min-eps << "," << l_max+eps << "] -> " <<  fLength << std::endl;
      isgood = 1;
    }

    vector<bool> hitskip(track.NPoints(),false);

    // Sort tracks by goodness fTrackType defined by criteria below
    for(size_t hit=0; hit<track.NPoints(); ++hit){
      if(isgood) break;
      if(track.LocationAtPoint(hit).X() == -999) {hitskip[hit] = true; continue;}

      if(!point_set){
        bot_edge.x = track.LocationAtPoint(hit).X() - x_min;
        bot_edge.y = sqrt(pow(track.LocationAtPoint(hit).Y() - y_min,2) + pow(track.LocationAtPoint(hit).Z() - z_min,2));
        point_set = true;
      }

      top_edge.x = track.LocationAtPoint(hit).X() - x_min;
      top_edge.y = sqrt(pow(track.LocationAtPoint(hit).Y() - y_min,2) + pow(track.LocationAtPoint(hit).Z() - z_min,2));

      float hit_spacing = 0;
      if(hit>0) hit_spacing = sqrt(pow(track.LocationAtPoint(hit).X()-track.LocationAtPoint(hit-1).X(),2)
                                 + pow(track.LocationAtPoint(hit).Y()-track.LocationAtPoint(hit-1).Y(),2)
                                 + pow(track.LocationAtPoint(hit).Z()-track.LocationAtPoint(hit-1).Z(),2));

      if(track.LocationAtPoint(hit).X() < x_min ||
         track.LocationAtPoint(hit).X() > x_max ||
         track.LocationAtPoint(hit).Y() < y_min-eps ||
         track.LocationAtPoint(hit).Y() > y_max+eps ||
         track.LocationAtPoint(hit).Z() < z_min-eps ||
         track.LocationAtPoint(hit).Z() > z_max+eps) {hitskip[hit] = true; hitout++;}

      if(hitout > 5 || hit_spacing > hs_thrs) {
        std::cout << "  BAD: " << hit << " " << hitout << " "
                          << track.LocationAtPoint(hit).X() << " "
                          << track.LocationAtPoint(hit).Y() << " "
                          << track.LocationAtPoint(hit).Z() << std::endl;
        isgood = 1;
        break;
      }
    }

    // Compute average track slope from the two track extremum points
    // later used to compute the minimal distance between hit point and line from extremums 

    if(top_edge.x-bot_edge.x == 0) isgood = 1;
    else if(!isgood){
      float a = 666.; a = (top_edge.y-bot_edge.y)/(top_edge.x-bot_edge.x);
      float b = 666.; b = bot_edge.y - a*bot_edge.x;

      for(size_t hit=0; hit<track.NPoints(); ++hit){
        if(hitskip[hit]) continue;

        float xpt = track.LocationAtPoint(hit).X() - x_min;
        float yzpt= sqrt(pow(track.LocationAtPoint(hit).Y() - y_min,2) + pow(track.LocationAtPoint(hit).Z() - z_min,2));

        float xx  = 1./(a*a+1)*(xpt+a*yzpt-a*b);
        float yy  = 1./a*(xpt-xx) + yzpt;

        if(sqrt(pow(xx-xpt,2) + pow(yy-yzpt,2)) > fTrackCurvature) {
          std::cout << "  BAD: " << hit << " distance: " << sqrt(pow(xx-xpt,2) + pow(yy-yzpt,2)) << " " << xpt << " " << yzpt << std::endl;
          isgood = 1;
          break;
        }
      }
    }

      if(fLogLevel >= 1)
        std::cout << "Event: #" << e.id().event() << ": track length: " << fLength << " cm / theta: " << fTheta/degTOrad << " deg found -> " << isgood << " status" << std::endl;

      ntrks[isgood]++;

      if(isgood) continue;

      // if the track passed selection criteria, loop over hits to retrieve pulse width
      const std::vector<art::Ptr<recob::Hit>>& Hits       = fmHits.at(trk);
      const std::vector<const recob::TrackHitMeta*>& thms = fmHits.data(trk);
 
      // follow GnocchiCalorimatery module to sort hits by plane and snippets
      std::vector<std::vector<unsigned>> hit_indices = OrganizeHitsSnippets(Hits, thms, track, fNplanes);

      vector<TH1F*> hWire(fNplanes); 

      // Loop over planes
      for (unsigned plane_i = 0; plane_i < fNplanes; plane_i++){
        // Store track average pulse to retrieve its RMS + number of wf used
        fNwf[plane_i]  = 0;
        fRMS[plane_i]  = -1.;
        fDqdx[plane_i] = -1.;
        hWire[plane_i] = new TH1F(Form("hWire_%d",plane_i),"",100,-50,50);

        // Loop over hits
        for (unsigned hit_i = 0; hit_i < hit_indices[plane_i].size(); hit_i++) {
          unsigned hit_index = hit_indices[plane_i][hit_i];

          // Find hits falling in the proper channel with multiplicity 1
          for (unsigned chan_i = 0; chan_i < Wires->size(); chan_i++){
            if(Hits[hit_index]->Channel() != Wires->at(chan_i).Channel()) continue;
            if(Hits[hit_index]->Multiplicity() != 1)                      continue;
            if(Hits[hit_index]->GoodnessOfFit() > 1)                      continue;

            float pitch  = GetPitch(track, Hits[hit_index], thms[hit_index]);
            float charge = Hits[hit_index]->ROISummedADC() * e_charge/calarea;
            fDqdx[plane_i] += charge/pitch;

            vector<float> wire_signal = Wires->at(chan_i).Signal();
            for (int iadc=Hits[hit_index]->StartTick();iadc<Hits[hit_index]->EndTick();++iadc){
              hWire[plane_i]->Fill(iadc-Hits[hit_index]->PeakTime(),wire_signal[iadc]);
              hAvgWire[plane_i]->Fill(iadc-Hits[hit_index]->PeakTime(),wire_signal[iadc]);
            }
            fNwf[plane_i]++;
            break;
          }
        } // end of hit loop

        // Compute FWHM and translate it to sigma (avoid fitting a gaussian)
        if(fNwf[plane_i] > 0){
          fRMS[plane_i] = GetWidth(hWire[plane_i],0.5)/(2.*sqrt(2.*log(2.)));
          fDqdx[plane_i] /= fNwf[plane_i]++;
        }

      } // end of plane loop
      fTree->Fill();
    } // end of tracks loop

    if(fLogLevel >= 1){
      std::cout << "Total # of " << fTrackType << " tracks: " << ntrks[0]+ntrks[1] << std::endl;
      std::cout << "  - good: " << ntrks[0] << std::endl;
      std::cout << "  - bad:  " << ntrks[1] << std::endl;
      std::cout << "  - out:  " << ntrks[2] << std::endl;
      std::cout << " " << std::endl;
    }
  } // end of track types loop
}

void pddpana::GetAnodeCathodePulses::beginJob()
{
  if(fLogLevel >= 1){
    std::cout << " #planes:       " << fWireReadoutGeom.Nplanes() << std::endl;
    std::cout << " TrackMinLen:   " << fTrackLenMin << std::endl;
    std::cout << " TrackMaxLen:   " << fTrackLenMax << std::endl;
    std::cout << " TrackThetaMin: " << fTrackThetaMin << std::endl;
    std::cout << " TrackThetaMax: " << fTrackThetaMax << std::endl;
  }

  hAvgWire.resize(fWireReadoutGeom.Nplanes());  
  fNwf.resize(fWireReadoutGeom.Nplanes());
  fRMS.resize(fWireReadoutGeom.Nplanes());

  art::ServiceHandle<art::TFileService> tfs;
  for (unsigned plane_i = 0; plane_i < fWireReadoutGeom.Nplanes(); plane_i++)
    hAvgWire[plane_i] = tfs->make<TH1F>(Form("hAvgWire_%d",plane_i),Form("Plane %d pulses",plane_i),100,-50,50);

  fTree = tfs->make<TTree>("DiffusionTree","Store anode and cathode pulses");
  fTree->Branch("EventNum",  &fEventNum,"EventNum/i"  );
  fTree->Branch("TrackId",   &fTrackId, "TrackId/i"   );
  fTree->Branch("TrackTyp",  &fTrackTyp,"TrackTyp/i"  );
  fTree->Branch("TrackLen",  &fLength,  "TrackLen/F"  );
  fTree->Branch("TrackTheta",&fTheta,   "TrackTheta/F");
  fTree->Branch("TrackPhi",  &fPhi,     "TrackPhi/F"  );
  fTree->Branch("TrackNwf",  &fNwf);
  fTree->Branch("TrackRMS",  &fRMS);
  fTree->Branch("TrackDqdx", &fDqdx);
}

void pddpana::GetAnodeCathodePulses::endJob()
{
  std::cout << "Done " << std::endl;
  std::cout << " " << std::endl;
}

void pddpana::GetAnodeCathodePulses::GetAngles(const recob::Track& track, float &theta, float &phi){

  float x_fst=-1, y_fst=-1, z_fst=-1;
  float x_lst=-1, y_lst=-1, z_lst=-1;

  for(size_t hit=0; hit<track.NPoints(); ++hit){
    if(track.LocationAtPoint(hit).X() == -999) continue;
    if(x_fst == -1){
      x_fst = track.LocationAtPoint(hit).X();
      y_fst = track.LocationAtPoint(hit).Y();
      z_fst = track.LocationAtPoint(hit).Z();
    }
    x_lst = track.LocationAtPoint(hit).X();
    y_lst = track.LocationAtPoint(hit).Y();
    z_lst = track.LocationAtPoint(hit).Z();
  }

  float x_len = x_lst-x_fst;
  float y_len = y_lst-y_fst;
  float z_len = z_lst-z_fst;

  phi = 0;
  if(     y_len >= 0 && z_len > 0) phi =                atan(y_len/z_len) / degTOrad;
  else if(y_len >= 0 && z_len < 0) phi = TMath::Pi()   -atan(y_len/z_len) / degTOrad;
  else if(y_len < 0  && z_len < 0) phi = TMath::Pi()   +atan(y_len/z_len) / degTOrad;
  else if(y_len < 0  && z_len > 0) phi = TMath::TwoPi()-atan(y_len/z_len) / degTOrad;

  theta = atan((x_len)/sqrt(pow(y_len,2)+pow(z_len,2))) / degTOrad;
}

bool pddpana::GetAnodeCathodePulses::HitIsValid(const art::Ptr<recob::Hit> hit,
                const recob::TrackHitMeta* thm,
                const recob::Track& track)
{
  if(!thm)                                     return false;
  if (thm->Index() == int_max_as_unsigned_int) return false;
  if (!track.HasValidPoint(thm->Index()))      return false;
  return true;
}


std::vector<std::vector<unsigned>> pddpana::GetAnodeCathodePulses::OrganizeHitsSnippets(
  const std::vector<art::Ptr<recob::Hit>>& hits,
  const std::vector<const recob::TrackHitMeta*>& thms,
  const recob::Track& track,
  unsigned nplanes)
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

  std::vector<std::vector<unsigned>> ret(nplanes);
  std::vector<std::vector<HitIdentifier>> hit_idents(nplanes);
  for (unsigned i = 0; i < hits.size(); i++) {
    if (HitIsValid(hits[i], thms[i], track)) {
      HitIdentifier this_ident(*hits[i]);

      // check if we have found a hit on this snippet before
      bool found_snippet = false;
      for (unsigned j = 0; j < ret[hits[i]->WireID().Plane].size(); j++) {
        if (this_ident == hit_idents[hits[i]->WireID().Plane][j]) {
          found_snippet = true;
          if (this_ident > hit_idents[hits[i]->WireID().Plane][j]) {
            ret[hits[i]->WireID().Plane][j] = i;
            hit_idents[hits[i]->WireID().Plane][j] = this_ident;
          }
          break;
        }
      }
      if (!found_snippet) {
        ret[hits[i]->WireID().Plane].push_back(i);
        hit_idents[hits[i]->WireID().Plane].push_back(this_ident);
      }
    }
  }
  return ret;
}

float pddpana::GetAnodeCathodePulses::GetWidth(TH1F* hist, float thrs = 0.5){
  int binmax = hist->GetMaximumBin();
  float max = hist->GetBinContent(binmax);

//  float lowamp1 = 1, highamp1 = 1;
//  float lowamp2 = 1, highamp2 = 1;
  int   lowbin1 = 0, highbin1 = 0;
  int   lowbin2 = 0, highbin2 = 0;

  for(int bin = binmax; bin > 0; bin--){
    if(hist->GetBinContent(bin) / max >= thrs) lowbin2 = bin;
    if(hist->GetBinContent(bin) / max < thrs) {lowbin1 = bin; break;}
  }
  for(int bin = binmax; bin < hist->GetNbinsX(); bin++){
    if(hist->GetBinContent(bin) / max >= thrs)  highbin1 = bin-1;
    if(hist->GetBinContent(bin) / max <  thrs) {highbin2 = bin-1; break;}
  }

  float lowedge = hist->GetBinCenter(lowbin1) + (hist->GetBinCenter(lowbin2)-hist->GetBinCenter(lowbin1))/(hist->GetBinContent(lowbin2)-hist->GetBinContent(lowbin1))*(thrs*max-hist->GetBinContent(lowbin1));

  float highedge = hist->GetBinCenter(highbin1) + (hist->GetBinCenter(highbin2)-hist->GetBinCenter(highbin1))/(hist->GetBinContent(highbin2)-hist->GetBinContent(highbin1))*(thrs*max-hist->GetBinContent(highbin1));

  return highedge-lowedge;
}

float pddpana::GetAnodeCathodePulses::GetPitch(const recob::Track& track,
                const art::Ptr<recob::Hit> hit,
                const recob::TrackHitMeta* meta)
{ 
  float angleToVert;
  if(hit->View() == 3) angleToVert = -0.5 * 3.14159;
  else                 angleToVert = 0;
  
  geo::Vector_t dir = track.DirectionAtPoint(meta->Index());

  
  float cosgamma = std::abs(std::sin(angleToVert) * dir.Y() + std::cos(angleToVert) * dir.Z());
  float pitch; 
  if (cosgamma) { pitch = 0.3125 / cosgamma; }
  else  pitch = 0.;

//  cout << "  " << hit->View() << " " << cosgamma << " " << pitch << " " << dir.X() << " " << dir.Y() << " " << dir.Z() << endl;
  
  geo::Point_t loc = track.LocationAtPoint(meta->Index());
  geo::Point_t locw_pdx_traj = loc + pitch * dir;
  
  return (locw_pdx_traj - loc).R();
}


DEFINE_ART_MODULE(pddpana::GetAnodeCathodePulses)
