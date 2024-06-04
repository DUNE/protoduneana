////////////////////////////////////////////////////////////////////////
// Class:       RecoEff
// Module Type: analyzer
// File:        RecoEff_module.cc
// Author:      D.Stefan and R.Sulej 
// Modified by: Jake Calcutt 
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"
#include "TTree.h"
#include "TVector3.h"


/*#ifndef PARTTREE
  #define PARTTREE
#endif

#ifndef IDETREE
  #define IDETREE
#endif
*/

namespace pdune
{
	class RecoEff;
}

class pdune::RecoEff : public art::EDAnalyzer {
public:

  enum EFilterMode { kCosmic, kBeam, kPrimary, kSecondary };

  struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> SimulationLabel { Name("SimulationLabel"), Comment("tag of simulation producer") };

      fhicl::Atom<art::InputTag> HitModuleLabel { Name("HitModuleLabel"), Comment("tag of hits producer") };

      fhicl::Atom<art::InputTag> TrackModuleLabel { Name("TrackModuleLabel"), Comment("tag of track producer") };
    
      fhicl::Atom<size_t> MinHitsPerPlane { Name("MinHitsPerPlane"), Comment("min hits per plane for a reconstructable MC particle") };

      fhicl::Sequence<std::string> Filters { Name("Filters"), Comment("on which particles the efficiency is calculated") };

      fhicl::Sequence<int> Pdg { Name("Pdg"), Comment("which PDG code, or 0 for all particles") };
      
      fhicl::Atom<size_t> EffHitMax { Name("EffHitMax"), Comment("max hits per MC particle in the track efficiency") };
      fhicl::Atom<size_t> EffHitBins { Name("EffHitBins"), Comment("number of bins in the track efficiency") };

      fhicl::Atom<size_t> EffEMax { Name("EffEMax"), Comment("max hits per MC particle in the track efficiency") };
      fhicl::Atom<size_t> EffEBins { Name("EffEBins"), Comment("number of bins in the track efficiency") };

      fhicl::Atom<size_t> EffDepMax { Name("EffDepMax"), Comment("max hits per MC particle in the track efficiency") };
      fhicl::Atom<size_t> EffDepBins { Name("EffDepBins"), Comment("number of bins in the track efficiency") };

      fhicl::Atom<size_t> EffLengthMax { Name("EffLengthMax"), Comment("max hits per MC particle in the track efficiency") };
      fhicl::Atom<size_t> EffLengthBins { Name("EffLengthBins"), Comment("number of bins in the track efficiency") };

      fhicl::Atom<size_t> EnableParticleTree { Name("EnableParticleTree"), Comment("Turn on saving Particle Tree") };
      fhicl::Atom<size_t> EnableIDETree { Name("EnableIDETree"), Comment("Turn on saving IDE Tree") };

  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit RecoEff(Parameters const& config);

  RecoEff(RecoEff const &) = delete;
  RecoEff(RecoEff &&) = delete;
  RecoEff & operator = (RecoEff const &) = delete;
  RecoEff & operator = (RecoEff &&) = delete;

  void analyze(art::Event const & evt) override;
  void beginJob() override;
  void endJob() override;

private:
  void ResetVars();
  double GetLengthInTPC(detinfo::DetectorClocksData const& clockData,
                        simb::MCParticle part);
  TVector3 GetPositionInTPC(detinfo::DetectorClocksData const& clockData,
                            simb::MCParticle part, int MCHit);
  double GetdEdX(simb::MCTrajectory traj, size_t this_step);
  
  TH1D* fDenominatorHist;
  TH1D* fNominatorHist;

  TH1D* fEnergyDenominatorHist;
  TH1D* fEnergyNominatorHist;

  TH1D* fLengthDenominatorHist;
  TH1D* fLengthNominatorHist;

  TH1D* fDepDenominatorHist;
  TH1D* fDepNominatorHist;

  TTree *fEventTree;
  TTree *fTrkTree;
  //TTree *fHitTree; // unused
  
  TTree *fTrkIDETree;
  TTree *fParticleTree;

  int fRun, fEvent;
  short fNRecoTracks;
  short fReconstructable;
  short fMatched;

  float fEnGen, fEkGen;
  float fT0;

  float fTrkPurityPerPlane[3], fTrkCompletnessPerPlane[3]; // fixed size, don't expect more than 3 planes, less is OK
  float fTrkPurity, fTrkCompletness; // all planes together
  int fTrkPid;
  short fTrkSize;
  short fTrkMatched;
  double fTrkLength;
  int fTrkID;
  double fTrkEnergy;
  double fTrkX[10000];
  double fTrkY[10000]; 
  double fTrkZ[10000];
  double fT0Corr[10000];
  int nT0s;

  std::vector<int> fTrkIDEs;
  int fTrkNIDEs;
  std::vector<double> fTrkContribs;

  double fTrueTrkE;
  double fTrueTrkLength;
  double fTrueTrkDep;
  int fTrueTrkPID;
  int fTrueTrkID;
  int fTrueTrkSize;
  int track_step;
  double fTrueTrkdEdX[10000];
  double fTrueTrkStepE[10000];
  double fTrueTrkX[10000];
  double fTrueTrkY[10000]; 
  double fTrueTrkZ[10000];
  int fTrueTrkOrigin;
  bool fTrueTrkCrossCPA;


  int fTrueTrkRecoGuess;
  double fTrueTrkEnergyGuess;
  std::vector<int> fTrueTrkRecos;
  int fTrueTrkNRecos;
  std::vector<double> fTrueTrkContribs;

  double MCTruthT0, TickT0;
  TH1D* fHitDx;
  TH1D* fHitDy;
  TH1D* fHitDz;
  TH1D* fHitDist3D;

  art::InputTag fSimulationLabel;
  art::InputTag fHitModuleLabel;
  art::InputTag fTrackModuleLabel;
  art::InputTag fTempTrackModuleLabel;

  size_t fMinHitsPerPlane;

  std::vector< EFilterMode > fFilters;
  std::vector< int > fPdg;

  size_t fEffHitMax, fEffHitBins;
  size_t fEffEMax, fEffEBins;
  size_t fEffDepMax, fEffDepBins;
  size_t fEffLengthMax, fEffLengthBins;
  bool  fEnableParticleTree, fEnableIDETree;

  
  int fPartStatus;
  int fPartPID;
  int fPartID;
  int fOrigin;
  double fPartLength;
  double fTotLength;
  std::string fProcess;
  std::string fEndProcess;

  art::ServiceHandle<geo::Geometry> geom;
  geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  double XDriftVelocity{};
  double WindowSize{};
};

pdune::RecoEff::RecoEff(Parameters const& config) : EDAnalyzer(config),
    fSimulationLabel(config().SimulationLabel()),
    fHitModuleLabel(config().HitModuleLabel()),
    fTrackModuleLabel(config().TrackModuleLabel()),


    fMinHitsPerPlane(config().MinHitsPerPlane()),
    fPdg(config().Pdg()),
    fEffHitMax(config().EffHitMax()),
    fEffHitBins(config().EffHitBins()),

    fEffEMax(config().EffEMax()),
    fEffEBins(config().EffEBins()),

    fEffDepMax(config().EffDepMax()),
    fEffDepBins(config().EffDepBins()),

    fEffLengthMax(config().EffLengthMax()),
    fEffLengthBins(config().EffLengthBins()),

    fEnableParticleTree(config().EnableParticleTree()),
    fEnableIDETree(config().EnableIDETree())
//    fEnableParticleTree(config().EnableParticleTree),
//    fEnableIDETree(config().EnableIDETree)

{
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);
  XDriftVelocity = detProp.DriftVelocity()*1e-3; //cm/ns
  WindowSize = detProp.NumberTimeSamples() * clockData.TPCClock().TickPeriod() * 1e3;

    auto flt = config().Filters();
    for (const auto & s : flt)
    {
        if (s == "cosmic")       { fFilters.push_back(pdune::RecoEff::kCosmic);  }
        else if (s == "beam")    { fFilters.push_back(pdune::RecoEff::kBeam);    }
        else if (s == "primary") { fFilters.push_back(pdune::RecoEff::kPrimary); }
        else if (s == "second")  { fFilters.push_back(pdune::RecoEff::kSecondary); }
        else { std::cout << "unsupported filter name: " << s << std::endl; }
    }
}

void pdune::RecoEff::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	
	fDenominatorHist = tfs->make<TH1D>("Denominator", "all reconstructable particles", fEffHitBins, 0., fEffHitMax);
	fNominatorHist = tfs->make<TH1D>("Nominator", "reconstructed and matched tracks", fEffHitBins, 0., fEffHitMax);

	fEnergyDenominatorHist = tfs->make<TH1D>("EnergyDenominator", "all reconstructable particles", fEffEBins, 0., fEffEMax);
	fEnergyNominatorHist = tfs->make<TH1D>("EnergyNominator", "reconstructed and matched tracks",  fEffEBins, 0., fEffEMax);

	fLengthDenominatorHist = tfs->make<TH1D>("LengthDenominator", "all reconstructable particles", fEffLengthBins, 0., fEffLengthMax);
	fLengthNominatorHist = tfs->make<TH1D>("LengthNominator", "reconstructed and matched tracks",  fEffLengthBins, 0., fEffLengthMax);

	fDepDenominatorHist = tfs->make<TH1D>("DepDenominator", "all reconstructable particles", fEffDepBins, 0., fEffDepMax);
	fDepNominatorHist = tfs->make<TH1D>("DepNominator", "reconstructed and matched tracks",  fEffDepBins, 0., fEffDepMax);

	fEventTree = tfs->make<TTree>("events", "summary tree");
	fEventTree->Branch("fRun", &fRun, "fRun/I");
	fEventTree->Branch("fEvent", &fEvent, "fEvent/I");
	fEventTree->Branch("fEnGen", &fEnGen, "fEnGen/F");
	fEventTree->Branch("fEkGen", &fEkGen, "fEkGen/F");
	fEventTree->Branch("fT0", &fT0, "fT0/F");
	fEventTree->Branch("fNRecoTracks", &fNRecoTracks, "fNRecoTracks/S");
	fEventTree->Branch("fReconstructable", &fReconstructable, "fReconstructable/S");
	fEventTree->Branch("fMatched", &fMatched, "fMatched/S");

	fTrkTree = tfs->make<TTree>("tracks", "track metrics");
	fTrkTree->Branch("fTrkPurity", &fTrkPurity, "fTrkPurity/F");
	fTrkTree->Branch("fTrkCompletness", &fTrkCompletness, "fTrkCompletness/F");
	fTrkTree->Branch("fTrkPurityPerPlane", &fTrkPurityPerPlane, "fTrkPurityPerPlane[3]/F");
	fTrkTree->Branch("fTrkCompletnessPerPlane", &fTrkCompletnessPerPlane, "fTrkCompletnessPerPlane[3]/F");
	fTrkTree->Branch("fTrkPid", &fTrkPid, "fTrkPid/I");
	fTrkTree->Branch("fTrkID", &fTrkID, "fTrkID/I");
	fTrkTree->Branch("fTrkSize", &fTrkSize, "fTrkSize/S");
	fTrkTree->Branch("fTrkMatched", &fTrkMatched, "fTrkMatched/S");
        fTrkTree->Branch("fTrkLength", &fTrkLength, "fTrkLength/D");
        fTrkTree->Branch("fTrkIDEs", &fTrkIDEs);
        fTrkTree->Branch("fTrkContribs", &fTrkContribs);
        fTrkTree->Branch("fTrkNIDEs", &fTrkNIDEs);
        fTrkTree->Branch("fEvent", &fEvent);
        fTrkTree->Branch("fTrkEnergy", &fTrkEnergy);
        fTrkTree->Branch("fTrkX", &fTrkX, "fTrkX[10000]/D");
        fTrkTree->Branch("fTrkY", &fTrkY, "fTrkY[10000]/D");
        fTrkTree->Branch("fTrkZ", &fTrkZ, "fTrkZ[10000]/D");
        fTrkTree->Branch("fT0Corr", &fT0Corr, "fT0Corr[10000]/D");
        fTrkTree->Branch("nT0s", &nT0s);

        if(fEnableParticleTree){
          fParticleTree = tfs->make<TTree>("particles","Particles");

          fParticleTree->Branch("fPID",&fPartPID);
          fParticleTree->Branch("fStatus",&fPartStatus);
          fParticleTree->Branch("fID",&fPartID);
          fParticleTree->Branch("fLength",&fPartLength);
          fParticleTree->Branch("fTotLength",&fTotLength);
          fParticleTree->Branch("fEvent",&fEvent);
          fParticleTree->Branch("fOrigin",&fOrigin);
          fParticleTree->Branch("fProcess",&fProcess);
          fParticleTree->Branch("fEndProcess",&fEndProcess);
        }
 
        if(fEnableIDETree){
          fTrkIDETree = tfs->make<TTree>("trackIDEs","trackIDE metrics");

          fTrkIDETree->Branch("fTrkE", &fTrueTrkE, "fTrkE/D");
          fTrkIDETree->Branch("fTrkLength", &fTrueTrkLength, "fTrkLength/D");
          fTrkIDETree->Branch("fTrkDep", &fTrueTrkDep, "fTrkDep/D");
          fTrkIDETree->Branch("fTrkPID", &fTrueTrkPID, "fTrkPID/I");
          fTrkIDETree->Branch("fTrkID", &fTrueTrkID, "fTrkID/I");
          fTrkIDETree->Branch("fTrkSize", &fTrueTrkSize, "fTrkSize/I");
          fTrkIDETree->Branch("fTrkX", &fTrueTrkX, "fTrkX[10000]/D");
          fTrkIDETree->Branch("fTrkY", &fTrueTrkY, "fTrkY[10000]/D");
          fTrkIDETree->Branch("fTrkZ", &fTrueTrkZ, "fTrkZ[10000]/D");
          fTrkIDETree->Branch("fT0Corr", &fT0Corr, "fT0Corr[10000]/D");
          fTrkIDETree->Branch("fTrkdEdX", &fTrueTrkdEdX, "fTrkdEdX[10000]/D");
          fTrkIDETree->Branch("fTrkStepE", &fTrueTrkStepE, "fTrkStepE[10000]/D");
          fTrkIDETree->Branch("fTrkSteps", &track_step, "fTrkSteps/D");
	  fTrkIDETree->Branch("fEvent", &fEvent, "fEvent/I");
          fTrkIDETree->Branch("fTrkOrigin",&fTrueTrkOrigin);
          fTrkIDETree->Branch("fTrkRecoGuess", &fTrueTrkRecoGuess);
          fTrkIDETree->Branch("fTrkEnergyGuess",&fTrueTrkEnergyGuess);
          fTrkIDETree->Branch("fTrkRecos", &fTrueTrkRecos);
          fTrkIDETree->Branch("fTrkNRecos", &fTrueTrkNRecos);
          fTrkIDETree->Branch("fTrkContribs", &fTrueTrkContribs);
          fTrkIDETree->Branch("fTrkCrossCPA", &fTrueTrkCrossCPA);
        }

    fHitDist3D = tfs->make<TH1D>("HitD3D", "MC-reco 3D distance", 400, 0., 10.0);
    fHitDx = tfs->make<TH1D>("HitDx", "MC-reco X distance", 400, 0., 10.0);
    fHitDy = tfs->make<TH1D>("HitDy", "MC-reco Y distance", 400, 0., 10.0);
    fHitDz = tfs->make<TH1D>("HitDz", "MC-reco Z distance", 400, 0., 10.0);

    //Need to correct for name in space point branch
    if(fTrackModuleLabel == "pandoraTrack"){
      fTempTrackModuleLabel = "pandora";
    }
    else{
       fTempTrackModuleLabel = fTrackModuleLabel; 
    }
}

void pdune::RecoEff::endJob()
{
    for (int i = 0; i < fDenominatorHist->GetNbinsX(); ++i)
    {
        double nom = fNominatorHist->GetBinContent(i);
        double den = fDenominatorHist->GetBinContent(i);
        std::cout << "nhits: " << i << " nom:" << nom << " den:" << den << std::endl;
        if (nom > den) { fNominatorHist->SetBinContent(i, den); }
    }
    art::ServiceHandle<art::TFileService> tfs;
    TEfficiency* fEfficiency = tfs->makeAndRegister<TEfficiency>("Efficiency", "tracking efficiency", *fNominatorHist, *fDenominatorHist);
    TEfficiency* fEnergyEfficiency = tfs->makeAndRegister<TEfficiency>("EnergyEfficiency", "tracking efficiency", *fEnergyNominatorHist, *fEnergyDenominatorHist);
    std::cout << "Efficiency created: " << fEfficiency->GetTitle() << " " << fEnergyEfficiency->GetTitle() << std::endl;
    TEfficiency* fLengthEfficiency = tfs->makeAndRegister<TEfficiency>("LengthEfficiency", "tracking efficiency", *fLengthNominatorHist, *fLengthDenominatorHist);
    std::cout << "Efficiency created: " << fEfficiency->GetTitle() << " " << fLengthEfficiency->GetTitle() << std::endl;
    TEfficiency* fDepEfficiency = tfs->makeAndRegister<TEfficiency>("DepEfficiency", "tracking efficiency", *fDepNominatorHist, *fDepDenominatorHist);
    std::cout << "Efficiency created: " << fEfficiency->GetTitle() << " " << fDepEfficiency->GetTitle() << std::endl;

}

void pdune::RecoEff::analyze(art::Event const & evt)
{
  ResetVars();
  
  fRun = evt.run();
  fEvent = evt.id().event();

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

  if(fEnableParticleTree){
    //Particle List diagnostics
    const sim::ParticleList& plist = pi_serv->ParticleList();  
    for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){ 

        simb::MCParticle * part = ipar->second;

        fPartPID = part->PdgCode();
        if(fPartPID > 2212) {continue;}

        fPartStatus = part->StatusCode();
        fPartID     = part->TrackId();
        fPartLength = GetLengthInTPC(clockData, *part);
        fTotLength  = part->Trajectory().TotalLength();
        fProcess    = part->Process();
        fEndProcess = part->EndProcess();
        fOrigin     = pi_serv->ParticleToMCTruth_P(part)->Origin();

        fParticleTree->Fill();
    }
  }

  // we are going to look only for these MC truth particles, which contributed to hits
  // and normalize efficiency to things which indeed generated activity and hits in TPC's:
  // - need a vector of hits associated to these MC particles,
  // - need the energy deposit corresponding to the particle part seen i  these hits.
  std::unordered_map<int, std::vector<recob::Hit> > mapTrackIDtoHits;
  std::unordered_map<int, double> mapTrackIDtoHitsEnergyPerPlane[3];
  std::unordered_map<int, double> mapTrackIDtoHitsEnergy;

  std::unordered_map<int, double> mapTrackIDtoLength;
  const auto hitListHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
  for (auto const & h : *hitListHandle)
  {
    std::unordered_map<int, double> particleID_E;

    for (auto const & id : bt_serv->HitToTrackIDEs(clockData,h)) // loop over std::vector< sim::TrackIDE > contributing to hit h
	{

        // select only hadronic and muon track, skip EM activity (electron's pdg, negative track id)
        if ((id.trackID > 0) && (abs((pi_serv->TrackIdToParticle_P(id.trackID))->PdgCode()) != 11))
        {
            particleID_E[id.trackID] += id.energy;
            const simb::MCParticle * part = pi_serv->TrackIdToParticle_P(id.trackID);
            mapTrackIDtoLength[id.trackID] = GetLengthInTPC(clockData, *part);
        }        
    }

    int best_id = 0;
    double max_e = 0;
    for (auto const & contrib : particleID_E)
    {
        if (contrib.second > max_e) // find particle TrackID corresponding to max energy contribution
        {
            max_e = contrib.second;
            best_id = contrib.first;
        }
    }

    if (max_e > 0)
    {
        mapTrackIDtoHits[best_id].push_back(h);
        mapTrackIDtoHitsEnergyPerPlane[h.WireID().Plane][best_id] += max_e;
        mapTrackIDtoHitsEnergy[best_id] += max_e;
    }
  }

  // ---------------------------------------------------------------------------------------


  // now lets map particles to their associated hits, but with some filtering of things
  // which have a minimal chance to be reconstructed: so at least some hits in two views
  std::unordered_map<int, std::vector< recob::Hit >> mapTrackIDtoHits_filtered;
  for (auto const & p : mapTrackIDtoHits)
  {
    bool skip = false;
    auto origin = pi_serv->TrackIdToMCTruth_P(p.first)->Origin();
    for (auto f : fFilters)
    {
        switch (f)
        {
            case pdune::RecoEff::kCosmic:  if (origin != simb::kCosmicRay) { skip = true; } break;
            case pdune::RecoEff::kBeam:    if (origin != simb::kSingleParticle) { skip = true; } break;

            case pdune::RecoEff::kPrimary:
                if (pi_serv->TrackIdToParticle_P(p.first)->Process() != "primary") { skip = true; }
                break;

            case pdune::RecoEff::kSecondary:
                {
                    int mId = pi_serv->TrackIdToParticle_P(p.first)->Mother();
                    const simb::MCParticle * mother = pi_serv->TrackIdToParticle_P(mId);
                    if ((mother == 0) || (mother->Process() != "primary")) { skip = true; }
                    break;
                }

            default: break;
        }
    }
    if (skip) { continue; } // skip if any condition failed

    if (!fPdg.empty())
    {
        skip = true;
        for (int code : fPdg) { if (pi_serv->TrackIdToParticle_P(p.first)->PdgCode() == code) { skip = false; break; } }
        if (skip) { continue; } // skip only if no PDG is matching
    }

    std::unordered_map<geo::View_t, size_t> hit_count;
    for (auto const & h : p.second) { hit_count[h.View()]++; }

    size_t nviews = 0;
    for (auto const & n : hit_count) {  if (n.second > fMinHitsPerPlane) { nviews++; } }
    if (nviews >= 2)
    {
        // passed all conditions, move to filtered maps
        mapTrackIDtoHits_filtered.emplace(p);
    }
  }
  fReconstructable = mapTrackIDtoHits_filtered.size();

  std::cout << "------------ event: " << fEvent << " has reconstructable: " << fReconstructable << std::endl;
  for (auto const &p : mapTrackIDtoHits_filtered)
  {

    //Want to get the initial energy of the particle
    int this_code = pi_serv->TrackIdToParticle_P(p.first)->PdgCode();
    double this_energy = pi_serv->TrackIdToParticle_P(p.first)->E(0);
    double this_length = mapTrackIDtoLength[p.first];

    std::cout << " : id " << p.first << " size " << p.second.size() << " en: " << mapTrackIDtoHitsEnergy[p.first] << " init en: " << this_energy << " PDG: " << this_code << " Length: " << this_length << std::endl;

    fDenominatorHist->Fill(p.second.size());
    fEnergyDenominatorHist->Fill(this_energy);
    fLengthDenominatorHist->Fill(this_length);
    fDepDenominatorHist->Fill(mapTrackIDtoHitsEnergy[p.first]);
  }
  // ---------------------------------------------------------------------------------------


   // match reconstructed tracks to MC particles
  const auto trkHandle = evt.getValidHandle< std::vector<recob::Track> >(fTrackModuleLabel);
  fNRecoTracks = trkHandle->size();
  
 
  art::FindManyP< recob::Hit > hitsFromTracks(trkHandle, evt, fTrackModuleLabel);
  art::FindManyP< recob::SpacePoint > spFromHits(hitListHandle, evt, fTempTrackModuleLabel);
//  art::FindMany< anab::T0 > t0sFromTracks(trkHandle, evt, fTempTrackModuleLabel ); 


  //Make map from G4Track (IDE) ID to track index
  std::unordered_map<int, int> mapIDEtoRecoTrack;
  std::unordered_map<int, double> mapIDEtoRecoEnergy;
  std::unordered_map<int, std::vector< std::pair<int, double > > > mapIDEtoTrackContribs;

  for (size_t t = 0; t < trkHandle->size(); ++t)     // loop over tracks
  {
    // *** here you should select if the reconstructed track is interesting, eg:
    // - associated PFParticle has particlular PDG code
    // - or just the recob::Track has interesting ParticleID

    std::unordered_map<int, double> trkID_E_perPlane[3]; // map filtered MC particles to their energy contributed to track t in each plane
    std::unordered_map<int, double> trkID_E;             // map filtered MC particles to their energy contributed to track t
    double totE_anyMC_inPlane[3] = { 0, 0, 0 };
    double totE_anyMC = 0;
    const auto & hits = hitsFromTracks.at(t);
    for (const auto & h : hits)                          // loop over hits assigned to track t
    {
        size_t plane = h->WireID().Plane;
        for (auto const & ide : bt_serv->HitToTrackIDEs(clockData, h))     // loop over std::vector< sim::TrackIDE >, for a hit h
        {
            if (mapTrackIDtoHits_filtered.find(ide.trackID) != mapTrackIDtoHits_filtered.end())
            {
                trkID_E_perPlane[plane][ide.trackID] += ide.energy; // ...and sum energies contributed to this reco track
                trkID_E[ide.trackID] += ide.energy;                 // by MC particles which we considered as reconstructable
            }
            totE_anyMC_inPlane[plane] += ide.energy;
            totE_anyMC += ide.energy;
        }
    }

    // find MC particle which cotributed maximum energy to track t
    int best_id = 0;
    double max_e = 0;
    std::array< double, 3 > e_inPlane = {{ 0, 0, 0 }};
    for (auto const & entry : trkID_E)
	{
        mapIDEtoTrackContribs[entry.first].push_back(std::make_pair(t, entry.second));
        fTrkIDEs.push_back(entry.first);
        fTrkContribs.push_back(entry.second);
        if (entry.second > max_e) // find track ID corresponding to max energy
        {
            max_e = entry.second;
            best_id = entry.first;

            for (size_t i = 0; i < 3; ++i)
            {
                e_inPlane[i] = trkID_E_perPlane[i][entry.first];
            }
        }
    }    
    //'Best guess' 
    if(max_e == 0){
      mapIDEtoRecoTrack[best_id] = -1;
    }
    else{
      mapIDEtoRecoTrack[best_id] = t;       
    }
    mapIDEtoRecoEnergy[best_id] = max_e;

    // check if reco track is matching to MC particle:
    if ((max_e > 0.5 * totE_anyMC) &&                    // MC particle has more than 50% energy contribution to the track
        (max_e > 0.5 * mapTrackIDtoHitsEnergy[best_id])) // track covers more than 50% of energy deposited by MC particle in hits
    {    
        fNominatorHist->Fill(mapTrackIDtoHits_filtered[best_id].size());

        double this_length = mapTrackIDtoLength[best_id];

        //Get Energy of this track
        std::cout << " id: " << best_id << " init energy: " << pi_serv->TrackIdToParticle_P(best_id)->E(0) << " Length: " << this_length << std::endl; 
        fEnergyNominatorHist->Fill(pi_serv->TrackIdToParticle_P(best_id)->E(0));
        fLengthNominatorHist->Fill(this_length);
        fDepNominatorHist->Fill(mapTrackIDtoHitsEnergy[best_id]);
        fTrkMatched = 1; fMatched++;
    }
    else { fTrkMatched = 0; }

    fTrkPurity = max_e / totE_anyMC;
    fTrkCompletness = max_e /  mapTrackIDtoHitsEnergy[best_id];
    for (size_t i = 0; i < 3; ++i)
    {
        if (totE_anyMC_inPlane[i] > 0) { fTrkPurityPerPlane[i] = e_inPlane[i] / totE_anyMC_inPlane[i]; }
        else { fTrkPurityPerPlane[i] = 0; }

        if (mapTrackIDtoHitsEnergyPerPlane[i][best_id] > 0) { fTrkCompletnessPerPlane[i] = e_inPlane[i] / mapTrackIDtoHitsEnergyPerPlane[i][best_id]; }
        else fTrkCompletnessPerPlane[i] = 0;
    }

    fTrkPid = (*trkHandle)[t].ParticleId();
    fTrkSize = hits.size();
    fTrkLength = (*trkHandle)[t].Length();
    fTrkID = t;
    fTrkEnergy = totE_anyMC;

    //T0 stuff
//    std::vector<const anab::T0*> T0s = t0sFromTracks.at(t);       
    std::vector<const anab::T0*> T0s;    
    nT0s = T0s.size();
    std::cout << "Got " << nT0s << " T0s" <<std::endl;

    if(nT0s > 0){
      std::cout << "T0s size: " << nT0s << std::endl;          
      MCTruthT0 = T0s[0]->Time();        
    }
    else{
      std::cout << "No T0s found" << std::endl;
      MCTruthT0 = 0;
    }        
    TickT0 = MCTruthT0 / sampling_rate(clockData);


    int track_step = 0;
    for(auto const & pos : (((*trkHandle)[t].Trajectory()).Trajectory()).Positions()){
      if (track_step > 9999) break;
//      fTrkX[track_step] = pos.X(); //need to correct t0
      auto const this_wire = hits[hits.size() - (1 + track_step)]->WireID();
      double t0Correction = detProp.ConvertTicksToX( TickT0, this_wire.Plane, this_wire.TPC, this_wire.Cryostat );
      std::cout << "Correction: " << t0Correction << " " << sampling_rate(clockData) << std::endl;
      fT0Corr[track_step] = t0Correction;
      fTrkX[track_step] = (pos.X()/* - t0Correction*/);
      fTrkY[track_step] = pos.Y();
      fTrkZ[track_step] = pos.Z();  
      track_step++;             
    }
   
    fTrkNIDEs = fTrkIDEs.size();
    fTrkTree->Fill();   
    for(int track_step = 0; track_step < 10000; ++track_step){
      fTrkX[track_step] = 0;
      fTrkY[track_step] = 0;
      fTrkZ[track_step] = 0;          
    }
    
    fTrkIDEs.clear();
    fTrkContribs.clear();


    if (fTrkMatched == 1)
    {
        for (const auto & h : hits) // loop over hits assigned to matched track
        {
            const auto & sps = spFromHits.at(h.key());
            if (sps.empty()) { continue; }
            const auto & sp = *sps.front();

            std::vector< sim::IDE > ides = bt_serv->HitToAvgSimIDEs(clockData, *h);

            std::array< double, 3 > hitpos = {{0, 0, 0}};
            double hitE = 0;
            for (auto const & ide : ides)
            {
                if (ide.trackID == best_id)
                {
                    hitpos[0] += ide.x * ide.energy;
                    hitpos[1] += ide.y * ide.energy;
                    hitpos[2] += ide.z * ide.energy;
                    hitE += ide.energy;
                }
            }
            if (hitE > 0)
            {
                hitpos[0] /= hitE; hitpos[1] /= hitE; hitpos[2] /= hitE;
                double dx = 0; // sp.XYZ()[0] - hitpos[0]; // --- need t0 correction
                double dy = sp.XYZ()[1] - hitpos[1];
                double dz = sp.XYZ()[2] - hitpos[2];
                //fHitDx->Fill(dx); // --- need t0 correction
                fHitDy->Fill(dy);
                fHitDz->Fill(dz);
                fHitDist3D->Fill(sqrt(dx*dx + dy*dy + dz*dz));
            }
        }
    }
  }

  if(fEnableIDETree){
    std::cout << "IDE Tree"<< std::endl;
    //Note: mapTrackIDtoHits -> filtered?
    //Add in branch that includes the 'max_e' from above
    //as well as 'best guess' reco track
    //keep thinking about this
    for (auto const & p : mapTrackIDtoHits_filtered)
    {
      fTrueTrkID = p.first;
      fTrueTrkSize = p.second.size();
      const simb::MCParticle * part = pi_serv->TrackIdToParticle_P(fTrueTrkID);
      fTrueTrkE = part->E(0);  
      fTrueTrkLength = mapTrackIDtoLength[fTrueTrkID];   
      fTrueTrkDep = mapTrackIDtoHitsEnergy[fTrueTrkID];
      fTrueTrkPID = part->PdgCode();

      fTrueTrkOrigin = pi_serv->TrackIdToMCTruth_P(fTrueTrkID)->Origin();

      fTrueTrkRecoGuess = mapIDEtoRecoTrack[p.first];
      fTrueTrkEnergyGuess = mapIDEtoRecoEnergy[p.first];

      std::vector<std::pair<int, double>> TrackContribs;
      TrackContribs = mapIDEtoTrackContribs[p.first];
      for (auto const & tc : TrackContribs){
        fTrueTrkRecos.push_back(tc.first);
        fTrueTrkContribs.push_back(tc.second);
      }

      auto const traj = part->Trajectory();

      track_step = 0;
      for(size_t step = 0; step < traj.size(); ++step){
        //Don't know if I need this break, but just to be sure...
        if (track_step > 9999) break;

        TVector3 pos = GetPositionInTPC(clockData, *part, step);

        if( (pos != TVector3(-1,-1,-1)) && (pos != TVector3(0,0,0))){       
          fTrueTrkX[track_step] = pos.X();
          fTrueTrkY[track_step] = pos.Y();
          fTrueTrkZ[track_step] = pos.Z();  
          fTrueTrkStepE[track_step] = part->E(step); 
          if(track_step > 0){
            fTrueTrkdEdX[track_step] = GetdEdX(traj, step); 
          }
          track_step++;
        }
      }

      if(fTrueTrkX[0] < 0 && fTrueTrkX[track_step - 1] > 0){
        fTrueTrkCrossCPA = true;
      }
      else if(fTrueTrkX[0] > 0 && fTrueTrkX[track_step - 1] < 0){
        fTrueTrkCrossCPA = true;       
      }
      else{fTrueTrkCrossCPA = false;}

      fTrueTrkNRecos = fTrueTrkRecos.size();
      fTrkIDETree->Fill(); 
      for(size_t step = 0; step < 10000; ++step){
        fTrueTrkX[step] = 0;
        fTrueTrkY[step] = 0;
        fTrueTrkZ[step] = 0;
        fTrueTrkStepE[step] = 0;
        fTrueTrkdEdX[step] = 0;
      }
      fTrueTrkRecos.clear();
      fTrueTrkContribs.clear();
    }
  }
  
  if (fMatched > fReconstructable)
  {
    std::cout << "WARNING: matched more reco tracks than selected MC particles: " << fMatched << " > " << fReconstructable << std::endl;
  }
  

  std::cout << "------------ reconstructed: " << fNRecoTracks << " and then matched to MC particles: " << fMatched
    << " (reconstructable: " << fReconstructable << ")" << std::endl;
  fEventTree->Fill();
}

void pdune::RecoEff::ResetVars()
{
	fRun = 0; fEvent = 0;

	fEnGen = 0; fEkGen = 0;
	fT0 = 0;

    for (size_t i = 0; i < 3; ++i)
    {
        fTrkPurityPerPlane[i] = 0;
        fTrkCompletnessPerPlane[i] = 0;
    }

	fTrkPurity = 0; fTrkCompletness = 0;
	fTrkPid = 0; fTrkSize = 0; fTrkMatched = 0;
        fTrkLength = 0;
	fMatched = 0;
	fNRecoTracks = 0;
	fReconstructable = 0;
        fTrkID = 0;
        fTrkIDEs.clear();
        fTrkContribs.clear();
        
        if(fEnableIDETree){
          fTrueTrkE      = 0; 
          fTrueTrkID     = 0;
          fTrueTrkLength = 0; 
          fTrueTrkDep    = 0; 
          fTrueTrkPID    = 0; 
          fTrueTrkOrigin = -1;
          for(size_t i = 0; i < 10000; ++i){
            fTrueTrkX[i] = 0;
            fTrueTrkY[i] = 0;
            fTrueTrkZ[i] = 0;
          }

          fTrueTrkRecoGuess = 0;
          fTrueTrkEnergyGuess = 0;
          fTrueTrkRecos.clear();
          fTrueTrkContribs.clear();
        }


        if(fEnableParticleTree){
          fPartStatus = 0;
          fPartPID = 0;
          fPartID = 0;
          fPartLength = 0;
          fTotLength = 0;
          fOrigin = -1;
          fProcess = "";
          fEndProcess = "";
        }
}

double pdune::RecoEff::GetLengthInTPC(detinfo::DetectorClocksData const& clockData,
                                      const simb::MCParticle part) {

  unsigned int nTrajectoryPoints = part.NumberTrajectoryPoints();
  std::vector<double> TPCLengthHits(nTrajectoryPoints,0);

  bool BeenInVolume = false;
  int FirstHit = 0, LastHit = 0;
  double TPCLength = 0;

  for(unsigned int MCHit=0; MCHit < nTrajectoryPoints; ++MCHit) {
    if (MCHit!=0) TPCLengthHits[MCHit] = pow ( pow( (part.Vx(MCHit-1)-part.Vx(MCHit)),2)
                                             + pow( (part.Vy(MCHit-1)-part.Vy(MCHit)),2)
    				         + pow( (part.Vz(MCHit-1)-part.Vz(MCHit)),2)
    				         , 0.5 );
    // check if in TPC
    auto const position = geo::vect::toPoint(part.Position(MCHit).Vect());

    geo::TPCID tpcid = geom->FindTPCAtPosition(position);
    if(tpcid.isValid){
      double XPlanePosition = wireReadout.Plane({tpcid, 0}).GetCenter().X();
      double DriftTimeCorrection = fabs( position.X() - XPlanePosition )/ XDriftVelocity;
      double TimeAtPlane = part.T() + DriftTimeCorrection;
      if( TimeAtPlane < trigger_offset(clockData) || TimeAtPlane > trigger_offset(clockData) + WindowSize){
        continue;}
      //Good hit in TPC
      LastHit = MCHit;
      if( !BeenInVolume ){
        BeenInVolume = true;
        FirstHit = MCHit;
      } 
    }
  }
  for (int Hit = FirstHit+1; Hit <= LastHit; ++Hit ) {
    TPCLength += TPCLengthHits[Hit];
  }
  return TPCLength;
}

TVector3 pdune::RecoEff::GetPositionInTPC(detinfo::DetectorClocksData const& clockData,
                                          const simb::MCParticle part, int MCHit) {
  auto const position = geo::vect::toPoint(part.Position(MCHit).Vect());

  // check if in TPC
  geo::TPCID tpcid = geom->FindTPCAtPosition(position);
  if(tpcid.isValid){
    double XPlanePosition = wireReadout.Plane({tpcid, 0}).GetCenter().X();
    double DriftTimeCorrection = fabs( position.X() - XPlanePosition )/ XDriftVelocity;
    double TimeAtPlane = part.T() + DriftTimeCorrection;
    if( TimeAtPlane < trigger_offset(clockData) || TimeAtPlane > trigger_offset(clockData) + WindowSize){
      return TVector3(-1,-1,-1);}
  }     
  else{
    return TVector3(-1,-1,-1);
  }
  
  return geo::vect::convertTo<TVector3>(position);
}

double pdune::RecoEff::GetdEdX(const simb::MCTrajectory traj, size_t this_step){
 size_t prev_step = this_step - 1;
 size_t next_step = this_step + 1;

 double dE_prev = fabs( traj.E(this_step) - traj.E(prev_step) );
 double dE_next = fabs( traj.E(this_step) - traj.E(next_step) );

 double dX_prev = sqrt( pow( (traj.X(this_step) - traj.X(prev_step)), 2 ) 
                 + pow( (traj.X(this_step) - traj.X(prev_step)), 2 )
                 + pow( (traj.X(this_step) - traj.X(prev_step)), 2 ));
                
 double dX_next = sqrt( pow( (traj.X(this_step) - traj.X(next_step)), 2 ) 
                 + pow( (traj.X(this_step) - traj.X(next_step)), 2 )
                 + pow( (traj.X(this_step) - traj.X(next_step)), 2 ));
 return (dE_prev + dE_next)/(dX_prev + dX_next);
}

DEFINE_ART_MODULE(pdune::RecoEff)
