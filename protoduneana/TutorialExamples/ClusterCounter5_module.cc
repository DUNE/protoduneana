////////////////////////////////////////////////////////////////////////
// Class:       ClusterCounter
// Module Type: analyzer
// File:        ClusterCounter_module.cc
//
// Same as #4, but use FHiCL parameters validation.
// Robert Sulej
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "TTree.h"

namespace tutorial {

class ClusterCounter : public art::EDAnalyzer {
public:

  struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> ClusterModuleLabel {
          Name("ClusterModuleLabel"),
          Comment("tag of cluster producer")
      };

      fhicl::Atom<size_t> MinSize {
          Name("MinSize"),
          Comment("minimum size of clusters")
      };
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit ClusterCounter(Parameters const& config);

  // Plugins should not be copied or assigned.
  ClusterCounter(ClusterCounter const &) = delete;
  ClusterCounter(ClusterCounter &&) = delete;
  ClusterCounter & operator = (ClusterCounter const &) = delete;
  ClusterCounter & operator = (ClusterCounter &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.  
  void beginJob() override;
  void endJob() override;

private:
  float sumAdc(const std::vector< art::Ptr<recob::Hit> > & hits) const;

  const simb::MCParticle* getTruthParticle(
    detinfo::DetectorClocksData const& clockData,
    const std::vector< art::Ptr<recob::Hit> > & hits,
    float & fraction, bool & foundEmParent) const;

  size_t fEvNumber;

  TTree *fEventTree;
  size_t fNClusters;

  TTree *fClusterTree;
  size_t fNHits;
  float fAdcSum;
  float fClean;

  // ******* fcl parameters *******
  art::InputTag fClusterModuleLabel;
  size_t fMinSize;
};

ClusterCounter::ClusterCounter(Parameters const& config) : art::EDAnalyzer(config),
    fClusterModuleLabel(config().ClusterModuleLabel()),
	fMinSize(config().MinSize())
{
}

void ClusterCounter::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs; // TTree's are created in the memory managed by ROOT (you don't delete them)

    fEventTree = tfs->make<TTree>("EventTree", "event by event info");
    fEventTree->Branch("event", &fEvNumber, "fEvNumber/I");
    fEventTree->Branch("nclusters", &fNClusters, "fNClusters/I");

    fClusterTree = tfs->make<TTree>("ClusterTree", "cluster by cluster info");
    fClusterTree->Branch("event", &fEvNumber, "fEvNumber/I");
    fClusterTree->Branch("nhits", &fNHits, "fNHits/I");
    fClusterTree->Branch("adcsum", &fAdcSum, "fAdcSum/F");
    fClusterTree->Branch("clean", &fClean, "fClean/F");
}

void ClusterCounter::analyze(art::Event const & evt)
{
    fEvNumber = evt.id().event();
    mf::LogVerbatim("ClusterCounter") << "ClusterCounter module on event #" << fEvNumber;

    // use auto to make the line shorter when you remember all art types,
    // the full type of the handle is: art::ValidHandle< std::vector<recob::Cluster> >
    auto cluHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);

    // this will let us look for all hits associated to selected cluster
    art::FindManyP< recob::Hit > hitsFromClusters(cluHandle, evt, fClusterModuleLabel);

    fNClusters = 0;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    // and here we go by index:
    for (size_t i = 0; i < cluHandle->size(); ++i)
    {
        fNHits = cluHandle->at(i).NHits();
        fAdcSum = sumAdc(hitsFromClusters.at(i));

        bool isEM = false;
        const simb::MCParticle* p = getTruthParticle(clockData, hitsFromClusters.at(i), fClean, isEM);
        if (p)
        {
            if (isEM) { mf::LogVerbatim("ClusterCounter") << "matched mother particle PDG: " << p->PdgCode(); }
            else { mf::LogVerbatim("ClusterCounter") << "matched particle PDG: " << p->PdgCode(); }
        }
        else { mf::LogWarning("ClusterCounter") << "No matcched particle??"; }

        fClusterTree->Fill();

        mf::LogVerbatim("ClusterCounter")
            << "NHits() = " << fNHits << ", assn size = " << hitsFromClusters.at(i).size()
            << " SummedADC() = " << cluHandle->at(i).SummedADC() << ", sum hits adc = " << fAdcSum;

        if (fNHits >= fMinSize) { ++fNClusters; }
    }
    fEventTree->Fill();
}

const simb::MCParticle* ClusterCounter::getTruthParticle(detinfo::DetectorClocksData const& clockData,
                                                         const std::vector< art::Ptr<recob::Hit> > & hits,
    float & fraction, bool & foundEmParent) const
{
  const simb::MCParticle* mcParticle =0;
    fraction = 0;
    foundEmParent = false;

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::unordered_map<int, double> trkIDE;
    for (auto const & h : hits)
    {
        for (auto const & ide : bt_serv->HitToTrackIDEs(clockData, h)) // loop over std::vector<sim::TrackIDE>
        {
            trkIDE[ide.trackID] += ide.energy; // sum energy contribution by each track ID
        }
    }

    int best_id = 0;
    double tot_e = 0, max_e = 0;
    for (auto const & contrib : trkIDE)
    {
        tot_e += contrib.second;     // sum total energy in these hits
        if (contrib.second > max_e)  // find track ID corresponding to max energy
        {
            max_e = contrib.second;
            best_id = contrib.first;
        }
    }

    if ((max_e > 0) && (tot_e > 0)) // ok, found something reasonable
    {
        if (best_id < 0)            // NOTE: negative ID means this is EM activity
        {                           // caused by track with the same but positive ID
            best_id = -best_id;     // --> we'll find mother MCParticle of these hits
            foundEmParent = true;
        }
        mcParticle = pi_serv->TrackIdToParticle_P(best_id); // MCParticle corresponding to track ID
        fraction = max_e / tot_e;
    }
    else { mf::LogWarning("ClusterCounter") << "No energy deposits??"; }

    return mcParticle;
}

float ClusterCounter::sumAdc(const std::vector< art::Ptr<recob::Hit> > & hits) const
{
    float sum = 0;
    for (auto const & h : hits)
    {
        sum += h->ROISummedADC();
    }
    return sum;
}

void ClusterCounter::endJob()
{
    mf::LogVerbatim("ClusterCounter") << "ClusterCounter finished job";
}

} // tutorial namespace

DEFINE_ART_MODULE(tutorial::ClusterCounter)
