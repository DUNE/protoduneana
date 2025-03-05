////////////////////////////////////////////////////////////////////////
// Class:       ClusterCounter3
// Module Type: analyzer
// File:        ClusterCounter3_module.cc
//
// Access clusters, find assigned hits, fill ROOT tree with.
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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TTree.h"

namespace tutorial {

class ClusterCounter3 : public art::EDAnalyzer {
public:
  explicit ClusterCounter3(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  ClusterCounter3(ClusterCounter3 const &) = delete;
  ClusterCounter3(ClusterCounter3 &&) = delete;
  ClusterCounter3 & operator = (ClusterCounter3 const &) = delete;
  ClusterCounter3 & operator = (ClusterCounter3 &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.  
  void beginJob() override;
  void endJob() override;

private:
  float sumAdc(const std::vector< art::Ptr<recob::Hit> > & hits) const;

  size_t fEvNumber;

  TTree *fEventTree;
  size_t fNClusters;

  TTree *fClusterTree;
  size_t fNHits;
  float fAdcSum;

  // ******* fcl parameters *******
  art::InputTag fClusterModuleLabel;
  size_t fMinSize;
};

ClusterCounter3::ClusterCounter3(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
    fClusterModuleLabel = p.get< std::string >("ClusterModuleLabel");
    fMinSize = p.get< size_t >("MinSize");
}

void ClusterCounter3::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs; // TTree's are created in the memory managed by ROOT (you don't delete them)

    fEventTree = tfs->make<TTree>("EventTree", "event by event info");
    fEventTree->Branch("event", &fEvNumber, "fEvNumber/I");
    fEventTree->Branch("nclusters", &fNClusters, "fNClusters/I");

    fClusterTree = tfs->make<TTree>("ClusterTree", "cluster by cluster info");
    fClusterTree->Branch("event", &fEvNumber, "fEvNumber/I");
    fClusterTree->Branch("nhits", &fNHits, "fNHits/I");
    fClusterTree->Branch("adcsum", &fAdcSum, "fAdcSum/F");
}

void ClusterCounter3::analyze(art::Event const & evt)
{
    fEvNumber = evt.id().event();
    mf::LogVerbatim("ClusterCounter3") << "ClusterCounter3 module on event #" << fEvNumber;

    // use auto to make the line shorter when you remember all art types,
    // the full type of the handle is: art::ValidHandle< std::vector<recob::Cluster> >
    auto cluHandle = evt.getValidHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);

    // this will let us look for all hits associated to selected cluster
    art::FindManyP< recob::Hit > hitsFromClusters(cluHandle, evt, fClusterModuleLabel);

    fNClusters = 0;

    // and here we go by index:
    for (size_t i = 0; i < cluHandle->size(); ++i)
    {
        fNHits = cluHandle->at(i).NHits();
        fAdcSum = sumAdc(hitsFromClusters.at(i));
        fClusterTree->Fill();

         mf::LogVerbatim("ClusterCounter3")
            << "NHits() = " << fNHits << ", assn size = " << hitsFromClusters.at(i).size()
            << " SummedADC() = " << cluHandle->at(i).SummedADC() << ", sum hits adc = " << fAdcSum;

        if (fNHits >= fMinSize) { ++fNClusters; }
    }
    fEventTree->Fill();
}

float ClusterCounter3::sumAdc(const std::vector< art::Ptr<recob::Hit> > & hits) const
{
    float sum = 0;
    for (auto const & h : hits)
    {
        sum += h->ROISummedADC();
    }
    return sum;
}

void ClusterCounter3::endJob()
{
    mf::LogVerbatim("ClusterCounter3") << "ClusterCounter3 finished job";
}

} // tutorial namespace

DEFINE_ART_MODULE(tutorial::ClusterCounter3)
