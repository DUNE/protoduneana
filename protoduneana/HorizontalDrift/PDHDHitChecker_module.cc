////////////////////////////////////////////////////////////////////////
// Class:       PDHDHitChecker
// Plugin Type: producer (Unknown Unknown)
// File:        PDHDHitChecker_module.cc
//
// Generated at Fri Aug 18 12:13:50 2023 by Jacob Calcutt using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/make_tool.h"
#include "TTree.h"
#include "art_root_io/TFileService.h"

#include "lardataobj/RecoBase/SpacePoint.h"

#include <memory>
namespace pdhd {

class PDHDHitChecker;

class PDHDHitChecker : public art::EDAnalyzer {
public:
  explicit PDHDHitChecker(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDHDHitChecker(PDHDHitChecker const&) = delete;
  PDHDHitChecker(PDHDHitChecker&&) = delete;
  PDHDHitChecker& operator=(PDHDHitChecker const&) = delete;
  PDHDHitChecker& operator=(PDHDHitChecker&&) = delete;

  // Required functions.
  void analyze(art::Event const &event) override;
  void beginJob() override;

private:
  void reset();
  bool good_sp(const recob::SpacePoint & sp);
  std::string fSpacePointTag;
  double fXLow, fYLow, fZHigh;
  TTree * fTree;
  int run, eventnum, subrun;

  std::vector<double> sp_xs, sp_ys, sp_zs;
};
}

void pdhd::PDHDHitChecker::beginJob() {


  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("hitcheck", "");
  fTree->Branch("event", &eventnum);
  fTree->Branch("run", &run);
  fTree->Branch("subrun", &subrun);
  fTree->Branch("sp_xs", &sp_xs);
  fTree->Branch("sp_ys", &sp_ys);
  fTree->Branch("sp_zs", &sp_zs);

}

void pdhd::PDHDHitChecker::reset() {
  eventnum = -999;
  subrun = -999;
  run = -999;
  sp_xs.clear();
  sp_ys.clear();
  sp_zs.clear();
}

pdhd::PDHDHitChecker::PDHDHitChecker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSpacePointTag(p.get<std::string>("SpacePointTag")),
    fXLow(p.get<double>("XLow", 0.)),
    fYLow(p.get<double>("YLow", 0.)),
    fZHigh(p.get<double>("ZHigh", 220.)) {
}

bool pdhd::PDHDHitChecker::good_sp(const recob::SpacePoint & sp) {
  return (sp.XYZ()[0] > fXLow && sp.XYZ()[1] > fYLow && sp.XYZ()[2] < fZHigh);
}

void pdhd::PDHDHitChecker::analyze(art::Event const& event) {
  //auto allHits = event.getValidHandle<std::vector<recob::Hit> >(fHitTag);
  reset();
  eventnum = event.id().event();
  subrun = event.subRun();
  run = event.run();
  auto space_pt_handle = event.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointTag);

  for (auto & space_point : (*space_pt_handle)) {
    if (!good_sp(space_point)) continue;
    //std::cout << space_point.XYZ()[0] << std::endl;
    sp_xs.push_back(space_point.XYZ()[0]);
    sp_ys.push_back(space_point.XYZ()[1]);
    sp_zs.push_back(space_point.XYZ()[2]);
  }
  fTree->Fill();
}

DEFINE_ART_MODULE(pdhd::PDHDHitChecker)
