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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/make_tool.h"
#include "TTree.h"
#include "art_root_io/TFileService.h"
#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

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
  bool good_hit(const recob::Hit & hit);
  std::string fSpacePointTag;
  std::string fHitTag;
  double fXLow, fYLow, fZHigh;
  TTree * fTree;
  std::vector<size_t> fValidTPCs;
  int run, eventnum, subrun;

  std::vector<double> sp_xs, sp_ys, sp_zs;
  std::vector<double> sp_inplane_zs, sp_inplane_ys, sp_inplane_time;
  std::vector<int> sp_inplane_tpc;
  std::vector<double> hit_plane, hit_time, hit_integral;
  std::vector<int> hit_channel, hit_tpc, hit_wire;

  art::ServiceHandle <geo::Geometry> fGeometryService;
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
  fTree->Branch("sp_inplane_zs", &sp_inplane_zs);
  fTree->Branch("sp_inplane_ys", &sp_inplane_ys);
  fTree->Branch("sp_inplane_tpc", &sp_inplane_tpc);
  fTree->Branch("sp_inplane_time", &sp_inplane_time);
  fTree->Branch("hit_channel", &hit_channel);
  fTree->Branch("hit_wire", &hit_wire);
  fTree->Branch("hit_plane", &hit_plane);
  fTree->Branch("hit_time", &hit_time);
  fTree->Branch("hit_tpc", &hit_tpc);
  fTree->Branch("hit_integral", &hit_integral);
}

void pdhd::PDHDHitChecker::reset() {
  eventnum = -999;
  subrun = -999;
  run = -999;
  sp_xs.clear();
  sp_ys.clear();
  sp_zs.clear();

  sp_inplane_ys.clear();
  sp_inplane_zs.clear();
  sp_inplane_tpc.clear();
  sp_inplane_time.clear();

  hit_channel.clear();
  hit_wire.clear();
  hit_integral.clear();
  hit_plane.clear();
  hit_time.clear();
  hit_tpc.clear();
}

pdhd::PDHDHitChecker::PDHDHitChecker(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fSpacePointTag(p.get<std::string>("SpacePointTag")),
    fHitTag(p.get<std::string>("HitTag")),
    fXLow(p.get<double>("XLow", -999.)),
    fYLow(p.get<double>("YLow", -999.)),
    fZHigh(p.get<double>("ZHigh", 999.)),
    fValidTPCs(p.get<std::vector<size_t>>("ValidTPCs", {2,6})) {
}

bool pdhd::PDHDHitChecker::good_sp(const recob::SpacePoint & sp) {
  return (sp.XYZ()[0] > fXLow && sp.XYZ()[1] > fYLow && sp.XYZ()[2] < fZHigh);
}

bool pdhd::PDHDHitChecker::good_hit(const recob::Hit & hit) {
  return (hit.WireID().TPC == 2 || hit.WireID().TPC == 6);
}

void pdhd::PDHDHitChecker::analyze(art::Event const& event) {
  reset();
  eventnum = event.id().event();
  subrun = event.subRun();
  run = event.run();
  const auto space_pt_handle = event.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointTag);
  const auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHitTag);

  //Check hit tag or sp tag?
  art::FindManyP<recob::Hit> associated_hits(space_pt_handle, event, fHitTag);
  //size_t iSP = 0; 
  //for (const auto & space_point : (*space_pt_handle)) {
  for (size_t iSP = 0; iSP < space_pt_handle->size(); ++iSP) {
    //if (!good_sp(space_point)) {
    //  ++iSP;
    //  continue;
    //}
    //std::cout << space_point.XYZ()[0] << std::endl;
    //sp_xs.push_back(space_point.XYZ()[0]);
    //sp_ys.push_back(space_point.XYZ()[1]);
    //sp_zs.push_back(space_point.XYZ()[2]);
    const auto & the_hits = associated_hits.at(iSP);
    //std::cout << "got " << the_hits.size() << " from SP " << iSP << std::endl;

    //Get the Wire IDs from the hit
    std::vector<geo::WireID> wires;
    bool any_excluded = false;
    for (const auto & hit : the_hits) {
      wires.push_back((*hit).WireID());
      if (std::find(fValidTPCs.begin(), fValidTPCs.end(), (*hit).WireID().TPC) == fValidTPCs.end())
        any_excluded = true;
    }
    if (any_excluded) continue;

    //Only have 1 intersection
    if (the_hits.size() == 2) {
      double y = 0., z = 0.;
      geo::Point_t point0;
      bool intersect = fGeometryService->WireIDsIntersect(
          wires[0], wires[1], point0);
      if (!intersect) {
        std::cout << "Warning no intersecion for single pair" << std::endl;
        //for (const auto & hit : the_hits) {
        //  std::cout << (*hit).View() << " ";
        //}

        continue;
      }
      y = point0.y();
      z = point0.z();
      //std::cout << "intersect at " << y << " " << z << std::endl;
      sp_inplane_ys.push_back(y);
      sp_inplane_zs.push_back(z);
    }
    else if (the_hits.size() == 3) {
      geo::Point_t point0;
      bool first_intersect = fGeometryService->WireIDsIntersect(
          wires[0], wires[1], point0);
      if (first_intersect) {
        sp_inplane_ys.push_back(point0.y());
        sp_inplane_zs.push_back(point0.z());
      }

      geo::Point_t point1;
      bool second_intersect = fGeometryService->WireIDsIntersect(
          wires[1], wires[2], point1);
      if (!first_intersect && second_intersect) {
        sp_inplane_ys.push_back(point1.y());
        sp_inplane_zs.push_back(point1.z());
      }
      else {
        std::cout << "Warning. Neither intersect" << std::endl;
      }
      //std::cout << "second pair intersect at " << y1 << " " << z1 << std::endl;
    }
    else {
      throw cet::exception("PDHDHitChecker_module.cc") <<
        "Got " << the_hits.size() << " hits from SP " << iSP;
    }
    //Just use first hit for tpc and time
    sp_inplane_tpc.push_back(the_hits[0]->WireID().TPC);
    sp_inplane_time.push_back(the_hits[0]->PeakTime());
    //++iSP;
  }

  /*for (const auto & hit : (*hit_handle)) {
    if (!good_hit(hit)) continue;
    hit_channel.push_back(hit.Channel());
    hit_wire.push_back(hit.WireID().Wire);
    hit_plane.push_back(hit.View());
    hit_time.push_back(hit.PeakTime());
    hit_tpc.push_back(hit.WireID().TPC);
    hit_integral.push_back(hit.Integral());
  }*/
  fTree->Fill();
}

DEFINE_ART_MODULE(pdhd::PDHDHitChecker)
