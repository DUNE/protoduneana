#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PointCharge.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larrecodnn/ImageMaker/ImageMaker.h"
#include "TTimeStamp.h"
#include "hep_hpc/hdf5/Ntuple.hpp"
#include "hep_hpc/hdf5/make_ntuple.hpp"

using namespace hep_hpc::hdf5;

namespace dnn {

  class Hit3D : public ImageMaker {
  public:
    explicit Hit3D(fhicl::ParameterSet const& ps);

    void saveImage(art::Event const& e, hep_hpc::hdf5::File& hdffile) override;

  private:
    art::InputTag fHitModuleLabel;
  };

  Hit3D::Hit3D(fhicl::ParameterSet const& ps)
    : fHitModuleLabel{ps.get<art::InputTag>("HitModuleLabel")}
  {}

  void Hit3D::saveImage(art::Event const& e, hep_hpc::hdf5::File& hdffile)
  {

    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

//    static Ntuple<unsigned int, unsigned int, unsigned int, double> evtids(
//      hdffile, "evtids", {"run", "subrun", "event", "evttime"});

    static auto graph =
      make_ntuple({hdffile, "graph", 1000},
                  make_scalar_column<int32_t>("run"),
                  make_scalar_column<int32_t>("subrun"),
                  make_scalar_column<int32_t>("event"),
                  make_scalar_column<int32_t>("tpc"),
                  make_scalar_column<int32_t>("plane"),
                  make_scalar_column<int32_t>("wire"),
                  make_scalar_column<int32_t>("channel"),
                  make_scalar_column<float>("charge"),
                  make_scalar_column<float>("peakt"),
                  make_scalar_column<float>("x"),
                  make_scalar_column<float>("y"),
                  make_scalar_column<float>("z"),
                  make_scalar_column<float>("energy")
                  );

    auto const & hits = e.getProduct<std::vector<recob::Hit>>(fHitModuleLabel);

    //Services
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    for (size_t i = 0; i<hits.size(); ++i){
      int tpc = hits[i].WireID().TPC;
      int plane = hits[i].WireID().Plane;
      int wire = hits[i].WireID().Wire;
      int channel = wireReadout.PlaneWireToChannel(hits[i].WireID());
      float charge = hits[i].Integral();
      float peakt = hits[i].PeakTime();
      float x = 0;
      float y = 0;
      float z = 0;
      float energy = 0;
      if (!e.isRealData()){
        auto const & ides = bt_serv->HitToAvgSimIDEs(clockData, hits[i]);
        for(size_t j = 0; j < ides.size(); ++j){
          if (ides[j].energy > energy){
            energy = ides[j].energy;
            x = ides[j].x;
            y = ides[j].y;
            z = ides[j].z;
          }
        }
      }
      //std::cout<<tpc<<" "<<plane<<" "<<wire<<" "<<channel<<" "<<charge<<" "<<peakt<<" "<<x<<" "<<y<<" "<<z<<" "<<energy<<std::endl;
      graph.insert(e.run(), e.subRun(), e.id().event(), tpc, plane, wire, channel, charge, peakt, x, y, z, energy);
    }
  
    return;
  }
}

DEFINE_ART_CLASS_TOOL(dnn::Hit3D)
