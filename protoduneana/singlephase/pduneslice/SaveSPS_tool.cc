#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/ToolMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
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
//#include "hep_hpc/hdf5/make_column.hpp"

using namespace hep_hpc::hdf5;

namespace dnn {

  class SaveSPS : public ImageMaker {
  public:
    explicit SaveSPS(fhicl::ParameterSet const& ps);

    void saveImage(art::Event const& e, hep_hpc::hdf5::File& hdffile) override;

  private:
    art::InputTag fSpacePointModuleLabel;
  };

  SaveSPS::SaveSPS(fhicl::ParameterSet const& ps)
    : fSpacePointModuleLabel{ps.get<art::InputTag>("SpacePointModuleLabel")}
  {}

  void SaveSPS::saveImage(art::Event const& e, hep_hpc::hdf5::File& hdffile)
  {

    art::ServiceHandle<geo::Geometry> geom;

//    static Ntuple<unsigned int, unsigned int, unsigned int, double> evtids(
//      hdffile, "evtids", {"run", "subrun", "event", "evttime"});

    static auto graph =
      make_ntuple({hdffile, "graph", 1000},
                  make_scalar_column<int32_t>("trackid"),
                  make_scalar_column<int32_t>("pdg"),
                  make_column<float, 1>(
                    "xyzq",                                // 1 means each element is a 1-d array
                    4,                             // extent of each array dimension
                    1024/ (4 * sizeof(float)), // chunk size
                    {PropertyList{H5P_DATASET_CREATE}(&H5Pset_shuffle)(&H5Pset_deflate, 6u)}));

//    // Get event information
//    double evttime;
//
//    art::Timestamp ts = e.time();
//    if (ts.timeHigh() == 0) {
//      TTimeStamp tts(ts.timeLow());
//      evttime = tts.AsDouble();
//    }
//    else {
//      TTimeStamp tts(ts.timeHigh(), ts.timeLow());
//      evttime = tts.AsDouble();
//    }

    // Get space point information
    auto const & sps = e.getProduct<std::vector<recob::SpacePoint>>(fSpacePointModuleLabel);
    auto const & pcs = e.getProduct<std::vector<recob::PointCharge>>(fSpacePointModuleLabel);

    auto spsHandle = e.getHandle< std::vector<recob::SpacePoint> >(fSpacePointModuleLabel);
    art::FindManyP<recob::Hit> fmhsp(spsHandle, e, fSpacePointModuleLabel);

    //Services
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

    for (size_t i = 0; i<sps.size(); ++i){
      float xyzq[4] = {0.};
      xyzq[0] = sps[i].XYZ()[0];
      xyzq[1] = sps[i].XYZ()[1];
      xyzq[2] = sps[i].XYZ()[2];
      xyzq[3] = pcs[i].charge();
      int trackid = 0;
      int pdg = 0;
      auto const& hits = fmhsp.at(i);
      if (!e.isRealData()){
        std::map<int,double> trkide;
        for (auto const & hit : hits){
          std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(clockData, hit);
          for(size_t j = 0; j < TrackIDs.size(); ++j){
            trkide[TrackIDs[j].trackID] += TrackIDs[j].energy;
          }
        }
        // Work out which IDE despoited the most charge in the hit if there was more than one.
        double maxe = -1;
        //double tote = 0;
        for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
          //tote += ii->second;
          if ((ii->second)>maxe){
            maxe = ii->second;
            trackid = ii->first;
          }
        }
        // Now have trackID, so get PdG code and T0 etc.
        const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(trackid);
        if (particle){
          pdg = particle->PdgCode();
        }
      }
      std::cout<<trackid<<" "<<pdg<<" "<<xyzq[0]<<" "<<xyzq[1]<<" "<<xyzq[2]<<" "<<xyzq[3]<<std::endl;
      graph.insert(trackid, pdg, &xyzq[0]);
    }
  
    return;
  }
}

DEFINE_ART_CLASS_TOOL(dnn::SaveSPS)
