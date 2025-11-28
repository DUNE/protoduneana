////////////////////////////////////////////////////////////////////////
// Class:       MichelTimingMcDeconHD
// Plugin Type: analyzer (art v3_05_01)
// File:        MichelTimingMcDeconHD_module.cc
//
// Generated at Feb 14, 2025 by Shuaixiang (Shu) Zhang
// Based on MichelTiming_module.cc developed by Tingjun
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/OpWaveform.h"

#include "TH1F.h"
#include "TTree.h"

#include <vector>

constexpr int kMaxWF = 4000;

using namespace std;

// Shu: Initial namespace is pdsp; pdhd indeed exists;
// Concern here is if pdhd includes enough modules available for this
// analysis---
namespace pdhd {
class MichelTimingMcDeconHD;
}

class pdhd::MichelTimingMcDeconHD : public art::EDAnalyzer {
public:
  explicit MichelTimingMcDeconHD(fhicl::ParameterSet const &p);

  // Plugins should not be copied or assigned.
  MichelTimingMcDeconHD(MichelTimingMcDeconHD const &) = delete;
  MichelTimingMcDeconHD(MichelTimingMcDeconHD &&) = delete;
  MichelTimingMcDeconHD &operator=(MichelTimingMcDeconHD const &) = delete;
  MichelTimingMcDeconHD &operator=(MichelTimingMcDeconHD &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  // art::InputTag fOpDetModuleLabel;

  TTree *anatree;
  int run;
  int event;
  vector<float> pandorat0;
  vector<int> trkid;
  vector<float> vtxx, vtxy, vtxz;
  vector<float> endx, endy, endz;
  vector<float> michelscore;
  vector<int> michelhits;
  vector<short> pdchannel;
  vector<float> pdt0;

  int nWF;

  // Test. 20241227---
  int nWFLong;
  int nWFShort;
  int nWFTotal;

  //  int waveform[kMaxWF][1024];
  float waveform[kMaxWF][1024];

  art::ServiceHandle<art::TFileService> tfs;
};

pdhd::MichelTimingMcDeconHD::MichelTimingMcDeconHD(fhicl::ParameterSet const &p)
    : EDAnalyzer{p} {}

void pdhd::MichelTimingMcDeconHD::analyze(art::Event const &e) {
  run = e.run();
  event = e.id().event();

  // Shu: Based on ChatGPT's explanation, this is related to PDS wvf...
  auto wfListHandle =
      e.getHandle<std::vector<recob::OpWaveform>>("opdec::Reco");
  pdchannel.clear();
  pdt0.clear();

  nWF = 0;

  nWFLong = 0;
  nWFShort = 0;
  nWFTotal = 0;

  for (int i = 0; i < kMaxWF; ++i) {
    for (int j = 0; j < 1024; ++j) {
      waveform[i][j] = 0;
    }
  }

  // Shu: All wvfs of internal trigger are considered for pdt0 calculation---
  for (const auto &wf : *wfListHandle) { // Use reference instead of copying
    if (wf.Channel() < 120) {
      pdchannel.push_back(wf.Channel());
      pdt0.push_back(wf.TimeStamp());
    }

    if (nWF < kMaxWF) {
      if (wf.Channel() < 120 && wf.Signal().size() <= 1024) {

        //        std::copy(wf.Signal().begin(), wf.Signal().end(),
        //        waveform[nWF]);//Replace Waveform() by Signal()
        std::copy_n(wf.Signal().begin(), 1024, waveform[nWF]);
        ++nWF;

        if (wf.Signal().size() > 1024) {
          std::cout << "Length of long Wvf (opch<120): " << wf.Signal().size()
                    << "; opch: " << wf.Channel() << "\n"
                    << std::endl;
          ++nWFLong;
        }

        if (wf.Signal().size() < 1024) {
          std::cout << "Length of normal Wvf (opch<120): " << wf.Signal().size()
                    << "; opch: " << wf.Channel() << "\n"
                    << std::endl;
          ++nWFShort;
        }
      }

      if (wf.Signal().size() > 1024) {
        std::cout << "Long Wvf: " << wf.Signal().size()
                  << "; opch: " << wf.Channel() << "\n"
                  << std::endl;
      }

      ++nWFTotal;
    }
  }

  std::cout << "\n\nOpch<120: Num of Long(>1024) wvf: " << nWFLong << "\n"
            << std::endl;
  std::cout << "Opch<120: Num of Short(<1024) wvf: " << nWFShort << "\n"
            << std::endl;
  std::cout << "Opch<120: Num of Normal(1024) & Short wvf: " << nWF << "\n"
            << std::endl;

  std::cout << "\nAll OpChs: Num of All([0, inf)) wvf: " << nWFTotal << "\n\n"
            << std::endl;

  pandorat0.clear();
  trkid.clear();
  vtxx.clear();
  vtxy.clear();
  vtxz.clear();
  endx.clear();
  endy.clear();
  endz.clear();
  michelscore.clear();
  michelhits.clear();

  // Shu: e (event) must contain both PDS and TPC info!
  // Shu: Here it seems that it retrieves T0 from pandora (TPC)
  // There is one T0 for each track---
  auto t0ListHandle = e.getHandle<std::vector<anab::T0>>("pandora");
  //  auto t0ListHandle = e.getHandle< std::vector<anab::T0>
  //  >("pandora::pdhdkeepupstage2");
  if (!t0ListHandle)
    return;

  // Shu: According to Tingjun, We should find the PFParticle (particle flow
  // ...) at first, Shu: then select the T0
  auto pfListHandle = e.getHandle<std::vector<recob::PFParticle>>("pandora");
  //  auto pfListHandle = e.getHandle< std::vector<recob::PFParticle>
  //  >("pandora::pdhdkeepupstage2");
  if (!pfListHandle)
    return;

  auto trkListHandle = e.getHandle<std::vector<recob::Track>>("pandoraTrack");
  //  auto trkListHandle = e.getHandle< std::vector<recob::Track>
  //  >("pandoraTrack::pdhdkeepupstage2");
  if (!trkListHandle)
    return;

  std::vector<art::Ptr<recob::Track>> trkList;
  if (trkListHandle) {
    art::fill_ptr_vector(
        trkList,
        trkListHandle); // Shu: Now trkList is filled with trkListHandle---
  }

  // Shu: link each element in t0ListHandle (T0!) to a single recob::PFParticle
  art::FindOneP<recob::PFParticle> assn3(t0ListHandle, e, "pandora");

  // Shu: Link PFParticle with pandora track
  art::FindOneP<recob::Track> assn4(pfListHandle, e, "pandoraTrack");

  std::vector<art::Ptr<recob::Hit>> hitList;
  auto hitListHandle = e.getHandle<std::vector<recob::Hit>>("hitpdune");
  //  auto hitListHandle = e.getHandle < std::vector < recob::Hit >
  //  >("hitpdune::pdhdkeepupstage2");

  if (hitListHandle) {
    art::fill_ptr_vector(
        hitList, hitListHandle); // Shu: hitList contains pointers to all the
                                 // recob::Hit objects from the "hitpdune"---
  }

  // Shu: The fmthm object allows you to retrieve all recob::Hit objects
  // associated with each recob::Track in trkListHandle. Shu: For example, if
  // you want to get all hits associated with a specific track, you could call
  // fmthm.at(trackIndex) to retrieve a vector of art::Ptr<recob::Hit> pointers
  // linked to that track.
  //  Get track-hit association
  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(
      trkListHandle, e, "pandoraTrack"); // to associate track and hits

  // Shu: The thass object allows you to retrieve all recob::Track objects
  // associated with each recob::Hit in hitListHandle. Shu: For example, if you
  // wanted to find out which tracks are associated with a specific hit, you
  // could call thass.at(hitIndex) to retrieve a vector of
  // art::Ptr<recob::Track> pointers linked to that hit. Shu: One hit can be
  // related to several tracks; For example, the joint hit between mu and e
  // track
  art::FindManyP<recob::Track> thass(
      hitListHandle, e, "pandoraTrack"); // to associate hit just trying

  // Shu: The hitResults object will be used to apply the MVA model to
  // recob::Hit objects in the current event, allowing for the classification or
  // scoring of hits based on the model's predictions. Shu: After initializing
  // hitResults, you can use it to process individual hits, retrieving the MVA
  // scores or classifications for further analysis.
  anab::MVAReader<recob::Hit, 4> hitResults(e, "emtrkmichelid:emtrkmichel");
  //  anab::MVAReader<recob::Hit,4> hitResults(e,
  //  "emtrkmichelid:emtrkmichel::pdhdkeepupstage2");

  // Shu: Here is the event-level (1 event ~ several tracks) cycle---
  for (size_t i = 0; i < t0ListHandle->size();
       ++i) { // Shu: t0ListHandle (T0 from pandora)---
    auto &t0 = (*t0ListHandle)[i];
    if (assn3.at(i).isAvailable()) { // Shu: Connect t0 with PFParticle---
      if (assn4.at(assn3.at(i).key())
              .isAvailable()) { // Shu: Connect pandora (reco) track with
                                // PFParticle
        auto &trk = assn4.at(assn3.at(i).key()); // Shu: The track!
        trkid.push_back(trk->ID());
        vtxx.push_back(trk->Vertex().X());
        vtxy.push_back(trk->Vertex().Y());
        vtxz.push_back(trk->Vertex().Z());
        endx.push_back(trk->End().X()); // Shu: track end!
        endy.push_back(trk->End().Y());
        endz.push_back(trk->End().Z());
        pandorat0.push_back(
            t0.Time() * 1e-3); // Shu: T0 here is provided by pandora (TPC reco)
        int endwire = -1;
        int endtpc = -1;
        double endpeakt = -1;
        std::vector<int> wirekeys;

        if (fmthm.isValid()) {
          float disend = 999999;
          auto vhit = fmthm.at(
              trk.key()); // Shu: hit vector (record all hits) of certain track

          auto vmeta = fmthm.data(
              trk.key()); // Shu: access the associated metadata for the hits of
                          // the specified track; the metadata could include
                          // additional information such as the quality of the
                          // association, hit weights, or other relevant
                          // characteristics that describe the relationship
                          // between the track and its associated hits. Here
                          // vmeta is also vector!

          for (size_t ii = 0; ii < vhit.size();
               ++ii) { // loop over all meta data hit
            bool fBadhit = false;
            if (vmeta[ii]->Index() ==
                static_cast<unsigned int>(std::numeric_limits<int>::max())) {
              fBadhit = true; // Shu: bad hit---
              continue;
            }

            if (vmeta[ii]->Index() >=
                trk->NumberTrajectoryPoints()) { // Shu:
                                                 // NumberTrajectoryPoints()
                                                 // returns the total number of
                                                 // trajectory points associated
                                                 // with the track, which could
                                                 // represent the discrete
                                                 // points along the particle’s
                                                 // path as it moves through the
                                                 // detector.

              throw cet::exception("Calorimetry_module.cc")
                  << "Requested track trajectory index " << vmeta[ii]->Index()
                  << " exceeds the total number of trajectory points "
                  << trk->NumberTrajectoryPoints() << " for track index "
                  << trkid.back()
                  << ". Something is wrong with the track reconstruction. "
                     "Please contact tjyang@fnal.gov!!";
            }

            if (!trk->HasValidPoint(
                    vmeta[ii]
                        ->Index())) { // Shu: HasValidPoint(...); It is expected
                                      // to take an index (likely representing a
                                      // specific point along the track's
                                      // trajectory) as an argument and return a
                                      // boolean value indicating whether that
                                      // point is valid.
              fBadhit = true;
              continue;
            }

            auto loc = trk->LocationAtPoint(
                vmeta[ii]->Index()); // Shu: retrieves the location of a
                                     // specific trajectory point on the track
                                     // trk using an index obtained from the
                                     // metadata object vmeta[ii].

            if (fBadhit)
              continue; // HY::If BAD hit, skip this hit and go next

            if (vhit[ii]->WireID().Plane == 2) { // Shu: 2: collection plane---
              wirekeys.push_back(vhit[ii].key());
              float dis = sqrt(pow(loc.X() - endx.back(), 2) +
                               pow(loc.Y() - endy.back(), 2) +
                               pow(loc.Z() - endz.back(), 2));
              if (dis < disend) { // Shu: the shortest (to 'end') dis is finally
                                  // kept at disend---
                disend = dis;     // Def: float disend = 999999;---
                endwire = vhit[ii]->WireID().Wire; // Shu: wire number is
                                                   // kept---
                endpeakt =
                    vhit[ii]
                        ->PeakTime(); // Shu: Peak time of the shortest hit---
                endtpc = vhit[ii]->WireID().TPC;
              }
            }
          }
        }

        float avgmichelscore = 0;
        int nhits = 0;
        for (size_t hitl = 0; hitl < hitList.size();
             hitl++) { // Shu: For hitList, refer to L215

          std::array<float, 4> cnn_out = hitResults.getOutput(
              hitList[hitl]); // Shu: For hitResults, refer to L228
          // Shu: cnn_out is declared as an std::array of four float values. The
          // size 4 suggests that the CNN produces four output values for each
          // hit, which could represent probabilities or scores for different
          // classification categories (e.g., track, shower, Michel electron,
          // background) based on a trained model.

          auto &tracks = thass.at(
              hitList[hitl].key()); // Shu: Refer to P223, for certain hit, we
                                    // can retrieve corresponding track---
          // Shu: tracks here is a reference to a
          // std::vector<art::Ptr<recob::Track>>. vector!!!

          // hit not on the track
          if (std::find(wirekeys.begin(), wirekeys.end(),
                        hitList[hitl].key()) != wirekeys.end())
            continue;
          // Shu: is checking if the key of a particular hit in hitList
          // (specifically hitList[hitl].key()) exists within the wirekeys
          // vector.---

          // hit not on a long track
          if (!tracks.empty() && int(tracks[0].key()) != trkid.back() &&
              trkList[tracks[0].key()]->Length() > 25)
            continue;
          // Shu: If tracks is non-empty, the condition proceeds to the next
          // check. Shu: 'continue' will be executed only if track not empty &
          // longer than 1 & longer than 25 Shu: It means only 'short' tracks be
          // considered in the following code block---

          int planeid = hitList[hitl]->WireID().Plane;
          if (planeid != 2)
            continue; // Shu: Only planeid == 2 is considered in following
                      // code---

          int tpcid = hitList[hitl]->WireID().TPC;
          if (tpcid != endtpc)
            continue; // Shu: refer to L287, only tpcid == tpcid(shortest hit to
                      // end) is considered

          float peakth1 = hitList[hitl]->PeakTime();
          int wireh1 = hitList[hitl]->WireID().Wire;

          // Shu: To judge hit belonging to Michel elctron (close to shortest
          // hit to end (Refer to L285))
          if (std::abs(wireh1 - endwire) < 15 &&
              std::abs(peakth1 - endpeakt) < 100 && tpcid == endtpc) {
            ++nhits;

            avgmichelscore += cnn_out[hitResults.getIndex(
                "michel")]; // Shu: Refer to L297, we know there are four values
                            // of cnn_out, one of each is corresponding to
                            // Michel (Ex: the left three are corresponding to
                            // track, EM shower, etc)
          }
        }
        if (nhits)
          avgmichelscore /= nhits; // Shu: Get the average score of hits---
        michelscore.push_back(avgmichelscore);
        michelhits.push_back(nhits);
      }
    }
  }

  anatree->Fill();
}

void pdhd::MichelTimingMcDeconHD::beginJob() {
  // Implementation of optional member function here.
  anatree = tfs->make<TTree>("anatree", "anatree");
  anatree->Branch("run", &run);
  anatree->Branch("event", &event);
  anatree->Branch("pandorat0", &pandorat0);
  anatree->Branch("trkid", &trkid);
  anatree->Branch("vtxx", &vtxx);
  anatree->Branch("vtxy", &vtxy);
  anatree->Branch("vtxz", &vtxz);
  anatree->Branch("endx", &endx);
  anatree->Branch("endy", &endy);
  anatree->Branch("endz", &endz);
  anatree->Branch("michelscore", &michelscore);
  anatree->Branch("michelhits", &michelhits);
  anatree->Branch("pdchannel", &pdchannel);
  anatree->Branch("pdt0", &pdt0);
  anatree->Branch("nWF", &nWF, "nWF/I");
  //  anatree->Branch("waveform", waveform, "waveform[nWF][1024]/I");
  anatree->Branch("waveform", waveform, "waveform[nWF][1024]/F");
}

DEFINE_ART_MODULE(pdhd::MichelTimingMcDeconHD)
