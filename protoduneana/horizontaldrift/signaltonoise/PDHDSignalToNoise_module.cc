////////////////////////////////////////////////////////////////////////
// Class:       PDHDSignalToNoise
// File:        PDHDSignalToNoise_module.cc
//
// Signal-to-noise measurement for ProtoDUNE-HD, after section 4.6 of the
// ProtoDUNE-SP performance paper.
//
// Signal: the maximum pulse height of the raw waveform, pedestal subtracted,
//         on each wire crossed by a selected cosmic muon track.
// Noise:  the sigma of a Gaussian fit to the raw ADC distribution in the
//         signal-free regions of that channel's waveform.
//
// Reads reconstruction (tracks, hits, deconvolved wires) from the primary
// input file and raw digits from a secondary input file; see the accompanying
// fcl for the source.secondaryFileNames configuration.
//
// Angle handling differs from the paper. "Perpendicular to the wire and also
// perpendicular to the electric field" is a per-plane condition, so instead of
// one global (theta_xz, theta_yz) cut this module stores, per hit:
//   theta_drift - angle of the track out of the anode plane (common to planes)
//   phi_wire    - in-plane angle to that plane's perpendicular-to-wire axis
// Both come from geometry, so no wire angle is hardcoded. Because collection
// and induction wires differ by 35.7 degrees, no single track is perpendicular
// to all three planes at once: the final selection is three separate samples,
// one per plane, sharing the common theta_drift cut. theta_xz / theta_yz are
// stored alongside only so the paper's simpler cut can be reproduced.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"

#include "protoduneana/horizontaldrift/signaltonoise/PDHDSignalToNoiseAlg.h"

#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {
  constexpr double kRad2Deg = 180.0 / M_PI;

  /// Track angles relevant to one wire plane, all in degrees.
  struct PlaneAngles {
    float theta_drift = -999.f; ///< angle out of the anode plane
    float phi_wire = -999.f;    ///< in-plane angle to the perpendicular-to-wire axis
    float pitch = -999.f;       ///< track length per wire crossed
  };

  /// Takes bare components rather than a vector type: recob::tracking::Vector_t
  /// and geo::Vector_t are distinct typedefs, so mixing them in a dot product
  /// is not portable across LArSoft versions.
  PlaneAngles ComputePlaneAngles(geo::PlaneGeo const& plane, double dx, double dy, double dz)
  {
    PlaneAngles a;

    const double mag = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!(mag > 0.)) return a;
    const geo::Vector_t d(dx / mag, dy / mag, dz / mag);

    auto const& n = plane.GetNormalDirection();          // drift direction
    auto const& p = plane.GetIncreasingWireDirection();  // in-plane, perpendicular to wires

    const double dn = d.Dot(n);
    a.theta_drift = std::asin(std::min(1.0, std::abs(dn))) * kRad2Deg;

    // In-plane projection, then its angle to the perpendicular-to-wire axis.
    const geo::Vector_t inplane = d - dn * n;
    const double ip = std::sqrt(inplane.Mag2());
    if (ip > 1e-6) {
      const double c = std::min(1.0, std::abs(inplane.Dot(p)) / ip);
      a.phi_wire = std::acos(c) * kRad2Deg;
    }

    // Track length per wire crossed. This folds both angles into exactly the
    // quantity that governs how much charge lands on one wire.
    const double cosgamma = std::abs(d.Dot(p));
    if (cosgamma > 1e-6) a.pitch = plane.WirePitch() / cosgamma;

    return a;
  }
}

namespace pdhd {
  class PDHDSignalToNoise;
}

class pdhd::PDHDSignalToNoise : public art::EDAnalyzer {
public:
  explicit PDHDSignalToNoise(fhicl::ParameterSet const& p);

  PDHDSignalToNoise(PDHDSignalToNoise const&) = delete;
  PDHDSignalToNoise(PDHDSignalToNoise&&) = delete;
  PDHDSignalToNoise& operator=(PDHDSignalToNoise const&) = delete;
  PDHDSignalToNoise& operator=(PDHDSignalToNoise&&) = delete;

  void analyze(art::Event const& e) override;
  void beginJob() override;
  void endJob() override;

private:
  void ResetHit();

  PDHDSignalToNoiseAlg fAlg;

  // Configuration
  std::string fTrackLabel;
  std::string fCaloLabel;
  art::InputTag fWireTag;
  std::vector<std::string> fRawLabels;
  std::vector<std::string> fRawNicks;
  float fMinTrackLength;
  float fMaxThetaDrift; ///< loose pre-cut, degrees
  float fMaxPhiInPlane; ///< loose pre-cut on the in-plane angle, degrees
  bool fUseChannelStatus;
  bool fCaloIndexIsHitKey; ///< true for GnocchiCalorimetry, false for Calorimetry
  unsigned fNLabels;

  // Trees
  TTree* fHitTree = nullptr;
  TTree* fChanTree = nullptr;

  // Event / track / hit scalars
  int fRun, fSubRun, fEvent;
  int fTrkID, fNTracks;
  float fTrkLen, fTrkThetaXZ, fTrkThetaYZ;
  int fChannel, fTPC, fAPA, fPlane, fWire;
  float fPeakT;
  float fHitThetaXZ, fHitThetaYZ, fThetaDrift, fPhiWire;
  float fPitchGeo, fPitchCalo;
  float fX, fY, fZ;
  int fRoiLo, fRoiHi;
  int fChanGood;

  // Per-raw-label quantities; one entry per element, branched individually so
  // the output tree stays flat. Never resized after beginJob.
  std::vector<float> fAmpHit, fRawHit, fAmpRoi, fPed, fSigma, fRms, fChi2Ndf;
  std::vector<float> fDTick; ///< (raw max tick, in wire frame) - PeakTime
  std::vector<int> fTickHit, fNFree, fStatus, fNSamp;

  // Channel tree
  int fCChannel, fCTPC, fCAPA, fCPlane, fCWire, fCGood, fCNFreeMask;
  std::vector<float> fCPed, fCSigma, fCRms, fCChi2Ndf;
  std::vector<int> fCNFree, fCStatus;

  // Counters
  long fNEvents = 0, fNTracksSeen = 0, fNTracksSel = 0, fNHits = 0;
  long fNNoiseOK = 0, fNNoiseFail = 0;
  long fNEventsWithRaw = 0, fNEventsNoRaw = 0;
  std::vector<long> fFreeTickTotal, fFreeTickChans; // per plane, for the mask sanity check

  geo::WireReadoutGeom const* fWireReadout = nullptr;
  lariov::ChannelStatusProvider const* fChanStatus = nullptr;
};

pdhd::PDHDSignalToNoise::PDHDSignalToNoise(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fAlg(p.get<fhicl::ParameterSet>("Alg", fhicl::ParameterSet()))
  , fTrackLabel(p.get<std::string>("TrackModuleLabel", "pandoraTrack"))
  , fCaloLabel(p.get<std::string>("CalorimetryModuleLabel", "pandoraGnocchiCalo"))
  , fWireTag(p.get<std::string>("WireModuleLabel", "wclsdatahd:gauss"))
  , fRawLabels(p.get<std::vector<std::string>>(
      "RawDigitLabels", {"tpcrawdecoder:daq", "wclsdatahdfilter:raw"}))
  , fRawNicks(p.get<std::vector<std::string>>("RawDigitNicknames", {"unfilt", "filt"}))
  , fMinTrackLength(p.get<float>("MinTrackLength", 100.))
  , fMaxThetaDrift(p.get<float>("MaxThetaDrift", 35.))
  , fMaxPhiInPlane(p.get<float>("MaxPhiInPlane", 60.))
  , fUseChannelStatus(p.get<bool>("UseChannelStatus", true))
  , fCaloIndexIsHitKey(p.get<bool>("CaloIndexIsHitKey", true))
{
  if (fRawLabels.size() != fRawNicks.size()) {
    throw cet::exception("PDHDSignalToNoise")
      << "RawDigitLabels (" << fRawLabels.size() << ") and RawDigitNicknames ("
      << fRawNicks.size() << ") must have the same length.\n";
  }
  fNLabels = fRawLabels.size();

  fAmpHit.assign(fNLabels, -999.f);
  fRawHit.assign(fNLabels, -999.f);
  fAmpRoi.assign(fNLabels, -999.f);
  fPed.assign(fNLabels, -999.f);
  fSigma.assign(fNLabels, -999.f);
  fRms.assign(fNLabels, -999.f);
  fChi2Ndf.assign(fNLabels, -999.f);
  fDTick.assign(fNLabels, -999.f);
  fTickHit.assign(fNLabels, -1);
  fNFree.assign(fNLabels, 0);
  fStatus.assign(fNLabels, -1);
  fNSamp.assign(fNLabels, 0);

  fCPed.assign(fNLabels, -999.f);
  fCSigma.assign(fNLabels, -999.f);
  fCRms.assign(fNLabels, -999.f);
  fCChi2Ndf.assign(fNLabels, -999.f);
  fCNFree.assign(fNLabels, 0);
  fCStatus.assign(fNLabels, -1);

  fFreeTickTotal.assign(3, 0);
  fFreeTickChans.assign(3, 0);
}

void pdhd::PDHDSignalToNoise::beginJob()
{
  fWireReadout = &art::ServiceHandle<geo::WireReadout>()->Get();
  if (fUseChannelStatus) {
    fChanStatus = &art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();
  }

  art::ServiceHandle<art::TFileService> tfs;

  fHitTree = tfs->make<TTree>("hits", "one entry per hit on a selected track");
  fHitTree->Branch("run", &fRun, "run/I");
  fHitTree->Branch("subrun", &fSubRun, "subrun/I");
  fHitTree->Branch("event", &fEvent, "event/I");
  fHitTree->Branch("ntracks", &fNTracks, "ntracks/I");
  fHitTree->Branch("trkid", &fTrkID, "trkid/I");
  fHitTree->Branch("trklen", &fTrkLen, "trklen/F");
  fHitTree->Branch("trk_thetaxz", &fTrkThetaXZ, "trk_thetaxz/F");
  fHitTree->Branch("trk_thetayz", &fTrkThetaYZ, "trk_thetayz/F");
  fHitTree->Branch("channel", &fChannel, "channel/I");
  fHitTree->Branch("tpc", &fTPC, "tpc/I");
  fHitTree->Branch("apa", &fAPA, "apa/I");
  fHitTree->Branch("plane", &fPlane, "plane/I");
  fHitTree->Branch("wire", &fWire, "wire/I");
  fHitTree->Branch("peakt", &fPeakT, "peakt/F");
  fHitTree->Branch("hit_thetaxz", &fHitThetaXZ, "hit_thetaxz/F");
  fHitTree->Branch("hit_thetayz", &fHitThetaYZ, "hit_thetayz/F");
  fHitTree->Branch("theta_drift", &fThetaDrift, "theta_drift/F");
  fHitTree->Branch("phi_wire", &fPhiWire, "phi_wire/F");
  fHitTree->Branch("pitch_geo", &fPitchGeo, "pitch_geo/F");
  fHitTree->Branch("pitch_calo", &fPitchCalo, "pitch_calo/F");
  fHitTree->Branch("x", &fX, "x/F");
  fHitTree->Branch("y", &fY, "y/F");
  fHitTree->Branch("z", &fZ, "z/F");
  fHitTree->Branch("roi_lo", &fRoiLo, "roi_lo/I");
  fHitTree->Branch("roi_hi", &fRoiHi, "roi_hi/I");
  fHitTree->Branch("chan_good", &fChanGood, "chan_good/I");

  for (unsigned i = 0; i < fNLabels; ++i) {
    const std::string& n = fRawNicks[i];
    fHitTree->Branch(("amp_hit_" + n).c_str(), &fAmpHit[i], ("amp_hit_" + n + "/F").c_str());
    fHitTree->Branch(("raw_hit_" + n).c_str(), &fRawHit[i], ("raw_hit_" + n + "/F").c_str());
    fHitTree->Branch(("tick_hit_" + n).c_str(), &fTickHit[i], ("tick_hit_" + n + "/I").c_str());
    fHitTree->Branch(("dtick_hit_" + n).c_str(), &fDTick[i], ("dtick_hit_" + n + "/F").c_str());
    fHitTree->Branch(("nsamp_" + n).c_str(), &fNSamp[i], ("nsamp_" + n + "/I").c_str());
    fHitTree->Branch(("amp_roi_" + n).c_str(), &fAmpRoi[i], ("amp_roi_" + n + "/F").c_str());
    fHitTree->Branch(("ped_" + n).c_str(), &fPed[i], ("ped_" + n + "/F").c_str());
    fHitTree->Branch(("sigma_" + n).c_str(), &fSigma[i], ("sigma_" + n + "/F").c_str());
    fHitTree->Branch(("rms_" + n).c_str(), &fRms[i], ("rms_" + n + "/F").c_str());
    fHitTree->Branch(("nfree_" + n).c_str(), &fNFree[i], ("nfree_" + n + "/I").c_str());
    fHitTree->Branch(("status_" + n).c_str(), &fStatus[i], ("status_" + n + "/I").c_str());
  }

  fChanTree = tfs->make<TTree>("channels", "one entry per channel touched by a selected track");
  fChanTree->Branch("run", &fRun, "run/I");
  fChanTree->Branch("subrun", &fSubRun, "subrun/I");
  fChanTree->Branch("event", &fEvent, "event/I");
  fChanTree->Branch("channel", &fCChannel, "channel/I");
  fChanTree->Branch("tpc", &fCTPC, "tpc/I");
  fChanTree->Branch("apa", &fCAPA, "apa/I");
  fChanTree->Branch("plane", &fCPlane, "plane/I");
  fChanTree->Branch("wire", &fCWire, "wire/I");
  fChanTree->Branch("chan_good", &fCGood, "chan_good/I");
  fChanTree->Branch("nfree_mask", &fCNFreeMask, "nfree_mask/I");
  for (unsigned i = 0; i < fNLabels; ++i) {
    const std::string& n = fRawNicks[i];
    fChanTree->Branch(("ped_" + n).c_str(), &fCPed[i], ("ped_" + n + "/F").c_str());
    fChanTree->Branch(("sigma_" + n).c_str(), &fCSigma[i], ("sigma_" + n + "/F").c_str());
    fChanTree->Branch(("rms_" + n).c_str(), &fCRms[i], ("rms_" + n + "/F").c_str());
    fChanTree->Branch(("chi2ndf_" + n).c_str(), &fCChi2Ndf[i], ("chi2ndf_" + n + "/F").c_str());
    fChanTree->Branch(("nfree_" + n).c_str(), &fCNFree[i], ("nfree_" + n + "/I").c_str());
    fChanTree->Branch(("status_" + n).c_str(), &fCStatus[i], ("status_" + n + "/I").c_str());
  }
}

void pdhd::PDHDSignalToNoise::ResetHit()
{
  for (unsigned i = 0; i < fNLabels; ++i) {
    fAmpHit[i] = fRawHit[i] = fAmpRoi[i] = -999.f;
    fPed[i] = fSigma[i] = fRms[i] = fChi2Ndf[i] = fDTick[i] = -999.f;
    fTickHit[i] = -1;
    fNFree[i] = 0;
    fStatus[i] = -1;
    fNSamp[i] = 0;
  }
  fPitchCalo = -999.f;
  fRoiLo = fRoiHi = -1;
}

void pdhd::PDHDSignalToNoise::analyze(art::Event const& e)
{
  fRun = e.run();
  fSubRun = e.subRun();
  fEvent = e.id().event();
  ++fNEvents;

  // ---- Products -----------------------------------------------------------
  art::Handle<std::vector<recob::Track>> trackHandle;
  if (!e.getByLabel(fTrackLabel, trackHandle)) {
    mf::LogWarning("PDHDSignalToNoise") << "no tracks under " << fTrackLabel;
    return;
  }
  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, trackHandle);
  fNTracks = tracks.size();

  art::Handle<std::vector<recob::Wire>> wireHandle;
  if (!e.getByLabel(fWireTag, wireHandle)) {
    mf::LogWarning("PDHDSignalToNoise") << "no wires under " << fWireTag.encode();
    return;
  }

  // Raw digits come from the secondary input file. getByLabel, not
  // getValidHandle: a missing secondary event is silent in art.
  std::vector<art::Handle<std::vector<raw::RawDigit>>> rawHandles(fNLabels);
  for (unsigned i = 0; i < fNLabels; ++i) {
    if (!e.getByLabel(art::InputTag(fRawLabels[i]), rawHandles[i])) {
      mf::LogWarning("PDHDSignalToNoise")
        << "no raw digits under " << fRawLabels[i] << " for run " << fRun << " event " << fEvent
        << " -- is the secondary input file configured and does it contain this event?";
      ++fNEventsNoRaw;
      return;
    }
  }
  ++fNEventsWithRaw;

  art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackHandle, e, fTrackLabel);
  art::FindManyP<anab::Calorimetry> fmcal(trackHandle, e, fCaloLabel);

  // ---- Per-event lookup tables -------------------------------------------
  std::map<raw::ChannelID_t, std::size_t> wireIndex;
  for (std::size_t i = 0; i < wireHandle->size(); ++i)
    wireIndex[(*wireHandle)[i].Channel()] = i;

  std::vector<std::map<raw::ChannelID_t, std::size_t>> rawIndex(fNLabels);
  for (unsigned l = 0; l < fNLabels; ++l)
    for (std::size_t i = 0; i < rawHandles[l]->size(); ++i)
      rawIndex[l][(*rawHandles[l])[i].Channel()] = i;

  // Computed lazily, only for channels a selected track actually touches.
  std::map<raw::ChannelID_t, WireMask> maskCache;
  std::vector<std::map<raw::ChannelID_t, std::vector<short>>> adcCache(fNLabels);
  std::vector<std::map<raw::ChannelID_t, NoiseResult>> noiseCache(fNLabels);
  std::map<raw::ChannelID_t, bool> chanWritten;

  // ---- Tracks -------------------------------------------------------------
  for (std::size_t itrk = 0; itrk < tracks.size(); ++itrk) {
    auto const& track = *tracks[itrk];
    ++fNTracksSeen;

    if (track.Length() < fMinTrackLength) continue;

    auto const dir0 = track.StartDirection();
    fTrkThetaXZ = std::atan2(dir0.X(), dir0.Z()) * kRad2Deg;
    fTrkThetaYZ = std::atan2(dir0.Y(), dir0.Z()) * kRad2Deg;
    fTrkLen = track.Length();
    fTrkID = track.ID();

    if (!fmthm.isValid()) continue;
    auto vhit = fmthm.at(itrk);
    auto vmeta = fmthm.data(itrk);
    if (vhit.empty()) continue;

    // Loose pre-cut, applied on the first valid trajectory point. It has to be
    // the union of the three per-plane acceptances, not a theta_yz cut, or the
    // induction samples would be thrown away here.
    bool anyPass = false;
    for (std::size_t ii = 0; ii < vhit.size() && !anyPass; ++ii) {
      if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max()))
        continue;
      if (vmeta[ii]->Index() >= track.NumberTrajectoryPoints()) continue;
      if (!track.HasValidPoint(vmeta[ii]->Index())) continue;
      auto const* plane = fWireReadout->PlanePtr(vhit[ii]->WireID().planeID());
      if (!plane) continue;
      auto const pd = track.DirectionAtPoint(vmeta[ii]->Index());
      const auto a = ComputePlaneAngles(*plane, pd.X(), pd.Y(), pd.Z());
      if (a.theta_drift >= 0.f && a.theta_drift < fMaxThetaDrift && a.phi_wire >= 0.f &&
          a.phi_wire < fMaxPhiInPlane)
        anyPass = true;
    }
    if (!anyPass) continue;
    ++fNTracksSel;

    // Track pitch from calorimetry, keyed by trajectory point, as a cross-check
    // on the geometric pitch. Absent if the file was made without calorimetry.
    // Keyed by (plane, index). Two traps here, both silent:
    //
    //  1. There is one Calorimetry object per plane and their index ranges
    //     overlap, so a map keyed on the index alone lets planes overwrite
    //     one another.
    //  2. What TpIndices() actually holds depends on the producer. The old
    //     Calorimetry module stores trajectory-point indices, as the name
    //     suggests, but GnocchiCalorimetry_module.cc:362-363 has that line
    //     commented out and pushes hits[i].key() instead -- i.e. hit keys.
    //     Matching the wrong index space still finds entries, because the two
    //     ranges overlap, and yields a plausible-looking wrong pitch.
    std::map<std::pair<unsigned, std::size_t>, float> caloPitch;
    if (fmcal.isValid()) {
      for (auto const& calo : fmcal.at(itrk)) {
        if (!calo->PlaneID().isValid) continue;
        const unsigned pl = calo->PlaneID().Plane;
        auto const& idx = calo->TpIndices();
        auto const& pitch = calo->TrkPitchVec();
        for (std::size_t k = 0; k < idx.size() && k < pitch.size(); ++k)
          caloPitch[{pl, idx[k]}] = pitch[k];
      }
    }

    // ---- Hits on this track ----------------------------------------------
    for (std::size_t ii = 0; ii < vhit.size(); ++ii) {
      // The three guards are all necessary; see michelremoving_module.cc.
      if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max()))
        continue;
      if (vmeta[ii]->Index() >= track.NumberTrajectoryPoints()) continue;
      if (!track.HasValidPoint(vmeta[ii]->Index())) continue;

      auto const& wid = vhit[ii]->WireID();
      auto const* plane = fWireReadout->PlanePtr(wid.planeID());
      if (!plane) continue;

      ResetHit();

      const std::size_t ipt = vmeta[ii]->Index();
      auto const loc = track.LocationAtPoint(ipt);
      auto const dir = track.DirectionAtPoint(ipt);
      const auto ang = ComputePlaneAngles(*plane, dir.X(), dir.Y(), dir.Z());

      fChannel = vhit[ii]->Channel();
      fTPC = wid.TPC;
      fAPA = wid.TPC / 2; // PDHD: two drift volumes per APA
      fPlane = wid.Plane;
      fWire = wid.Wire;
      fPeakT = vhit[ii]->PeakTime();
      fHitThetaXZ = std::atan2(dir.X(), dir.Z()) * kRad2Deg;
      fHitThetaYZ = std::atan2(dir.Y(), dir.Z()) * kRad2Deg;
      fThetaDrift = ang.theta_drift;
      fPhiWire = ang.phi_wire;
      fPitchGeo = ang.pitch;
      fX = loc.X();
      fY = loc.Y();
      fZ = loc.Z();

      // See the caloPitch comment: which index to look up depends on the
      // calorimetry producer. pandoraGnocchiCalo indexes by hit key.
      const std::size_t calokey = fCaloIndexIsHitKey ? vhit[ii].key() : ipt;
      auto cp = caloPitch.find({static_cast<unsigned>(wid.Plane), calokey});
      if (cp != caloPitch.end()) fPitchCalo = cp->second;

      fChanGood = (fChanStatus == nullptr) ? 1 : (fChanStatus->IsGood(fChannel) ? 1 : 0);

      auto wi = wireIndex.find(fChannel);
      if (wi == wireIndex.end()) continue;
      auto const& wire = (*wireHandle)[wi->second];
      const std::size_t nticksWire = wire.NSignal();

      // Signal-free mask, cached per channel.
      auto mit = maskCache.find(fChannel);
      if (mit == maskCache.end()) {
        mit = maskCache.emplace(fChannel, fAlg.MakeWireMask(wire, nticksWire)).first;
        if (fPlane < 3) {
          fFreeTickTotal[fPlane] += mit->second.nfree;
          ++fFreeTickChans[fPlane];
        }
      }
      auto const& mask = mit->second;

      const auto roi = fAlg.RoiContaining(mask, static_cast<int>(std::round(fPeakT)));
      fRoiLo = roi.first;
      fRoiHi = roi.second;

      for (unsigned l = 0; l < fNLabels; ++l) {
        auto ri = rawIndex[l].find(fChannel);
        if (ri == rawIndex[l].end()) continue;
        auto const& digit = (*rawHandles[l])[ri->second];

        auto ait = adcCache[l].find(fChannel);
        if (ait == adcCache[l].end()) {
          std::vector<short> adcs(digit.Samples());
          raw::Uncompress(digit.ADCs(), adcs, static_cast<int>(digit.GetPedestal()),
                          digit.Compression());
          ait = adcCache[l].emplace(fChannel, std::move(adcs)).first;
        }
        auto const& adcs = ait->second;
        if (adcs.empty() || nticksWire == 0) continue;

        // The raw products are not all on the same time base. PDHD delivers
        // tpcrawdecoder:daq at 5859 ticks of 512 ns, while WireCell resamples
        // to 5999 ticks of 500 ns over the same 3 ms; hit times and the ROI
        // mask live in the latter frame. Scale rather than assume.
        const double scale = double(adcs.size()) / double(nticksWire);
        fNSamp[l] = static_cast<int>(adcs.size());

        auto nit = noiseCache[l].find(fChannel);
        if (nit == noiseCache[l].end()) {
          nit = noiseCache[l].emplace(fChannel, fAlg.FitNoise(adcs, mask, scale)).first;
          if (nit->second.ok)
            ++fNNoiseOK;
          else
            ++fNNoiseFail;
        }
        auto const& noise = nit->second;

        // Fall back to the digit's own pedestal if the fit failed, so the
        // amplitude is still meaningful; status_* records which happened.
        const float ped = noise.ok ? noise.ped : digit.GetPedestal();

        fPed[l] = noise.ped;
        fSigma[l] = noise.sigma;
        fRms[l] = noise.rms;
        fChi2Ndf[l] = noise.chi2ndf;
        fNFree[l] = noise.nfree;
        fStatus[l] = noise.ok ? 0 : 1;

        const int c = PDHDSignalToNoiseAlg::ToDigitTick(fPeakT, scale);
        const int hw =
          std::max(1, static_cast<int>(std::lround(fAlg.HitHalfWindow() * scale)));
        const auto mh = fAlg.MaxInWindow(adcs, c - hw, c + hw + 1, ped);
        fAmpHit[l] = mh.amp;
        fRawHit[l] = mh.raw;
        fTickHit[l] = mh.tick;
        // Reported back in the wire frame so the alignment check is uniform
        // across products on different time bases.
        if (mh.tick >= 0) fDTick[l] = static_cast<float>(mh.tick / scale - fPeakT);

        if (roi.first >= 0) {
          const auto mr = fAlg.MaxInWindow(adcs,
                                           PDHDSignalToNoiseAlg::ToDigitTick(roi.first, scale),
                                           PDHDSignalToNoiseAlg::ToDigitTick(roi.second, scale),
                                           ped);
          fAmpRoi[l] = mr.amp;
        }
      }

      fHitTree->Fill();
      ++fNHits;

      // One channels-tree entry per channel per event.
      if (!chanWritten[fChannel]) {
        chanWritten[fChannel] = true;
        fCChannel = fChannel;
        fCTPC = fTPC;
        fCAPA = fAPA;
        fCPlane = fPlane;
        fCWire = fWire;
        fCGood = fChanGood;
        fCNFreeMask = mask.nfree;
        for (unsigned l = 0; l < fNLabels; ++l) {
          auto nit = noiseCache[l].find(fChannel);
          if (nit == noiseCache[l].end()) continue;
          fCPed[l] = nit->second.ped;
          fCSigma[l] = nit->second.sigma;
          fCRms[l] = nit->second.rms;
          fCChi2Ndf[l] = nit->second.chi2ndf;
          fCNFree[l] = nit->second.nfree;
          fCStatus[l] = nit->second.ok ? 0 : 1;
        }
        fChanTree->Fill();
      }
    }
  }
}

void pdhd::PDHDSignalToNoise::endJob()
{
  // The commonest way to misconfigure this job is to pass the primary file with
  // -s while secondaryFileNames sits in the fcl: the "a" key then no longer
  // byte-matches source.fileNames, art silently supplies no secondaries, and
  // the job runs to completion writing an empty tree. Refuse to exit quietly.
  if (fNEvents > 0 && fNEventsWithRaw == 0) {
    throw cet::exception("PDHDSignalToNoise")
      << "processed " << fNEvents << " events and never found raw digits under any of the "
      << fNLabels << " configured labels.\n"
      << "The secondary input file is almost certainly not being consulted. Check that\n"
      << "source.secondaryFileNames is set and that its 'a' key is byte-identical to the\n"
      << "source.fileNames entry -- passing the primary with -s breaks that match.\n";
  }

  mf::LogInfo log("PDHDSignalToNoise");
  log << "\n=== PDHDSignalToNoise summary ===\n"
      << "  events              " << fNEvents << "\n"
      << "  events with raw     " << fNEventsWithRaw << " (missing in " << fNEventsNoRaw << ")\n"
      << "  tracks seen         " << fNTracksSeen << "\n"
      << "  tracks selected     " << fNTracksSel << "\n"
      << "  hits written        " << fNHits << "\n"
      << "  noise fits ok/fail  " << fNNoiseOK << " / " << fNNoiseFail << "\n";

  // Mask sanity check. If the free-tick fraction is ~0 the zero-run assumption
  // about the dense WireCell output is wrong and every S/N number downstream is
  // meaningless, so say so loudly rather than let it pass.
  bool bad = false;
  for (int pl = 0; pl < 3; ++pl) {
    const double mean =
      fFreeTickChans[pl] ? double(fFreeTickTotal[pl]) / fFreeTickChans[pl] : 0.;
    log << "  plane " << pl << ": " << fFreeTickChans[pl] << " channels, mean free ticks "
        << mean << "\n";
    if (fFreeTickChans[pl] > 0 && mean < double(fAlg.MinFreeTicks())) bad = true;
  }
  if (bad) {
    log << "  WARNING: mean free-tick count is below MinFreeTicks on at least one plane.\n"
        << "  The signal-free mask is not finding quiet regions -- check that the wire\n"
        << "  product really is the dense WireCell output before trusting any noise value.\n";
  }
}

DEFINE_ART_MODULE(pdhd::PDHDSignalToNoise)
