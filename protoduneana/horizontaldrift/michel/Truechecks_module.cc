////////////////////////////////////////////////////////////////////////
// Class:       Truechecks
// Plugin Type: analyzer
// File:        Truechecks_module.cc
//
// Generated at Tue Feb 4 15:23:48 2025 by Jeremy Quelin Lechevranton
//
// Updated by Shuaixiang (Shu) Zhang at Nov 24, 2025
// The updated module applies stricter cuts to select Michel signals
////////////////////////////////////////////////////////////////////////

// Art includes
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

// LArSoft includes
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "larcorealg/Geometry/Exceptions.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
// #include "lardata/ArtDataHelper/TrackUtils.h"

// ProtoDUNE includes
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"

// ROOT includes
#include "TCanvas.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TTree.h"

// #include "ROOT/RVec.hxx"
// #include "RVec.hxx"

// std includes
#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ana {
class Truechecks;
struct Binning {
  int n;
  float min, max;
};
} // namespace ana

class ana::Truechecks : public art::EDAnalyzer {
public:
  explicit Truechecks(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Truechecks(Truechecks const &) = delete;
  Truechecks(Truechecks &&) = delete;
  Truechecks &operator=(Truechecks const &) = delete;
  Truechecks &operator=(Truechecks &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Utilities
  art::ServiceHandle<art::TFileService> tfs;

  const geo::GeometryCore *asGeo;
  const geo::WireReadoutGeom *asWire;
  const detinfo::DetectorPropertiesService *asDetProp;
  const detinfo::DetectorClocksService *asDetClocks;

  protoana::ProtoDUNETruthUtils truthUtil;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  // Conversion factors
  //float feltoMeV =
      23.6 * 1e-6 / 0.7; // 23.6 eV/e- * 1e-6 MeV/eV / 0.7 recombination factor
  float fADCtoMeV =
      200 * 23.6 * 1e-6 / 0.7; // 200 e-/ADC.tick * 23.6 eV/e- * 1e-6 MeV/eV /
                               // 0.7 recombination factor
  float fChannelPitch;
  float fDriftVelocity; // cm/µs
  float fSamplingRate;  // µs/tick

  // Verbosity
  int iLogLevel;
  enum EnumFlag { kImportant, kBasics, kInfos, kDetails };

  // Diagnostic Variables
  std::map<std::string, unsigned> map_mup_endproc;
  std::map<std::string, unsigned> map_mum_endproc;
  unsigned n_mup = 0, n_mum = 0;
  // unsigned n_mep = 0, n_mem = 0;
  unsigned n_cme_wh = 0, n_cme_nh = 0;
  float mean_cme_h = 0;

  // Products
  std::vector<std::vector<std::string>> vvsProducts;
  art::InputTag tag_mcp, tag_sed, tag_wir, tag_hit, tag_clu, tag_trk, tag_spt;

  bool Log(bool cond, int flag, int tab, std::string msg, std::string succ,
           std::string fail);
  std::string GetParticleName(int pdg);
};

ana::Truechecks::Truechecks(fhicl::ParameterSet const &p)
    : EDAnalyzer{p}, iLogLevel(p.get<int>("LogLevel")),
      vvsProducts(p.get<std::vector<std::vector<std::string>>>("Products")) {
  if (iLogLevel >= kBasics)
    std::cout << "\n\n\033[93m"
              << "Truechecks::Truechecks: ================="
              << "\033[0m" << std::endl;
  // Basic Utilities
  asGeo = &*art::ServiceHandle<geo::Geometry>();
  asWire = &art::ServiceHandle<geo::WireReadout>()->Get();
  asDetProp = &*art::ServiceHandle<detinfo::DetectorPropertiesService>();
  asDetClocks = &*art::ServiceHandle<detinfo::DetectorClocksService>();

  geo::WireGeo const wiregeo1 =
      asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 0});
  geo::WireGeo const wiregeo2 =
      asWire->Wire(geo::WireID{geo::PlaneID{geo::TPCID{0, 0}, geo::kW}, 1});
  fChannelPitch = geo::WireGeo::WirePitch(wiregeo1, wiregeo2);

  auto const clockData = asDetClocks->DataForJob();
  auto const detProp = asDetProp->DataForJob(clockData);

  // Retrieving product tags
  for (std::vector<std::string> prod : vvsProducts) {

    const std::string process = prod[0], label = prod[1], instance = prod[2],
                      type = prod[3];

    const art::InputTag tag = art::InputTag(label, instance);

    if (type == "simb::MCParticle")
      tag_mcp = tag;
    else if (type == "sim::SimEnergyDeposit")
      tag_sed = tag;
    else if (type == "recob::Hit")
      tag_hit = tag;
    else if (type == "recob::Wire")
      tag_wir = tag;
    else if (type == "recob::Cluster")
      tag_clu = tag;
    else if (type == "recob::Track")
      tag_trk = tag;
    else if (type == "recob::SpacePoint")
      tag_spt = tag;
  }

  // fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
  // fDriftVelocity = detProp.DriftVelocity();

  if (iLogLevel >= kBasics)
    std::cout << "\033[93m"
              << "End of Truechecks::Truechecks =========="
              << "\033[0m\n"
              << std::endl;
}

void ana::Truechecks::analyze(art::Event const &e) {

  if (iLogLevel >= kBasics) {
    std::cout << "\n\n\n\033[93m"
              << "Truechecks::analyze: Initialization evt#" << std::setw(5)
              << e.id().event() << " ===================================="
              << "\033[0m" << std::endl;
  }

  auto const clockData = asDetClocks->DataFor(e);
  auto const detProp = asDetProp->DataFor(e, clockData);
  fSamplingRate = detinfo::sampling_rate(clockData) * 1e-3;
  fDriftVelocity = detProp.DriftVelocity();

  auto const &vh_trk = e.getValidHandle<std::vector<recob::Track>>(tag_trk);
  std::vector<art::Ptr<recob::Track>> vp_trk;
  art::fill_ptr_vector(vp_trk, vh_trk);

  std::vector<unsigned> n_dau;
  std::vector<unsigned> n_muioni;
  std::vector<std::vector<int>> dau_pdg;
  std::vector<std::vector<std::string>> dau_process;

  // Print MC truth========================================================
  std::cout << "True-level Michel scan for current event!" << std::endl;
  std::cout << "Within TPC (-356<x<356, 0<y<610, 0<z<465)" << std::endl;

  // Get full list of MC particles in event
  auto const &mcps = *e.getValidHandle<std::vector<simb::MCParticle>>(tag_mcp);

  // Loop over all particles to find candidate Michel electrons
  for (auto const &part : mcps) {
    // Step 1: Must be e-/e+
    if (std::abs(part.PdgCode()) != 11)
      continue;

    // Step 2: Parent must be a muon
    int mother_id = part.Mother();
    if (mother_id <= 0)
      continue; // <-- prevents PI (print) warnings
    const simb::MCParticle *parent = pi_serv->TrackIdToParticle_P(mother_id);
    if (!parent)
      continue;
    if (std::abs(parent->PdgCode()) != 13)
      continue;

    // Step 3: Require Michel e⁻ to be within TPC
    double x = part.Vx();
    double y = part.Vy();
    double z = part.Vz();
    double t = part.T();
    if (x < -356 || x > 356)
      continue;
    if (y < 0 || y > 610)
      continue;
    if (z < 0 || z > 465)
      continue;

    // Step 4: Michel e⁻ must originate near muon's end point
    double dx = x - parent->EndX();
    double dy = y - parent->EndY();
    double dz = z - parent->EndZ();
    double dist2 = dx * dx + dy * dy + dz * dz;
    if (dist2 > 25.0)
      continue; // >5 cm^2 → not from muon decay point

    // Step 5: Must have at least ONE neutrino (ν_e or ν_μ) from the SAME vertex
    // (space + time) Use squared distance in cm^2 to avoid sqrt.
    const double r2Match = 25.0; // (5 cm)^2
    const double tMatch = 5.0;   // ns

    bool hasNuSibling = false; // at least one of |PDG| == 12 or 14
    int nNuAtVertex = 0;
    for (auto const &sibling : mcps) {
      if (sibling.Mother() != parent->TrackId())
        continue;
      int apdg = std::abs(sibling.PdgCode());
      if (apdg != 12 && apdg != 14)
        continue; // only ν_e / ν_μ (and anti)

      double ddx = sibling.Vx() - x;
      double ddy = sibling.Vy() - y;
      double ddz = sibling.Vz() - z;
      double ddt = std::abs(sibling.T() - t);
      if (ddx * ddx + ddy * ddy + ddz * ddz <= r2Match && ddt <= tMatch) {
        hasNuSibling = true;
        ++nNuAtVertex; // keep count for diagnostics
      }
    }
    if (!hasNuSibling)
      continue;

    // Step 6: Simply reject delta-like electrons outright
    const std::string eproc = part.Process();
    auto isDeltaLike = [&](const std::string &p) {
      return (p == "muIoni" || p == "hIoni" || p == "eIoni" || p == "muBrems" ||
              p == "compt" || p == "phot" || p == "conv" || p == "annihil");
    };
    if (isDeltaLike(eproc))
      continue;

    // Step 7: Request Kinetic Energy > 1 MeV
    double KE = (part.E() - part.Mass()) * 1e3; // kinetic energy in MeV
    if (KE < 1.0)
      continue;

    // If we reach here → robust true Michel electron!
    //        double lifetime_end_us = (parent->EndT() - parent->T()) * 1e-3; //
    //        official: muon birth → muon end
    double lifetime_e_us =
        (part.T() - parent->T()) * 1e-3; // muon birth → Michel birth

    std::cout << "Michel (PDG = " << part.PdgCode()
              << ") (trackID = " << part.TrackId() << "): Kinetic = " << KE
              << " MeV, lifetime = " << lifetime_e_us << " us"
              << "\n  Start (x,y,z,t) = (" << x << ", " << y << ", " << z
              << ", " << part.T() << " ns)"
              << "\n    e⁻ Process = " << part.Process()
              << ", muon end process = " << parent->EndProcess()
              << "\n  Mu(x0,y0,z0,t0) = (" << parent->Vx() << ", "
              << parent->Vy() << ", " << parent->Vz() << ", " << parent->T()
              << " ns)"
              << "\n  Mu(x1,y1,z1,t1) = (" << parent->EndX() << ", "
              << parent->EndY() << ", " << parent->EndZ() << ", "
              << parent->EndT() << " ns)"
              << "\n---------------" << std::endl;
  }

  //======================================================================
  if (iLogLevel >= kInfos)
    std::cout << "\nLOOPING over " << vp_trk.size() << " tracks..."
              << std::endl;

  // Shu: Track-level processing, 20250703---
  // Here are reconstructed tracks---
  for (art::Ptr<recob::Track> const &p_trk : vp_trk) {
    if (iLogLevel >= kInfos)
      std::cout << "trk#" << p_trk->ID() << "\r" << std::flush;

    if (p_trk->Length() < 40)
      continue; // Shu: remove short tracks (expect muon long enough),
                // 20250703---

    // Print reco track's reconstructed endpoint
    std::cout << "\tMuonRecoEnd (" << p_trk->End().X() << ", "
              << p_trk->End().Y() << ", " << p_trk->End().Z() << ") "
              << std::endl;

    // Shu: The key, all info of track particle can be acquired here,
    // 20250703--- Shu: (ChatGPT) For a certain Pandora reco track, find the
    // most likely MC truth particle--- Shu: This method is convenient, but not
    // always accurate, especially in dense or noisy events. It only gives the
    // MCParticle that contributed the most charge via BackTrackerService, not
    // necessarily the one that aligns best in 3D
    //         simb::MCParticle const * mcp =
    //         truthUtil.GetMCParticleFromRecoTrack(clockData, *p_trk, e,
    //         tag_trk.label());

    // MCPartcile finding for current RECO
    // track-------------------------------------------------------
    //  Step 1: Try truthUtil match first
    //  Retrieve single best-matching MC truth particle for the given reco track
    //  Match is based on which MCParticle contributed the most ionization
    //  energy to the hits (inherent space constraint) used to reconstruct the
    //  track (via BackTracker IDEs) Does NOT require matching of track
    //  endpoints or similar energy
    simb::MCParticle const *mcp = truthUtil.GetMCParticleFromRecoTrack(
        clockData, *p_trk, e, tag_trk.label());

    if (!mcp) {
      std::cout << "\t[No usable MC match]" << std::endl;
      continue;
    }
    if (abs(mcp->PdgCode()) !=
        13) { // If the track is not mu^- / mu^+, jump; 20250703---
      std::cout << "\t[No e match]" << std::endl;
      continue;
    }
    if (mcp->EndProcess() == "Transportation") {
      std::cout << "\t[Transportation]" << std::endl;
      continue;
    }

    auto endpoint_distance_cm = [&](const recob::Track &trk,
                                    const simb::MCParticle &tru) {
      double dxE = trk.End().X() - tru.EndX();
      double dyE = trk.End().Y() - tru.EndY();
      double dzE = trk.End().Z() - tru.EndZ();
      double dE = std::sqrt(dxE * dxE + dyE * dyE + dzE * dzE);

      double dxS =
          trk.Start().X() - tru.EndX(); // larsoft confuses start and end points
                                        // during reconstruction
      double dyS = trk.Start().Y() - tru.EndY();
      double dzS = trk.Start().Z() - tru.EndZ();
      double dS = std::sqrt(dxS * dxS + dyS * dyS + dzS * dzS);

      return std::min(dS, dE);
    };

    double dist_truthutil = endpoint_distance_cm(*p_trk, *mcp);

    if (dist_truthutil > 15.0) {
      std::cout << "\t[Far distance (Reco; True)]" << std::endl;
      continue;
    }
    //------------------------------------------------------------------------------------------------

    if (mcp->PdgCode() > 0) { // Shu: count mu^-; 20250703---
      n_mum++;
      map_mum_endproc[mcp->EndProcess()]++;
    } else {
      n_mup++;
      map_mup_endproc[mcp->EndProcess()]++;
    }

    // Michel electron
    // finding----------------------------------------------------------------------
    // Shu: For real decay of muon, at least 3 daughter particles, 20250703---
    // continue: end current track processing
    if (mcp->NumberDaughters() < 3)
      continue;

    // Not only consider last three daughter particles: In ideal physics, decay
    // particles are the last daughter particles, while in real simulation,
    // "late delta rays", 'finals-tate interactions', 'trakcing artificts' are
    // likely to be added to daughter() list after decay---
    simb::MCParticle const *mcp_mich = nullptr;

    // Tunables
    const double r2Stop = 25.0;  // (5 cm)^2: electron must be near muon stop
    const double r2ENu = 25.0;   // (5 cm)^2: e–ν spatial match
    const double tENu = 5.0;     // 5 ns:     e–ν time match
    const double KEminMeV = 1.0; // MeV:      minimum Michel KE

    auto inTPC = [&](double x, double y, double z) {
      return (x >= -356 && x <= 356) && (y >= 0 && y <= 610) &&
             (z >= 0 && z <= 465);
    };
    auto isDeltaLike = [&](const std::string &p) {
      return (p == "muIoni" || p == "hIoni" || p == "eIoni" || p == "muBrems" ||
              p == "compt" || p == "phot" || p == "conv" || p == "annihil");
    };

    // Cache muon stop position
    const double mx = mcp->EndX();
    const double my = mcp->EndY();
    const double mz = mcp->EndZ();

    // Search for best Michel electron among daughters
    double best_ke = -1.0;

    for (int i = 0; i < mcp->NumberDaughters(); ++i) {
      const simb::MCParticle *d =
          pi_serv->TrackIdToParticle_P(mcp->Daughter(i));
      if (!d)
        continue; // No decay daughter
      if (std::abs(d->PdgCode()) != 11)
        continue; // e± only

      // Electron vertex
      double x = d->Vx(), y = d->Vy(), z = d->Vz(), t = d->T();
      if (!inTPC(x, y, z))
        continue;

      // Must be close to muon stop
      double dx = x - mx, dy = y - my, dz = z - mz;
      if (dx * dx + dy * dy + dz * dz > r2Stop)
        continue;

      // Reject obvious delta-like processes
      if (isDeltaLike(d->Process()))
        continue;

      // KE cut
      double KE = (d->E() - d->Mass()) * 1e3; // MeV
      if (KE < KEminMeV)
        continue;

      // Require ≥1 neutrino sibling co-vertexed (space + time)
      bool hasNuSibling = false;
      for (int j = 0; j < mcp->NumberDaughters(); ++j) {
        const simb::MCParticle *s =
            pi_serv->TrackIdToParticle_P(mcp->Daughter(j));
        if (!s)
          continue;
        int apdg = std::abs(s->PdgCode());
        if (apdg != 12 && apdg != 14)
          continue;

        double ddx = s->Vx() - x, ddy = s->Vy() - y, ddz = s->Vz() - z;
        double ddt = std::abs(s->T() - t);
        if (ddx * ddx + ddy * ddy + ddz * ddz <= r2ENu && ddt <= tENu) {
          hasNuSibling = true;
          break;
        }
      }
      if (!hasNuSibling)
        continue;

      // Keep the highest-KE candidate
      // At this point, mcp_mich is nullptr if no Michel daughter passed the
      // cuts, or points to the best (highest KE) Michel electron/positron.
      if (KE > best_ke) {
        best_ke = KE;
        mcp_mich = d;
      }
    }

    // If not a full Michel signature, skip this muon
    if (!(mcp_mich))
      continue;

    // All reconstructed hits (recob::Hit) in the event that were primarily
    // caused by the given MCParticle
    std::vector<const recob::Hit *> v_hit_michel =
        truthUtil.GetMCParticleHits(clockData, *mcp_mich, e, tag_hit.label());
    if (mcp->EndProcess() == "muMinusCaptureAtRest") {
      if (v_hit_michel.size()) { // If we observe michel, capture here is not
                                 // true capture---
        n_cme_wh++;
        mean_cme_h += v_hit_michel.size();
      } else {
        n_cme_nh++;
      }
    }

    double lifetime_eM_us =
        (mcp_mich->T() - mcp->T()) * 1e-3; // muon birth → Michel birth

    std::cout << "\t[Michel observed!] (PDG: " << mcp_mich->PdgCode()
              << ") / mcp::TrackID: " << mcp_mich->TrackId() << " / "
              << v_hit_michel.size() << " hits" << std::endl;
    // Modified by Shu, 20250703---
    std::cout << "\tMichel True K-energy : "
              << (mcp_mich->E() - mcp_mich->Mass()) * 1e3
              << " MeV  /  lifetime [us]: " << lifetime_eM_us << std::endl;
    std::cout << "\tMichel origin: (x, y, z, t)[cm, ns]        = ("
              << mcp_mich->Vx() << ", " << mcp_mich->Vy() << ", "
              << mcp_mich->Vz() << ", " << mcp_mich->T() << ")" << std::endl;
    std::cout << "\tMCTruth muon End (x1, y1, z1, t1)[cm, ns]  = ("
              << mcp->EndX() << ", " << mcp->EndY() << ", " << mcp->EndZ()
              << ", " << mcp->EndT() << ")" << std::endl;
    std::cout << "\tMCTruth muon Start (x0, y0, z0, t0)[cm, ns]= (" << mcp->Vx()
              << ", " << mcp->Vy() << ", " << mcp->Vz() << ", " << mcp->T()
              << ")" << std::endl;
    std::cout << "\tMichel Process = " << mcp_mich->Process()
              << " / muon end process = " << mcp->EndProcess() << std::endl;

    float mich_ide_energy = 0;
    float mich_hit_energy = 0;
    for (const recob::Hit *hit_michel : v_hit_michel) {
      if (hit_michel->View() != geo::kW)
        continue;

      mich_hit_energy += hit_michel->Integral() * fADCtoMeV;

      std::vector<sim::TrackIDE> v_tid =
          bt_serv->HitToTrackIDEs(clockData, *hit_michel);
      for (const sim::TrackIDE &tid : v_tid) {
        if (tid.trackID == mcp_mich->TrackId()) {
          mich_ide_energy += tid.energy;
        }
      }
    }
  }

  if (iLogLevel >= kBasics)
    std::cout << "\033[93m"
              << "End of Truechecks::analyze "
                 "======================================================="
              << "\033[0m" << std::endl;
} // end analyze

void ana::Truechecks::beginJob() {
  if (iLogLevel >= kBasics)
    std::cout << "\033[93m"
              << "Truechecks::beginJob: "
                 "============================================================"
              << "\033[0m" << std::endl;
  if (iLogLevel >= kBasics)
    std::cout << "\033[93m"
              << "End of Truechecks::beginJob "
                 "======================================================"
              << "\033[0m" << std::endl;
} // end beginJob

void ana::Truechecks::endJob() {
  if (iLogLevel >= kBasics)
    std::cout
        << "\033[93m"
        << "Truechecks::endJob: "
           "=============================================================="
        << "\033[0m" << std::endl;

  // std::cout << "µ+ decay rate: " << 100.*n_mep / n_mup << "% (" << n_mup <<
  // ")" << std::endl; for (auto const& [key, val] : map_mup_endproc) {
  //     std::cout << "  " << key << ": " << val << std::endl;
  // }
  // std::cout << "µ- decay rate: " << 100.*n_mem / n_mum << "% (" << n_mum <<
  // ")" << std::endl; for (auto const& [key, val] : map_mum_endproc) {
  //     std::cout << "  " << key << ": " << val << std::endl;
  // }
  //    std::cout << "µ- decaying after capture: " << 100.*(n_cme_wh + n_cme_nh)
  //    / map_mum_endproc["muMinusCaptureAtRest"] << "% (" << (n_cme_wh +
  //    n_cme_nh) << ")" << std::endl; std::cout << "  w/ hits: " << n_cme_wh <<
  //    " (~" << mean_cme_h/n_cme_wh << " hits/michel)" << std::endl; std::cout
  //    << "  w/o hit: " << n_cme_nh << std::endl;

  if (iLogLevel >= kBasics)
    std::cout << "\033[93m"
              << "End of Truechecks::endJob "
                 "========================================================"
              << "\033[0m\n\n"
              << std::endl;
} // end endJob

bool ana::Truechecks::Log(bool cond, int flag, int tab, std::string msg,
                          std::string succ, std::string fail) {
  if (iLogLevel >= flag) {
    std::cout << std::string(tab, '\t') << msg << " ";
    if (cond)
      std::cout << "\033[92m" << succ << "\033[0m" << std::endl;
    else
      std::cout << "\033[91m" << fail << "\033[0m" << std::endl;
  }
  return cond;
}

std::string ana::Truechecks::GetParticleName(int pdg) {

  std::vector<std::string> periodic_table = {
      "",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
      "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc",
      "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
      "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
      "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe"};

  switch (pdg) {
  case 11:
    return "e-";
  case -11:
    return "e+";
  case 12:
    return "ve";
  case -12:
    return "-ve";
  case 13:
    return "µ-";
  case -13:
    return "µ+";
  case 14:
    return "vµ";
  case -14:
    return "-vµ";
  case 22:
    return "γ";
  case 2212:
    return "p";
  case 2112:
    return "n";
  }

  if (pdg > 1000000000) {
    unsigned ex = pdg % 10;
    unsigned A = (pdg / 10) % 1000;
    unsigned Z = (pdg / 10000) % 1000;
    unsigned L = (pdg / 10000000);
    if (L == 100 && Z && Z < periodic_table.size())
      return Form("%u%s%s", A, periodic_table[Z].c_str(), ex ? "*" : "");
  }

  return Form("%d", pdg);
}

DEFINE_ART_MODULE(ana::Truechecks)
