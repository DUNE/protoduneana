////////////////////////////////////////////////////////////////////////
// Class:       peakBackTrack
// Plugin Type: analyzer (Unknown Unknown)
// File:        peakBackTrack_module.cc
//
// Generated at Sun Jan 25 10:54:24 2026 by Shuaixiang Zhang using cetskelgen
// from cetlib version 3.18.02.
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

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RawData/RDTimeStamp.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"


#include "nusimdata/SimulationBase/MCParticle.h"


#include <unordered_map>
#include <cmath>

#include <vector>
#include <string>
#include <limits>
#include <optional>
#include <algorithm>
#include <iomanip>


namespace michelPDHD {
  class peakBackTrack;
}


class michelPDHD::peakBackTrack : public art::EDAnalyzer {
public:
  explicit peakBackTrack(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  peakBackTrack(peakBackTrack const&) = delete;
  peakBackTrack(peakBackTrack&&) = delete;
  peakBackTrack& operator=(peakBackTrack const&) = delete;
  peakBackTrack& operator=(peakBackTrack&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // -----------------------
  // FHiCL configuration
  // -----------------------
  art::InputTag opWfTag_;
  int selectChannel_; 
  double selectWfTimeStamp_;
  int windowStartTick_;
  int windowEndTick_;

  // -----------------------
  // Helper: construct a fake OpHit from (channel, tick-window)
  // Times will be in microseconds
  // -----------------------
  recob::OpHit makeFakeOpHit_(int channel, double wfTimeStampUs, int startTick, int endTick) const;

  // convert pdg value to string name---
  std::string pdgName_(int pdg) const;

};


michelPDHD::peakBackTrack::peakBackTrack(fhicl::ParameterSet const& p)
  :EDAnalyzer{p},
   opWfTag_(p.get<art::InputTag>("OpWaveformTag", "opdec::Reco")),
   selectChannel_(p.get<int>("SelectChannel", 37)),         // REQUIRED, default value to avoid breakdown
   selectWfTimeStamp_(p.get<double>("SelectWfTimeStamp")),
   windowStartTick_(p.get<int>("WindowStartTick", 115)),     // REQUIRED
   windowEndTick_(p.get<int>("WindowEndTick", 120))         // REQUIRED
  { }



// detailed definition of makeFakeOpHit_ function
recob::OpHit michelPDHD::peakBackTrack::makeFakeOpHit_(int channel,
                                                      double wfTimeStamp,
                                                      int startTick,
                                                      int endTick) const
{
    // Absolute tick positions for the interval edges
    double const startTime = wfTimeStamp + (static_cast<double>(startTick)) * 16 * 1e-3; //us
    double const endTime = wfTimeStamp + (static_cast<double>(endTick)) * 16 * 1e-3;

    // Peak time in us
    double const peakTimeUs = (startTime + endTime) / 2.0;

    // Width in us
    double const widthUs = (static_cast<double>(endTick - startTick)) * 16 * 1e-3;

    // In this OpHit definition, we don't have explicit "start time".
    // PhotonBackTracker uses (channel, peakTime, width) as the time interval descriptor.
    double const peakTimeAbsUs = peakTimeUs; // keep consistent unless you have a true absolute reference

    unsigned short const frame = 0;
    double const area = 1.0;
    double const peakheight = 0.0;
    double const pe = 1.0;
    double const fasttototal = 0.0;

    return recob::OpHit(channel,
                        peakTimeUs,
                        peakTimeAbsUs,
                        frame,
                        widthUs,
                        area,
                        peakheight,
                        pe,
                        fasttototal);
}


std::string michelPDHD::peakBackTrack::pdgName_(int pdg) const
{
    switch (pdg) {
        case 13:   return "mu-";
        case -13:  return "mu+";
        case 11:   return "e-";
        case -11:  return "e+";
        case 22:   return "gamma";
        case 211:  return "pi+";
        case -211: return "pi-";
        case 111:  return "pi0";
        case 2212: return "p";
        case 2112: return "n";
        case 12:   return "nu_e";
        case -12:  return "nubar_e";
        case 14:   return "nu_mu";
        case -14:  return "nubar_mu";
        default:   return "PDG=" + std::to_string(pdg);
    }
}






void michelPDHD::peakBackTrack::analyze(art::Event const& e)
{
  // Fetch inputs
  auto wfListHandle = e.getHandle<std::vector<recob::OpWaveform>>(opWfTag_);

  if (wfListHandle->empty()) {
    mf::LogWarning("peakBackTrack") << "opdec::Reco wvf collection is empty; skip event.";
    return;
  }

  // PhotonBackTracker service
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;

  bool foundSelected = false;

  for (auto const& wf : *wfListHandle) {

    // Only process the requested channel
    if (static_cast<int>(wf.Channel()) != selectChannel_) continue;

    // Only process wvf with the same timestamp
    auto pdt0_us = wf.TimeStamp(); // In MC, it is automatically in unit of us---
    double tolerance_us = 0.010;
    if (std::fabs(pdt0_us - selectWfTimeStamp_) > tolerance_us) continue;

    // Sanity: ensure waveform is long enough for the requested tick window
    auto const& sig = wf.Signal();
    int const n = static_cast<int>(sig.size());
    if (windowStartTick_ < 0 || windowStartTick_ >= n || windowEndTick_ > n || windowEndTick_ <= windowStartTick_) {
      mf::LogWarning("peakBackTrack")
        << "Invalid tick window for channel " << selectChannel_
        << ": [" << windowStartTick_ << ", " << windowEndTick_
        << "), waveform size=" << n << ". Skip.";
      return;
    }

    // Sanity: to confirm the current wvf is the one we want
    std::cout << "\n[double check] Selected waveform with opch: " << wf.Channel()
              <<", and shifted timestamp: " << pdt0_us << "us\n" << std::endl;


    // Build a fake OpHit using only channel + time interval (derived from timestamps + ticks)
    recob::OpHit const fakeHit =
      makeFakeOpHit_(selectChannel_,
                     pdt0_us,
                     windowStartTick_,
                     windowEndTick_);


    // printing from OpHitToTrackIds--------------------------------------------------------------
    // Backtrack: simplest output = list of contributing G4 track IDs
    std::vector<int> const trackIds = pbt->OpHitToTrackIds(fakeHit);

    art::ServiceHandle<cheat::ParticleInventoryService> pi;

    mf::LogInfo("peakBackTrack")
        << "[OpHitToTrackIds] (Unique contributing trackIDs): n=" << trackIds.size();

    for (int const tid : trackIds) {
        simb::MCParticle const* p = pi->TrackIdToParticle_P(tid);
        std::string pname = "UNKNOWN";
        int pdg = 0;

        if (p) {
            pdg = p->PdgCode();
            pname = pdgName_(pdg);
        }

        mf::LogVerbatim("peakBackTrack")
            << "  [ID] trackID=" << tid
            << " PDG=" << pdg
            << " (" << pname << ")";
    }

    // print detailed info of mu^+/mu^-----------
    std::cout << "\nInfo of muons:" << std::endl;

    for (int const tid : trackIds) {

      simb::MCParticle const* p = pi->TrackIdToParticle_P(tid);
      if (!p) continue;

      int pdg = p->PdgCode();

      // -----------------------------------
      // Only keep mu+ / mu-
      // -----------------------------------
      if (pdg != 13 && pdg != -13) continue;

      std::string pname = pdgName_(pdg);

      // Mother
      int mother = p->Mother();

      // GEANT4 processes
      std::string startProc = p->Process();
      std::string endProc   = p->EndProcess();

      // Creation (generation) info
      double startE = p->E() * 1000;        // total energy MeV

      // End info
      const TLorentzVector& end = p->EndPosition();
      double endE = p->EndE() * 1000;       // total energy MeV

      std::cout
          << "---------------------\n"
            << "TrackID = " << tid
            << ";  PDG = " << pdg << " (" << pname << ")"
            << ";  motherID = " << mother
            << "\nstartProcess = " << startProc
            << ";  start(x,y,z,t) = ("
            << p->Vx()<< ", "
            << p->Vy() << ", "
            << p->Vz() << ", "
            << p->T() << ")"
            << ";  start E = " << startE << " MeV"
            << "\nEndProcess = " << endProc
            << ";  end(x,y,z,t) = ("
            << end.X() << ", "
            << end.Y() << ", "
            << end.Z() << ", "
            << end.T() << ")"
            << ";  end E = " << endE << " MeV\n"
            << std::endl;
    }

 

    // ---- Summarize contributions per trackID using energyFrac ---------------------------------
    std::vector<sim::TrackSDP> const trackSDPs = pbt->OpHitToTrackSDPs(fakeHit);

    // energyFrac is a fraction (0..1-ish) of this OpHit "energy"/area attributed to that track.
    // It is NOT "photon count", but it is a good "importance" metric for ranking.
    std::unordered_map<int, double> fracPerTrack;
    double totalFrac = 0.0;

    for (auto const& tsdp : trackSDPs) {
        int const tid = tsdp.trackID;
        double const f = static_cast<double>(tsdp.energyFrac);

        fracPerTrack[tid] += f;
        totalFrac += f;
    }

    // Sort by summed energyFrac
    std::vector<std::pair<int, double>> ranked;
    ranked.reserve(fracPerTrack.size());
    for (auto const& kv : fracPerTrack) {
        ranked.emplace_back(kv.first, kv.second);
    }

    std::sort(ranked.begin(), ranked.end(),
              [](auto const& a, auto const& b) {
                  return a.second > b.second;
              });

    // ---- Print summary ----
    std::cout
        << "[OpHitToTrackSDPs] run=" << e.run()
        << " subrun=" << e.subRun()
        << " event=" << e.event()
        << " ch=" << selectChannel_
        << " nTrackSDPs=" << trackSDPs.size()
        << " sumEnergyFrac=" << totalFrac << std::endl;

    
    int nPrint = 0;
    for (auto const& [tid, fsum] : ranked) {

      simb::MCParticle const* p = pi->TrackIdToParticle_P(tid);
      int pdg = 0;
      std::string pname = "UNKNOWN";
      int motherTid = 0;
      int motherPdg = 0;
      std::string motherName = "NONE";
      std::string proc = "UNKNOWN";
      double vx = 0.0, vy = 0.0, vz = 0.0, t_gen = 0.0, energy = 0.0;

      if (p) {
        pdg = p->PdgCode();
        pname = pdgName_(pdg);
        motherTid = p->Mother();
        proc = p->Process();
        vx = p->Vx();
        vy = p->Vy();
        vz = p->Vz();
        t_gen = p->T();
        energy = p->E() * 1000;

        if (motherTid != 0) {
          simb::MCParticle const* mom = pi->TrackIdToParticle_P(motherTid);
          if (mom) {
            motherPdg = mom->PdgCode();
            motherName = pdgName_(motherPdg);
          }
        }
      }

      std::cout
        << "  [Rank] trackID=" << tid
        << " PDG=" << pdg << " (" << pname << ")"
        << " mother=" << motherTid
        << " PDG=" << motherPdg << " (" << motherName << ")"
        << " process=" << proc
        << " vtx=("
        << std::fixed << std::setprecision(1)
        << vx << "," << vy << "," << vz << "," << t_gen
        << ")"
        << " E=" << std::fixed << std::setprecision(2)
        << energy << "MeV"
        << " sumEnergyFrac=" << std::setprecision(4) << fsum << std::endl;


      if (++nPrint >= 20) break;
    }
    std::cout << "\n" << std::endl;


    foundSelected = true;
    break; // since SelectChannel is required, stop after handling it
  }

  if (!foundSelected) {
    mf::LogWarning("peakBackTrack")
      << "Selected channel " << selectChannel_
      << " not found in OpWaveform collection for this event.";
  }
}




void michelPDHD::peakBackTrack::beginJob()
{
  mf::LogInfo("peakBackTrack")
    << "peakBackTrack beginJob\n"
    << "[Input of the module]\n"
    << "  SelectChannel           = " << selectChannel_ << "\n"
     << "  SelectWfTimeStamp       = " << selectWfTimeStamp_ << "\n"
    << "  WindowStartTick         = " << windowStartTick_ << "\n"
    << "  WindowEndTick           = " << windowEndTick_;
   
}



void michelPDHD::peakBackTrack::endJob()
{
  mf::LogInfo("peakBackTrack") << "peakBackTrack endJob";
}




DEFINE_ART_MODULE(michelPDHD::peakBackTrack)
