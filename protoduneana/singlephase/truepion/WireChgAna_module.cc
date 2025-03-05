////////////////////////////////////////////////////////////////////////
// Class:       WireChgAna
// Plugin Type: analyzer (art v3_05_01)
// File:        WireChgAna_module.cc
//
// Generated at Tue Jul 21 14:37:27 2020 by Tingjun Yang using cetskelgen
// from cetlib version v3_10_00.
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
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "TTree.h"
#include "TH1D.h"

namespace pdsp {
  class WireChgAna;
}


class pdsp::WireChgAna : public art::EDAnalyzer {
public:
  explicit WireChgAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireChgAna(WireChgAna const&) = delete;
  WireChgAna(WireChgAna&&) = delete;
  WireChgAna& operator=(WireChgAna const&) = delete;
  WireChgAna& operator=(WireChgAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TTree *ftree;
  TH1D *rawwf[3][800];
  TH1D *deconwf[3][800];
  double charge[3][800];
  double deconchg[3][800];
  double hitchg[3][800];
  double hitsumadc[3][800];
  double nelec[3][800];

};


pdsp::WireChgAna::WireChgAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::WireChgAna::analyze(art::Event const& e)
{

  for (int i = 0; i<3; ++i){
    for (int j = 0; j<800; ++j){
      charge[i][j] = 0;
      nelec[i][j] = 0;
      deconchg[i][j] = 0;
      hitchg[i][j] = 0;
      hitsumadc[i][j] = 0;
    }
  }

  // Implementation of required member function here.
  art::InputTag itag("tpcrawdecoder","daq");
  auto rawListHandle = e.getHandle< std::vector<raw::RawDigit> >(itag);

  art::InputTag itag2("wclsdatanfsp","gauss");
  auto wireListHandle = e.getHandle< std::vector < recob::Wire > >(itag2);

  auto hitListHandle = e.getHandle<std::vector<recob::Hit> >("gaushit");

  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

  if (rawListHandle) {
    for (raw::RawDigit const& digit : *rawListHandle) {
      const auto & wireid = wireReadout.ChannelToWire(digit.Channel());
      if (wireid[0].TPC != 1) continue;
      //if (wireid[0].Plane != 2) continue;
      int wire = wireid[0].Wire;
      int plane = wireid[0].Plane;
      double this_charge = 0;
      std::vector<short> rawadc(6000);
      raw::Uncompress(digit.ADCs(), rawadc, digit.GetPedestal(), digit.Compression());
      for (size_t itck = 0; itck < rawadc.size(); ++itck){
        this_charge += rawadc[itck] - digit.GetPedestal();
        rawwf[plane][wire]->SetBinContent(itck, rawadc[itck] - digit.GetPedestal());
      }
      charge[plane][wire] = this_charge;
    }
  }

  auto const& simChannels = e.getProduct<std::vector<sim::SimChannel>>("tpcrawdecoder:simpleSC");
  for ( auto const& channel : simChannels ){
    const auto & wireid = wireReadout.ChannelToWire(channel.Channel());
    if (wireid[0].TPC != 1) continue;
    //if (wireid[0].Plane != 2) continue;
    int wire = wireid[0].Wire;
    int plane = wireid[0].Plane;
    double this_nelec = 0;
    auto const& timeSlices = channel.TDCIDEMap();
    for ( auto const& timeSlice : timeSlices ){
      auto const& energyDeposits = timeSlice.second;
      for ( auto const& energyDeposit : energyDeposits ){
        this_nelec += energyDeposit.numElectrons;
      }
    }
    nelec[plane][wire] = this_nelec;
  }

  if (wireListHandle) {
    for (recob::Wire const& wire : *wireListHandle) {
      const auto & wireid = wireReadout.ChannelToWire(wire.Channel());
      if (wireid[0].TPC != 1) continue;
      //if (wireid[0].Plane != 2) continue;
      int wire_num = wireid[0].Wire;
      int plane = wireid[0].Plane;
      double this_charge = 0;
      const auto & signal = wire.Signal();

      for (size_t itck = 0; itck < signal.size(); ++itck){
        this_charge += signal[itck];
        deconwf[plane][wire_num]->SetBinContent(itck, signal[itck]);
      }
      deconchg[plane][wire_num] = this_charge;
    }
  }

  if (hitListHandle) {
    for (recob::Hit const & hit : *hitListHandle){
      if (hit.WireID().TPC != 1) continue;
      int wire = hit.WireID().Wire;
      int plane = hit.WireID().Plane;
      hitchg[plane][wire] += hit.Integral();
      hitsumadc[plane][wire] += hit.ROISummedADC();
    }
  }

  ftree->Fill();
}

void pdsp::WireChgAna::beginJob()
{
  art::ServiceHandle<art::TFileService> fileServiceHandle;
  ftree = fileServiceHandle->make<TTree>("ftree", "raw digit info");
  ftree->Branch("charge", charge, "charge[3][800]/D");
  ftree->Branch("nelec", nelec, "nelec[3][800]/D");
  ftree->Branch("deconchg", deconchg, "deconchg[3][800]/D");
  ftree->Branch("hitchg", hitchg, "hitchg[3][800]/D");
  ftree->Branch("hitsumadc", hitsumadc, "hitsumadc[3][800]/D");

  for (int i = 0; i<3; ++i){
    for (int j = 0; j<800; ++j){
      rawwf[i][j] = fileServiceHandle->make<TH1D>(Form("rawwf_%d_%d",i,j),Form("Plane %d Wire %d",i,j),6000,0,6000);
      deconwf[i][j] = fileServiceHandle->make<TH1D>(Form("deconwf_%d_%d",i,j),Form("Plane %d Wire %d",i,j),6000,0,6000);
    }
  }

}

void pdsp::WireChgAna::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(pdsp::WireChgAna)
