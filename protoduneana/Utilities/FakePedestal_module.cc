////////////////////////////////////////////////////////////////////////
// Class:       FakePedestal
// Plugin Type: analyzer (art v2_10_03)
// File:        FakePedestal_module.cc
//
// Produce fake pedestal values to test database interface
//
// Generated at Fri Mar  9 22:58:52 2018 by Tingjun Yang using cetskelgen
// from cetlib version v3_02_00.
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

#include "larcore/Geometry/WireReadout.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

#include <fstream>

namespace pdune {
  class FakePedestal;
}


class pdune::FakePedestal : public art::EDAnalyzer {
public:
  explicit FakePedestal(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FakePedestal(FakePedestal const &) = delete;
  FakePedestal(FakePedestal &&) = delete;
  FakePedestal & operator = (FakePedestal const &) = delete;
  FakePedestal & operator = (FakePedestal &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.
  const lariov::DetPedestalProvider& m_PedestalProvider;
};


pdune::FakePedestal::FakePedestal(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  m_PedestalProvider(art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider())
{}

void pdune::FakePedestal::analyze(art::Event const & e)
{
  std::ofstream outfile("pdune_fakeped.txt");
  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  for (size_t chan = 0; chan < wireReadout.Nchannels(); ++chan){
    float ped_mean    = m_PedestalProvider.PedMean(chan);
    float ped_rms     = m_PedestalProvider.PedRms(chan);
    float ped_meanerr = m_PedestalProvider.PedMeanErr(chan);
    float ped_rmserr  = m_PedestalProvider.PedRmsErr(chan);
    outfile<<chan<<" "<<ped_mean<<" "<<ped_rms<<" "<<ped_meanerr<<" "<<ped_rmserr<<std::endl;
  }
  outfile.close();
}

DEFINE_ART_MODULE(pdune::FakePedestal)
