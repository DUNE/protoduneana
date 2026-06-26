#include <iostream>
#include <utility>
#include <set>

#include "art/Framework/Core/EDAnalyzer.h" 
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "art_root_io/TFileService.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h" 
#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/RawData/RawDigit.h"

#include <bitset>

namespace pdhd {

using RawDigitVector = std::vector<raw::RawDigit>;

class PDHDChannelSaturationFilter : public art::EDFilter {
  public:
    explicit PDHDChannelSaturationFilter(fhicl::ParameterSet const & pset);
    virtual ~PDHDChannelSaturationFilter() {};
    virtual bool filter(art::Event& e);

  private:
    std::string fRawDigitLabel;
    int fThresholdTicksForSaturatedChannel;
    int fThresholdChannelsForRejection;
};

PDHDChannelSaturationFilter::PDHDChannelSaturationFilter(fhicl::ParameterSet const & pset):
  EDFilter(pset),
  fRawDigitLabel(pset.get<std::string>("RawDigitLabel")),
  fThresholdTicksForSaturatedChannel(pset.get<int>("ThresholdTicksForSaturatedChannel")),
  fThresholdChannelsForRejection(pset.get<int>("ThresholdChannelsForRejection")) {}


bool PDHDChannelSaturationFilter::filter(art::Event & evt) {


  auto &rawdigits = *evt.getValidHandle<RawDigitVector>(fRawDigitLabel);
  if (rawdigits.empty()) {
    std::cout << "WARNING: no RawDigit found." << std::endl;
    return false;
  }

  geo::WireReadoutGeom const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

  const int nticks = rawdigits[0].Samples();
  const int nchans = rawdigits.size();
  int Event_ = evt.id().event();
  std::cout << "Saturation Filter evt = " << Event_ << std::endl;
  std::cout << "nticks = " << nticks << std::endl;
  std::cout << "nchans = " << nchans << std::endl;

  int n_saturated_channels = 0;
  
  for (const auto &rd : rawdigits) {
    unsigned int plane = wireReadout.ChannelToWire(rd.Channel()).front().Plane;
    if (plane < 2) continue;
    
    int n_negative_ticks_in_a_row = 0;
    for (int j = 0; j < nticks; j++) {
      float current_adc = rd.ADC(j) - rd.GetPedestal();
      
      if (current_adc < 0) {
        n_negative_ticks_in_a_row++;
      } else {
        n_negative_ticks_in_a_row = 0;
      }

      if (n_negative_ticks_in_a_row > fThresholdTicksForSaturatedChannel) {
        n_saturated_channels++;
        break;
      }
    }
  } // end of rawdigits

  std::cout << "number of saturated channels = " << n_saturated_channels << std::endl;

  // If there are enough saturated channels filter out, else pass
  if (n_saturated_channels > fThresholdChannelsForRejection) {
    std::cout << "[FILTER] Too Many Saturated Channels -> Remove." << std::endl; 
    return false;
  } else {
    return true;
  }
}



DEFINE_ART_MODULE(PDHDChannelSaturationFilter)

}
