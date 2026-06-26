/* New Daphne Interface, adapted from DAPHNEInterface2 to deal with
 a vector of subdet_label required for PDVD */

#include "DAPHNEInterfaceBase.h"
#include "dunecore/DuneObj/DUNEHDF5FileInfo2.h"
#include "dunecore/HDF5Utils/HDF5RawFile3Service.h"
#include "detdataformats/daphne/DAPHNEFrame.hpp"
#include "detdataformats/daphne/DAPHNEStreamFrame.hpp"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"

#include "daqdataformats/v4_4_0/Fragment.hpp"
#include "daqdataformats/v4_4_0/SourceID.hpp"

#include "DAPHNEUtils.h"

namespace daphne {
using dunedaq::daqdataformats::SourceID;
using dunedaq::daqdataformats::FragmentType;
using dunedaq::daqdataformats::Fragment;
using dunedaq::daqdataformats::FragmentHeader;
using DAPHNEStreamFrame = dunedaq::fddetdataformats::Daphnestreamframe2;
using DAPHNEFrame = dunedaq::fddetdataformats::Daphneframe2;

class DAPHNEInterface3 : public DAPHNEInterfaceBase {

 private:
  unsigned count_warnings;

 private:
  static const size_t FragmentHeaderSize = sizeof(FragmentHeader);
  static const size_t FrameSize = sizeof(DAPHNEFrame);
  static const size_t StreamFrameSize = sizeof(DAPHNEStreamFrame);

    const detinfo::DetectorClocksData* fClocksData;

  template <class T>
  size_t GetNFrames(size_t frag_size, size_t frag_header_size) {
    return (frag_size - FragmentHeaderSize)/sizeof(T);
  }

  //Determine if we're streaming, and pick the corresponding frame type
  //and processing method
  void UnpackFragment(
      std::unique_ptr<Fragment> & frag,
      std::unordered_map<unsigned int, std::vector<raw::OpDetWaveform>> & wf_map,
      utils::DAPHNETree * daphne_tree) {

    bool is_stream = (frag->get_fragment_type() != FragmentType::kDAPHNE);

    if (!is_stream) {
      ProcessFrames(frag, wf_map, daphne_tree);
    }
    else {
      ProcessStreamFrames(frag, wf_map, daphne_tree);
    }
  }

  void ProcessFrames(
      std::unique_ptr<Fragment> & frag,
      std::unordered_map<unsigned int, std::vector<raw::OpDetWaveform>> & wf_map,
      utils::DAPHNETree * daphne_tree) {

    auto n_frames = GetNFrames<DAPHNEFrame>(frag->get_size(),
                                            FragmentHeaderSize);
    for (size_t i = 0; i < n_frames; ++i) {
      auto frame
          = reinterpret_cast<DAPHNEFrame*>(
              static_cast<uint8_t*>(frag->get_data()) + i*FrameSize);
      ProcessFrame(frame, wf_map, daphne_tree);
    }
  }

  //Get number of streaming Frames then loop over them and process each one
  void ProcessStreamFrames(
      std::unique_ptr<Fragment> & frag,
      std::unordered_map<unsigned int, std::vector<raw::OpDetWaveform>> & wf_map,
      utils::DAPHNETree * daphne_tree) {

    auto n_frames = GetNFrames<DAPHNEStreamFrame>(frag->get_size(),
                                                  FragmentHeaderSize);
    for (size_t i = 0; i < n_frames; ++i) {
      auto frame
          = reinterpret_cast<DAPHNEStreamFrame*>(
              static_cast<uint8_t*>(frag->get_data()) + i*StreamFrameSize);
      ProcessStreamFrame(frame, wf_map, daphne_tree);
    }
  }

  void ProcessStreamFrame(
      DAPHNEStreamFrame * frame,
      std::unordered_map<unsigned int, std::vector<raw::OpDetWaveform>> & wf_map,
      utils::DAPHNETree * daphne_tree) {

    art::ServiceHandle<dune::DAPHNEChannelMapService> channel_map;
    auto b_link = frame->daq_header.link_id;
    auto b_slot = frame->daq_header.slot_id;


    //Each streaming frame comes with data from 4 channels
    std::array<size_t, 4> frame_channels = {
      frame->header.channel_0,
      frame->header.channel_1,
      frame->header.channel_2,
      frame->header.channel_3};
    // Loop over channels
    for (size_t i = 0; i < frame->s_channels_per_frame; ++i) {
      auto offline_channel = -1;

      try {
        offline_channel = channel_map->GetOfflineChannel(
          b_slot, b_link, frame_channels[i]);
      }
      catch (const std::range_error & err) {
        if (count_warnings < 100) {
          std::cout << "WARNING: Could not find offline channel for " <<
                       b_slot << " " << b_link << " " << frame_channels[i] << std::endl;
          ++count_warnings;
        }
        // Don't process data from unknown offline channels!
        continue;
      }

      //Make output
      auto & waveform = daphne::utils::MakeWaveform(
            offline_channel,
            frame->s_adcs_per_channel,
            (frame->get_timestamp() & 0xffffffffff) // mask out most significant bits. Live with only 40 least significant bits ~ 1,099,511,627,776 ticks = ~4.9 h
            *fClocksData->OpticalClock().TickPeriod(), // making sure the timestamp is in microseconds; assumed it's in ticks in DAQ
            wf_map,
            true);

      // Loop over ADC values in the frame for channel i
      for (size_t j = 0; j < static_cast<size_t>(frame->s_adcs_per_channel); ++j) {
        waveform.push_back(frame->get_adc(j, i));
        if (daphne_tree != nullptr)
          daphne_tree->fADCValue[j] = frame->get_adc(j, i);
      }

      if (daphne_tree != nullptr) {
        daphne_tree->fSlot = b_slot;
        daphne_tree->fDaphneChannel = frame_channels[i];
        daphne_tree->fOfflineChannel = offline_channel;
        daphne_tree->fFrameTimestamp = frame->get_timestamp();
        daphne_tree->fTriggerSampleValue = 0;
        daphne_tree->fThreshold = 0;
        daphne_tree->fBaseline = 0;
        daphne_tree->Fill();
      }
    }
  }

  void ProcessFrame(
      DAPHNEFrame * frame,
      std::unordered_map<unsigned int, std::vector<raw::OpDetWaveform>> & wf_map,
      utils::DAPHNETree * daphne_tree) {

    art::ServiceHandle<dune::DAPHNEChannelMapService> channel_map;

    int b_channel_0 = frame->get_channel();
    int b_link = frame->daq_header.link_id;
    int b_slot = frame->daq_header.slot_id;
    auto offline_channel=-1;
    try {
       offline_channel = channel_map->GetOfflineChannel(
          b_slot, b_link, b_channel_0);
    }
    catch (const std::range_error & err) {
      //Just throw a warning so users can check out the rest of the data
      //maybe we can configure this to crash for keepup reco
      if (count_warnings < 100) {
        std::cout << "WARNING: Could not find offline channel for " <<
                     b_slot << " " << b_link << " " << b_channel_0 << std::endl;
        ++count_warnings;
      }
      // Don't process data from unknown offline channels!
      return;
    }

    //Make output waveform and fill
    auto & waveform = daphne::utils::MakeWaveform(
        offline_channel,
        static_cast<size_t>(frame->s_num_adcs),
        (frame->get_timestamp() & 0xffffffffff) // mask out most significant bits. Live with only 40 least significant bits ~ 1,099,511,627,776 ticks = ~4.9 h
        *fClocksData->OpticalClock().TickPeriod(), // making sure the timestamp is in microseconds; assumed it's in ticks in DAQ
        wf_map);
    for (size_t j = 0; j < static_cast<size_t>(frame->s_num_adcs); ++j) {
      waveform.push_back(frame->get_adc(j));
      if (daphne_tree != nullptr)
        daphne_tree->fADCValue[j] = frame->get_adc(j);
    }

    if (daphne_tree != nullptr) {
      daphne_tree->fSlot = b_slot;
      daphne_tree->fDaphneChannel = b_channel_0;
      daphne_tree->fOfflineChannel = offline_channel;
      daphne_tree->fFrameTimestamp = frame->get_timestamp();
      daphne_tree->fTriggerSampleValue = frame->header.trigger_sample_value;
      daphne_tree->fThreshold = frame->header.threshold;
      daphne_tree->fBaseline = frame->header.baseline;
      daphne_tree->Fill();
    }

  }

  bool CheckIsDetReadout(const SourceID & source_id) {
    return (source_id.subsystem == SourceID::Subsystem::kDetectorReadout);
  }

  bool CheckFragSize(std::unique_ptr<Fragment> & frag) {
    // Large enough to have header
    return (frag->get_size() > FragmentHeaderSize);
  }


 public:

  DAPHNEInterface3(fhicl::ParameterSet const& p) :
    count_warnings(0)
  {};

  void Process2(
      art::Event &evt,
      std::string inputlabel,
      std::vector<std::string> subdet_label,
      std::unordered_map<unsigned int, WaveformVector> & wf_map,
      utils::DAPHNETree * daphne_tree) override {

      // this will live through out the event being processed
    const auto clocksData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    fClocksData = &clocksData;

    //Get the HDF5 file to be opened
    auto infoHandle = evt.getHandle<raw::DUNEHDF5FileInfo2>(inputlabel);
    size_t   evtno = infoHandle->GetEvent();
    size_t   seqno = infoHandle->GetSequence();

    dunedaq::hdf5libs::HDF5RawDataFile::record_id_t record_id
        = std::make_pair(evtno, seqno);

    //Open the HDF5 file and get source ids
    art::ServiceHandle<dune::HDF5RawFile3Service> rawFileService;
    auto raw_file = rawFileService->GetPtr();
    auto source_ids = raw_file->get_source_ids(record_id);
    //Loop over source ids
    for (const auto & source_id : source_ids)  {
      // only want detector readout data (i.e. not trigger info)
      if (!CheckIsDetReadout(source_id)) continue;

      //Loop over geo ids
      auto geo_ids = raw_file->get_geo_ids_for_source_id(record_id, source_id);
      for (const auto &geo_id : geo_ids) {
        //Check that it's photon detectors
        if (!utils::CheckSubdet(geo_id, subdet_label)) continue;

        //Get the fragment
        auto frag = raw_file->get_frag_ptr(record_id, geo_id);

        // Too small to even have a header
        if (!CheckFragSize(frag)) continue;

        //Process it
        UnpackFragment(frag, wf_map, daphne_tree);
      }
    }
  };

  void Process(
      art::Event &evt,
      std::string inputlabel,
      std::string subdet_label,
      std::unordered_map<unsigned int, std::vector<raw::OpDetWaveform>> & wf_map,
      utils::DAPHNETree * daphne_tree) override {};

};
}

DEFINE_ART_CLASS_TOOL(daphne::DAPHNEInterface3)
