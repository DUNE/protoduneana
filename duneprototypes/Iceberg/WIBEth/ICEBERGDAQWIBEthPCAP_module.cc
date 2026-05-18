////////////////////////////////////////////////////////////////////////
// Class:       ICEBERGDAQWIBEthPCAP
// Plugin Type: analyzer (Unknown Unknown)
// File:        ICEBERGDAQWIBEthPCAP_module.cc
//
//   Module to emulate DAQ-formatted writing of raw::RawDigits in 
//     PCAP (Packet Capture) format
//   By Ron Rechenmacher based on FDHDDAQWIBEthBinary by Thomas Junk
//
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

#include <net/ethernet.h>		// for struct ether_header
#include <netinet/ip.h>			// struct ip
#include <netinet/udp.h>		// struct udphdr
#include <arpa/inet.h>			// ntohs, ntohl

#include <hdf5.h>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <ios>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include "detdataformats/wibeth/WIBEthFrame.hpp"
#include "detdataformats/DAQEthHeader.hpp"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "duneprototypes/Protodune/hd/ChannelMap/PD2HDChannelMapService.h"

#pragma pack(push, 1)
typedef struct pcap_hdr {
    uint32_t magic_number;
    uint16_t version_major;
    uint16_t version_minor;
    int32_t  thiszone;
    uint32_t sigfigs;
    uint32_t snaplen;
    uint32_t network;
} pcap_hdr_t;

typedef struct pcaprec_hdr {
    uint32_t ts_sec;   // timestamp seconds
    uint32_t ts_usec;  // timestamp microseconds or nanoseconds
    uint32_t incl_len; // number of bytes of packet saved in file
    uint32_t orig_len; // original length of packet (on the wire)
} pcaprec_hdr_t;

struct pkt_data_hdrs {
  struct ether_header ehdr;
  struct ip           ihdr;
  struct udphdr       uhdr;
};
#pragma pack(pop)

class ICEBERGDAQWIBEthPCAP : public art::EDAnalyzer {
public:
  explicit ICEBERGDAQWIBEthPCAP(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  ICEBERGDAQWIBEthPCAP(ICEBERGDAQWIBEthPCAP const&) = delete;
  ICEBERGDAQWIBEthPCAP(ICEBERGDAQWIBEthPCAP&&) = delete;
  ICEBERGDAQWIBEthPCAP& operator=(ICEBERGDAQWIBEthPCAP const&) = delete;
  ICEBERGDAQWIBEthPCAP& operator=(ICEBERGDAQWIBEthPCAP&&) = delete;

  void analyze(art::Event const& e) override;

private:

  std::string fOutfilenameBase;
  std::string fRawDigitLabel;
  int fCollectionPedestalOffset;
  int fInductionPedestalOffset;
  int fCrate;
};

ICEBERGDAQWIBEthPCAP::ICEBERGDAQWIBEthPCAP(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
    // More initializers here.
{
  fOutfilenameBase = p.get<std::string>("filename","icebergwibethbinary");
  fRawDigitLabel = p.get<std::string>("rawdigitlabel","tpcrawdecoder:daq");
  fCollectionPedestalOffset = p.get<int>("CollectionPedestalOffset",900);
  fInductionPedestalOffset = p.get<int>("InductionPedestalOffset",2000);
  fCrate = p.get<int>("Crate",8);
}

// Detect endianness of host at runtime
int is_little_endian() {
    uint16_t test = 0x1;
    return *(uint8_t *)&test == 0x1;
}
// Swap bytes for endian conversion
uint16_t swap16(uint16_t x) { return (x >> 8) | (x << 8); }

uint32_t swap32(uint32_t x) {
    return  (x >> 24) |
           ((x >> 8) & 0x0000FF00) |
           ((x << 8) & 0x00FF0000) |
           (x << 24);
}

void ICEBERGDAQWIBEthPCAP::analyze(art::Event const& e)
{
  art::ServiceHandle<dune::PD2HDChannelMapService> electronicsMap;

  auto runno = e.run();
  auto subrun = e.subRun();
  auto evtno = e.event();

  bool warnedNegative = false;  // warn just once per event

  // construct output filename from fOutfilenameBase_run<nnnnnn>_sr<nnn>_event<nnnnnn>.dat

  std::string trgname = "_run";
  std::ostringstream ofm1;
  ofm1 << fOutfilenameBase << "_run" << std::internal << std::setfill('0') << std::setw(6) << runno;
  ofm1 << "_sub" << std::internal << std::setfill('0') << std::setw(3) << subrun;
  ofm1 << "_evt" << std::internal << std::setfill('0') << std::setw(6) << evtno << ".pcap";
  std::string outputfilename = ofm1.str();
  FILE *ofile = fopen(outputfilename.c_str(),"w");

  auto hdr = pcap_hdr();
  hdr.magic_number = 0xa1b2c3d4;
  hdr.version_major = 2;
  hdr.version_minor = 4;
  hdr.thiszone = 0;
  hdr.sigfigs = 0;
  hdr.snaplen = 262144;
  hdr.network = 1;

  // the file will be majority little endian...
  if (!is_little_endian()) {
    // Convert all fields to host byte order
    hdr.version_major = swap16(hdr.version_major);
    hdr.version_minor = swap16(hdr.version_minor);
    hdr.thiszone      = (int32_t)swap32((uint32_t)hdr.thiszone);
    hdr.sigfigs       = swap32(hdr.sigfigs);
    hdr.snaplen       = swap32(hdr.snaplen);
    hdr.network       = swap32(hdr.network);
  }

  fwrite(&hdr, sizeof(hdr), 1, ofile);

  // this will throw an exception if the raw digits cannot be found.

  auto const& RawDigits = e.getProduct< std::vector<raw::RawDigit> >(fRawDigitLabel);

  // need a map of indices to raw digits -- does not need to be sorted
  // check in the same loop that all raw digits have the same number of samples

  std::unordered_map<uint32_t,uint32_t> rdmap;
  size_t nSamples = 0;
  for (uint32_t iptn=0; iptn<RawDigits.size(); ++iptn)
    {
      rdmap[RawDigits[iptn].Channel()] = iptn;
      size_t nSc = RawDigits[iptn].Samples();
      if (nSamples == 0)
	{
	  nSamples = nSc;
	}
      if (nSamples != nSc)
	{
	  throw cet::exception("ICEBERGDAQWIBEthPCAP") << "raw digits have different numbers of samples: " 
							 << nSamples << " " <<  nSc;
	}
    }
  std::vector<short> uncompressed(nSamples);

  // WIBEth frames save 64 channels' worth of data in 64 time samples for each frame
  // hard-code these for now.  dunedaq::fddetdataformats::WIBEthFrame::s_num_channels
  // dunedaq::fddetdataformats::WIBEthFrame::s_time_samples_per_frame  are both 64

  // Loop over the detector components -- if we do not have all the channels' raw digits to make a 64-channel WIB frame,
  // just skip making it.

  uint32_t chanlist[64];  // a list of channel IDs to go into the frame
  int planelist[64];      // which plane each of these channels is in

  int det_id = 3;  // HD_TPC
  int crate = fCrate;  // no need to loop over crates
  for (int slot=0; slot<5; ++slot)
    {
      for (int streamindex=0; streamindex<8; ++streamindex)
	{
	  // make a list of raw digit indices for the 64 channels that go into this stream's frames
	  int stream = streamindex < 4 ? streamindex : streamindex+60;    // 0, 1, 2, 3, 64, 65, 66, 67
	  int locstream = stream & 0x3;
	  int link = (stream >> 6) & 1;
	  bool skip=false;
	  for (int streamchan=0; streamchan<64; ++streamchan)
	    {
	      int wibframechan = streamchan + 64*locstream; 
	      auto cinfo = electronicsMap->GetChanInfoFromWIBElements(crate,slot,link,wibframechan);
	      if (!cinfo.valid)  // shouldn't happen -- map should have everything we can loop over
		{
		  skip = true;
		  break;
		  throw cet::exception("ICEBERGDAQWIBEthPCAP") << "electronics map does not have "  << crate << " " << slot << " " << stream << " " << streamchan << "\n";
		}
	      auto cmi = rdmap.find(cinfo.offlchan);
	      if (cmi == rdmap.end())  // MC may not have all channels.  Skip entire stream if a channel in it is missing
		{
		  skip = true;
		  break;
		}
	      chanlist[streamchan] = cmi->second;   // raw digit index
	      planelist[streamchan] = cinfo.plane;
	    } // end loop over stream chan to fill chanlist
	  if (skip) continue;  // didn't find all the channels we need to make this frame
	      
	  dunedaq::fddetdataformats::WIBEthFrame frame;
	  std::vector<std::vector<short>> uncompressed(64, std::vector<short>(nSamples));
	  for (int streamchan=0; streamchan<64; ++streamchan)
	    {
	      int ch = chanlist[streamchan];
	      int pedestal = (int) (RawDigits[ch].GetPedestal() + 0.5); // nearest integer
	      raw::Uncompress(RawDigits[ch].ADCs(),uncompressed.at(streamchan),pedestal,RawDigits[ch].Compression());
	    }

	  int iseq = 0;
	  for (size_t isample=0; isample<nSamples; ++isample)
	    {
	      int isf = isample % 64;
	      for (int streamchan=0; streamchan<64; ++streamchan)
		{
		  int d = uncompressed.at(streamchan).at(isample);
		  int pedestaloffset = (planelist[streamchan] == 2) ? fCollectionPedestalOffset : fInductionPedestalOffset;
		  d += pedestaloffset;
		  if (d<0)
		    {
		      d = 0;
		      if (!warnedNegative)
			{
			  MF_LOG_WARNING("ICEBERGDAQWIBEthPCAP_module") << "Negative ADC value in raw::RawDigit.  Setting to zero to put in WIB frame\n";
			  warnedNegative = true;
			}
		    }
		  uint16_t d16 = d & 0xFFFF;
		  frame.set_adc(streamchan,isf,d16);
		}
		  
	      if ( (isample+1) % 64 == 0)
		{
		  frame.daq_header.det_id = det_id;
		  frame.daq_header.crate_id = crate;
		  frame.daq_header.slot_id = slot;
		  frame.daq_header.stream_id = stream;
		  frame.daq_header.seq_id = iseq;
		  frame.daq_header.block_length = 899;
		  frame.daq_header.timestamp = iseq*2048;
		  // maybe some better values here for WIBEthHeader.
		  frame.header.colddata_timestamp_0 = iseq*2048;
		  frame.header.pad_0 = 0;
		  frame.header.colddata_timestamp_1 = 0;
		  frame.header.pad_1 = 0;
		  frame.header.cd = 0;
		  frame.header.crc_err = 0;
		  frame.header.link_valid = 1;
		  frame.header.lol = 0;
		  frame.header.wib_sync = 1;
		  frame.header.femb_sync = 1;
		  frame.header.pulser = 0;
		  frame.header.calibration = 0;
		  frame.header.ready = 1;
		  frame.header.context = 0;
		  frame.header.version = 0;
		  frame.header.channel = 0;
		  frame.header.extra_data = 0;

		  pcaprec_hdr_t       pchdr;
		  struct timeval tv;
		  gettimeofday(&tv,NULL);
		  pchdr.ts_sec = tv.tv_sec;
		  pchdr.ts_usec = tv.tv_usec;
		  pchdr.incl_len =
		    sizeof(pkt_data_hdrs)
		    + sizeof(frame)+ 1;   // ethernet trailer in practice is 1 byte
		  pchdr.orig_len = pchdr.incl_len;
		  fwrite(&pchdr, sizeof(pchdr), 1, ofile);
		  printf("incl_len=%d pkt_data_hdrs=%zd frame=%zd\n",
		    	 pchdr.incl_len, sizeof(pkt_data_hdrs), sizeof(frame) );

		  struct pkt_data_hdrs hdrs;
		  unsigned char dst[ETH_ALEN] = {0xff,0xff,0xff,0xff,0xff,0xff};
		  unsigned char src[ETH_ALEN] = {0x00,0x11,0x22,0x33,0x44,0x55};
		  memcpy(&(hdrs.ehdr.ether_dhost),&dst,sizeof(dst));
		  memcpy(&(hdrs.ehdr.ether_shost),&src,sizeof(src));
		  hdrs.ehdr.ether_type = htons(ETH_P_IP);
		  hdrs.ihdr.ip_hl	  = 5;
		  hdrs.ihdr.ip_v	  = 4;
		  hdrs.ihdr.ip_tos	  = 0;
		  hdrs.ihdr.ip_len	  = htons(7228);
		  hdrs.ihdr.ip_id	  = htons(56064);
		  hdrs.ihdr.ip_off	  = htons(0);
		  hdrs.ihdr.ip_ttl	  = 128;
		  hdrs.ihdr.ip_p	  = 17;
		  hdrs.ihdr.ip_sum	  = htons(0);
		  hdrs.ihdr.ip_src.s_addr = htonl(0xc0a80101);
		  hdrs.ihdr.ip_dst.s_addr = htonl(0xc0a80102);
		  hdrs.uhdr.source = htons(2048);
		  hdrs.uhdr.dest   = htons(2049);
		  hdrs.uhdr.len    = htons(7208);
		  hdrs.uhdr.check  = htons(0);
		  fwrite(&hdrs, sizeof(hdrs), 1, ofile);

		  fwrite(&frame, sizeof(frame), 1, ofile);
		  unsigned char trailer_byte = 0xc0;
		  fwrite(&trailer_byte, sizeof(trailer_byte), 1, ofile);
		  ++iseq;
		}
	    }
	      
	} // end loop over stream index
    } // end loop over slots
  fclose(ofile);
} // end analyze method

DEFINE_ART_MODULE(ICEBERGDAQWIBEthPCAP)
