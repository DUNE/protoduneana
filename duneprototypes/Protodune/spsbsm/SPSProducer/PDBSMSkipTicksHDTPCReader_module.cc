////////////////////////////////////////////////////////////////////////
// Class:       PDBSMSkipTicksHDTPCReader
// Plugin Type: producer (Unknown Unknown)
// File:        PDBSMSkipTicksHDTPCReader_module.cc
//
//   Module to exercise the PDHDDataInterfaceWIB3 or WIBEth tools.
//    Read raw::RawDigits into the event for a hardcoded list of APAs
//    
// Generated at Thu Nov 17 17:05:55 2022 by Thomas Junk using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dunecore/DuneObj/RDStatus.h"
#include "dunecore/DuneObj/PDSPTPCDataInterfaceParent.h"
#include "dunecore/DuneObj/DUNEHDF5FileInfo2.h"
#include "TTree.h"
#include "art_root_io/TFileService.h"

#include <memory>

class PDBSMSkipTicksHDTPCReader;


class PDBSMSkipTicksHDTPCReader : public art::EDProducer {
public:
  explicit PDBSMSkipTicksHDTPCReader(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDBSMSkipTicksHDTPCReader(PDBSMSkipTicksHDTPCReader const&) = delete;
  PDBSMSkipTicksHDTPCReader(PDBSMSkipTicksHDTPCReader&&) = delete;
  PDBSMSkipTicksHDTPCReader& operator=(PDBSMSkipTicksHDTPCReader const&) = delete;
  PDBSMSkipTicksHDTPCReader& operator=(PDBSMSkipTicksHDTPCReader&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginJob() override;

private:

  std::string m_InputLabel;
  std::string m_OutputInstance;
  std::vector<int> m_APAList;
  bool m_OutputStatusTree;
  std::unique_ptr<PDSPTPCDataInterfaceParent> m_DecoderTool;
  size_t m_SkipFirstNTPCTicks;

  TTree *m_StatusTree;
  int m_Event, m_Run, m_Subrun;
  std::vector<unsigned int> m_StatWord;

  void SetRDTSFlags(
      const std::vector<raw::RawDigit> & raw_digits,
      std::vector<raw::RDTimeStamp> & rd_timestamps);
  void FillTree(const art::Event& e,
                const std::vector<raw::RDStatus> & rdstatuscol);
};


PDBSMSkipTicksHDTPCReader::PDBSMSkipTicksHDTPCReader(fhicl::ParameterSet const& p)
  : EDProducer{p},
  m_InputLabel(p.get<std::string>("InputLabel","daq")),
  m_OutputInstance(p.get<std::string>("OutputInstance","daq")),
  m_APAList(p.get<std::vector<int>>("APAList")),
  m_OutputStatusTree(p.get<bool>("OutputStatusTree")),
  m_DecoderTool{art::make_tool<PDSPTPCDataInterfaceParent>(p.get<fhicl::ParameterSet>("DecoderToolParams"))},
  m_SkipFirstNTPCTicks(p.get<size_t>("SkipFirstNTPCTicks", 0))
{
  produces<std::vector<raw::RawDigit>>(m_OutputInstance);
  produces<std::vector<raw::RDStatus>>(m_OutputInstance);
  produces<std::vector<raw::RDTimeStamp>>(m_OutputInstance);
  produces<art::Assns<raw::RawDigit,raw::RDTimeStamp>>(m_OutputInstance);
  consumes<raw::DUNEHDF5FileInfo2>(m_InputLabel);  // the tool actually does the consuming of this product
}

void PDBSMSkipTicksHDTPCReader::FillTree(const art::Event& e,
                             const std::vector<raw::RDStatus> & rdstatuscol) {
  if (!m_OutputStatusTree) return;

  m_Run = e.run();
  m_Subrun = e.subRun();
  m_Event = e.id().event();

  m_StatWord.clear();
  for (const auto & rdstat : rdstatuscol) {
    m_StatWord.push_back(rdstat.GetStatWord());
  }

  m_StatusTree->Fill();
}

void PDBSMSkipTicksHDTPCReader::SetRDTSFlags(
    const std::vector<raw::RawDigit> & raw_digits,
    std::vector<raw::RDTimeStamp> & rd_timestamps) {
  //Needed for FEMBFilter when Raw Digits get dropped
  for (size_t i = m_SkipFirstNTPCTicks; i < raw_digits.size(); ++i) {
    rd_timestamps[i].SetFlags(raw_digits[i].Channel());
  }
}

void PDBSMSkipTicksHDTPCReader::produce(art::Event& e)
{
  std::vector<raw::RawDigit> rawdigitcol;
  std::vector<raw::RDStatus> rdstatuscol;
  std::vector<raw::RDTimeStamp> rdtscol;
  art::Assns<raw::RawDigit,raw::RDTimeStamp> rdtacol;

  m_DecoderTool->retrieveDataForSpecifiedAPAs(e, rawdigitcol, rdtscol, rdstatuscol, m_APAList);

  // Cannot have the start tick exceeding the window
  if (m_SkipFirstNTPCTicks > rawdigitcol.size()) {
    std::cout << "[WARNING] Start tick exceeds readout window. Setting start to zero." << std::endl;
    m_SkipFirstNTPCTicks = 0;
  }

  std::cout << "[INFO] Skipping first " << m_SkipFirstNTPCTicks << 
    " TPC ticks out of " << rawdigitcol.size() << std::endl; 
  
  SetRDTSFlags(rawdigitcol, rdtscol);

  if (m_OutputStatusTree) FillTree(e, rdstatuscol);

  // make associations between raw digits and RDTimestamps

  art::PtrMaker<raw::RawDigit> rdpm(e,m_OutputInstance);
  art::PtrMaker<raw::RDTimeStamp> tspm(e,m_OutputInstance);

  for (size_t i = m_SkipFirstNTPCTicks; i < rawdigitcol.size(); ++i) {
    auto const rawdigitptr = rdpm(i);
    auto const rdtimestampptr = tspm(i);
    rdtacol.addSingle(rawdigitptr,rdtimestampptr);
  }

  e.put(std::make_unique<decltype(rawdigitcol)>(std::move(rawdigitcol)),m_OutputInstance);
  e.put(std::make_unique<decltype(rdtscol)>(std::move(rdtscol)),m_OutputInstance);
  e.put(std::make_unique<decltype(rdtacol)>(std::move(rdtacol)),m_OutputInstance);
  e.put(std::make_unique<decltype(rdstatuscol)>(std::move(rdstatuscol)),m_OutputInstance);

}

void PDBSMSkipTicksHDTPCReader::beginJob() {
  if (m_OutputStatusTree) {
    art::ServiceHandle<art::TFileService> tfs;
    m_StatusTree = tfs->make<TTree>("tree","RDStatus Tree");

    m_StatusTree->Branch("statword", &m_StatWord);
    m_StatusTree->Branch("event", &m_Event);
    m_StatusTree->Branch("run", &m_Run);
    m_StatusTree->Branch("subrun", &m_Subrun);
  }
}
DEFINE_ART_MODULE(PDBSMSkipTicksHDTPCReader)
