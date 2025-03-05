////////////////////////////////////////////////////////////////////////
// Class:       HitNormCheck
// Plugin Type: analyzer (art v3_04_00)
// File:        HitNormCheck_module.cc
//
// Generated at Tue Feb 18 11:11:40 2020 by Vyacheslav Galymov using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <algorithm>
#include <iterator>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "art_root_io/TFileService.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" 
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"

// ROOT
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"

// boost
#include <boost/format.hpp> 

using std::cerr;
using std::cout;
using std::endl;
using std::string;

namespace test {
  class HitNormCheck;
}


class test::HitNormCheck : public art::EDAnalyzer {

public:
  explicit HitNormCheck(fhicl::ParameterSet const& p);
  
  // Plugins should not be copied or assigned.
  HitNormCheck(HitNormCheck const&) = delete;
  HitNormCheck(HitNormCheck&&) = delete;
  HitNormCheck& operator=(HitNormCheck const&) = delete;
  HitNormCheck& operator=(HitNormCheck&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  int fLogLevel;
  std::string fHitModuleLabel;
  std::string fWireModuleLabel;
  int fPedSamples;
  int fPadRawRight;
  int fMaxPadRoi;
  int fWirePadLeft;
  int fWirePadRight;
  int fDump;
  //int fHitSepTicks;
  
  typedef struct hitinfo_t
  {
    unsigned evenum;
    unsigned hitid;
    unsigned chan;
    int  start;
    int  end;
    int  multiplicity;
    float  peak;
    float  hitint;
    float  adcsum;
    float  adcsumpad;
    float  rawsum;
    float  rawped;
    float  rawpedrms;
    float  rawsumnosig;
  } hitinfo_t ;
  
  float Pedestal( const raw::RawDigit::ADCvector_t &digi, 
		  //const std::vector<float> &wire,
		  int last_hitend,
		  int this_hitstart );
  
  float PedestalRMS( const raw::RawDigit::ADCvector_t &digi, 
		  int last_hitend, int this_hitstart );

  float PulseSumNoSignal( const raw::RawDigit::ADCvector_t &digi, 
			  int last_hitend, int this_hitstart, float pedestal );

  void PedestalProperties( const raw::RawDigit::ADCvector_t &digi, 
			   const std::vector<float> &wire, 
			   int this_hitstart,
			   hitinfo_t &hinfo );

  float PulseSum( const raw::RawDigit::ADCvector_t &digi, 
		  int pstart, int pend, float pedestal );
  
  template< typename Iterator >
  void DumpRoi( const std::string roi_type, 
		unsigned eve, unsigned chan, unsigned hit, 
		int roi_start, float ped,
		Iterator begin, Iterator end );
  
  string makeHistoName( const std::string roi_type, 
			unsigned evid, unsigned chid, unsigned hitid ){
    boost::format fmt( "%s_ev%06u_ch%04u_hit%05u" );
    fmt % roi_type; fmt % evid; fmt % chid; fmt % hitid;
    return fmt.str();
  }
  
  
  //
  TTree   *fTree;
  unsigned fEveNum;
  unsigned fChanId;
  unsigned fHitId;
  int      fStartTick;
  int      fEndTick;
  int      fMultiplicity;
  float    fPeak;
  float    fAdcSum;
  float    fAdcSumPad;
  float    fHitInt;
  float    fRawSum;
  float    fRawPed;
  float    fRawPedRMS;
  float    fRawSumNoSig;
};


test::HitNormCheck::HitNormCheck(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fLogLevel( p.get< int >("LogLevel") ),
  fHitModuleLabel( p.get< std::string  >("HitModuleLabel") ),
  fWireModuleLabel( p.get< std::string  >("WireModuleLabel") ),
  fPedSamples( p.get< int >("PedSamples") ),
  fPadRawRight( p.get< int >("PadRawRight") ),
  fMaxPadRoi( p.get< int >("MaxPadRoi") ),
  fWirePadLeft( p.get< int >("WirePadLeft") ),
  fWirePadRight( p.get< int >("WirePadRight") ),
  fDump( p.get< int >("Dump") )
  
  {}

//
void test::HitNormCheck::analyze(art::Event const& e)
{
  const string myname = "test::HitNormCheck::analyze: ";

  // get hits ValidHandle< std::vector<recob::Hits> >
  auto Hits  = e.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);

  // get wires
  auto Wires = e.getValidHandle<std::vector<recob::Wire>>(fWireModuleLabel);
  
  if( fLogLevel >= 2 ){
    cout<<myname<<"The event contains "<< Hits->size() <<" hits\n";
  }

  // get wires associated with hits
  art::FindOneP<recob::Wire>   HitToWires(Hits, e, fHitModuleLabel);
  art::FindOneP<raw::RawDigit> WireToDigits(Wires, e, fWireModuleLabel);

  //
  int dumped = 0;
  //TFile *fout = NULL;
  //if( fDump > 0 ) {
  //fout = TFile::Open("dump_hit_rois.root", "RECREATE");
  //}
  
  //
  // main loop
  raw::ChannelID_t curChannel = raw::InvalidChannelID;
  unsigned iSkip  = 0;
  auto wireIt     = Wires->begin();
  //int last_hitend = 0;
  
  std::vector< hitinfo_t > hitInfoVec;

  //for (const recob::Hit& hit: *Hits) {
  for( size_t hitIter = 0; hitIter < Hits->size(); hitIter++ ){
    
    art::Ptr<recob::Hit>   hit(Hits, hitIter);
    
    if( fLogLevel >= 3 ){
	cout<<myname<<"Hit "<<hitIter<<*hit;
    }
    
    //
    hitinfo_t hInfo;
    hInfo.evenum = e.id().event();
    hInfo.chan   = hit->Channel();
    hInfo.hitid  = hitIter;
    hInfo.start  = hit->StartTick();
    hInfo.end    = hit->EndTick();
    hInfo.multiplicity = hit->Multiplicity();
    hInfo.adcsum  = hit->ROISummedADC();
    hInfo.hitint  = hit->Integral();
    hInfo.peak    = hit->PeakAmplitude();
    
    if( hInfo.multiplicity > 1 ){
      if( fLogLevel >= 3 ){
	cout<<myname<<"Process multi hit group"<<endl;
      }
      size_t hitIdx = 1;
      while( hitIdx < (unsigned)hInfo.multiplicity ){
	hitIter++;
	art::Ptr<recob::Hit> multiHit(Hits, hitIter);
	if( multiHit->LocalIndex() != (int)hitIdx ){
	  cerr<<myname<<"Hit index is wrong: "<<multiHit->LocalIndex() <<" != "<<hitIdx <<endl;
	  cout<<hitIter<<" "<<hitIdx<<" "<<hInfo.multiplicity<<" "<<multiHit->LocalIndex()<<" "<<multiHit->Channel()<<" "<<hit->Channel()<<" "<<hInfo.chan<<endl;
	  break;
	}
	hInfo.hitint += multiHit->Integral();
	hitIdx++;
      }
      
      
    }

    // multi-hits
    if( hit->LocalIndex() > 0 ) {
      cout<<myname<<"This should not happen"<<endl;
      continue;
    }


    // get recob::Wire for this hit
    auto wire = HitToWires.at(hitIter);
    
    // get recob::Wire zero padded to full length
    auto signal = wire->Signal();
    int wire_start = hit->StartTick() - fWirePadLeft;
    if( wire_start < 0 ) wire_start = 0;
    int wire_end  = hit->EndTick() + fWirePadRight;
    if( wire_end > (int)signal.size() ) wire_end = (int)signal.size() - 1;
    hInfo.adcsumpad = std::accumulate( signal.begin() + wire_start, 
				       signal.begin() + wire_end, 0. );
    //if( fLogLevel >= 3 ){
    //cout<<myname<<"Signal size "<<signal.size()<<" summed ADC in hit window "<<sumADC<<endl;
    //}
    
    // get raw digit
    raw::ChannelID_t wireChannel = wire->Channel();
    if( wireChannel != hit->Channel() ){
      cerr<<myname<<"ERROR mismatch in hit and wire channel"<<endl;
      break;
    }

    // loop over wires
    for( ;wireIt != Wires->end(); ++wireIt )
      {
	if( wireIt->Channel() == wireChannel ){
	  //cout<<"found it "<<endl; 
	  break;
	}
      }

    if (wireIt == Wires->end()){
      cerr<<myname<<"Could not find recob::Wire for channel for the hit"<<wireChannel<<endl;
      break;
    }
    
    auto digit = WireToDigits.at( std::distance(Wires->begin(), wireIt) );
    if( digit->Channel() != wireChannel ) {
      cout<<myname<<"There is a channel mismatch with the RawDigit "<<wireChannel
	  <<" "<<digit->Channel()<<endl;
    }
    
    // check compression
    if( digit->Compression() != raw::kNone ){
      cerr<<myname<<"ERROR: compression flag is set on raw data"<<endl;
      break;
    }

    // get ADCs
    auto adc = digit->ADCs();

    //
    if( curChannel != wireChannel ){
      curChannel  = wireChannel;
      //last_hitend = 0;
    }
    
    // calculate pedestal for the raw data
    //hInfo.rawped      = Pedestal( adc, last_hitend, hit->StartTick() );
    //hInfo.rawpedrms   = PedestalRMS( adc, last_hitend, hit->StartTick() );
    //hInfo.rawsumnosig = PulseSum( adc, last_hitend, hit->StartTick(), hInfo.rawped );
    //last_hitend     = hit->EndTick();
    PedestalProperties(adc, signal, hit->StartTick(), hInfo );
        
    //
    if( hInfo.rawped < 0 ){
      if( fLogLevel >= 2 ){ 
	cout<<myname<<"Could not calculate pedestal for the hit group "<<endl;
      }
    }
    
    if( fLogLevel >= 3 ){
      cout<<myname<<"Pedestal for this hit group in channel "
	  <<curChannel<<" : "<<hInfo.rawped<<" ADC"<<endl;
    }
    
    int pstart = hit->StartTick();
    int pend   = hit->EndTick() + fPadRawRight;
    if( pend > (int)digit->NADC() ) {
      if( fLogLevel >= 2 ){ 
	cout<<myname<<"Skipping hit group at the end of the readout window"<<endl;
      }
      ++iSkip; 
      continue;
    }
    
    hInfo.rawsum = PulseSum( adc, pstart, pend, hInfo.rawped );
    if( fLogLevel >= 3 ){
      cout<<myname<<"Raw pulse sum "<<fRawSum << " ADC"<<endl;
    }

    if( dumped < fDump ){
      //float tmp = (hInfo.adcsum - hInfo.rawsum)/(hInfo.rawsum);
      if( hInfo.rawped > 0 ){ //&& tmp < -0.5){
	// not all ROIs are correct though as we will remove some which are too close later
	DumpRoi<decltype(signal.begin())>( "wire", hInfo.evenum, 
					   hInfo.chan, hInfo.hitid, 
					   pstart, 0, 
					   signal.begin() + pstart,
					   signal.begin() + hit->EndTick() );
	
	DumpRoi<decltype(adc.begin())>( "digi", hInfo.evenum, 
					hInfo.chan, hInfo.hitid, 
					pstart, hInfo.rawped, 
					adc.begin() + pstart,
					adc.begin() + pend );
	dumped++;
      }
    }
    
    hitInfoVec.push_back( hInfo );

    //
    //if( tmp < -0.4 && hInfo.rawped > 0 )
    //cout<<myname<<"to check "<<tmp<<" "<<hInfo.rawsum<<*hit;
    
  }// end looping over hits

    
  // use only hits that do not follow each other too closely
  // and then for which we calculated pedestal
  int hitSep = fPadRawRight * 2;
  for( auto it = hitInfoVec.begin(); it!=hitInfoVec.end(); ++it )
    {
      auto nx = std::next( it, 1 );
      
      int dt = hitSep;
      if( nx != hitInfoVec.end() )
	{
	  if( it->chan == nx->chan ){
	    dt = nx->start - it->start;
	  }
	}
      if( dt < hitSep )	{
	if( fLogLevel >= 2 ){
	  cout<<myname<<"Skip hits that are separated "<<dt<<" ticks"<<endl;
	}
	++iSkip;
	continue;
      }

      if( it->rawped < 0 ) {
	++iSkip; continue;
      }

      //cout<<it->chan<<" / "<<nx->chan<<" "<<it->end<<" "<<nx->start<<" "<<dt<<endl;
      
      //
      fEveNum    = it->evenum;
      fChanId    = it->chan;
      fHitId     = it->hitid;
      fStartTick = it->start;
      fEndTick   = it->end;
      fMultiplicity = it->multiplicity;
      fAdcSum    = it->adcsum;
      fAdcSumPad = it->adcsumpad;
      fHitInt    = it->hitint;
      fRawSum    = it->rawsum;
      fRawPed    = it->rawped;
      fRawPedRMS = it->rawpedrms;
      fPeak      = it->peak;
      fRawSumNoSig = it->rawsumnosig;
      fTree->Fill();
    }
  
  if( fLogLevel >= 1 ){
    cout<<myname<<"Processed "<<Hits->size()<<" hits and skipped "<<iSkip<<endl;
  }
  
}

void test::HitNormCheck::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("hitCheck","Check hit charge");
  fTree->Branch("fEveNum", &fEveNum, "fEveNum/i");
  fTree->Branch("fChanId", &fChanId, "fChanId/i");
  fTree->Branch("fHitId", &fHitId, "fHitId/i");
  fTree->Branch("fStartTick", &fStartTick, "fStartTick/I");
  fTree->Branch("fEndTick", &fEndTick, "fEndTick/I");
  fTree->Branch("fMultiplicity", &fMultiplicity, "fMultiplicity/I");
  fTree->Branch("fPeak",   &fPeak,   "fPeak/F");
  fTree->Branch("fHitInt", &fHitInt, "fHitInt/F");
  fTree->Branch("fAdcSum", &fAdcSum, "fAdcSum/F");
  fTree->Branch("fAdcSumPad", &fAdcSumPad, "fAdcSumPad/F");
  fTree->Branch("fRawSum", &fRawSum, "fRawSum/F");
  fTree->Branch("fRawSumNoSig", &fRawSumNoSig, "fRawSumNoSig/F");
  fTree->Branch("fRawPed", &fRawPed, "fRawPed/F");
  fTree->Branch("fRawPedRMS", &fRawPedRMS, "fRawPedRMS/F");

}

void test::HitNormCheck::endJob()
{}

// calculate pedestal
float test::HitNormCheck::Pedestal( const raw::RawDigit::ADCvector_t &digi,
				    int last_hitend, int this_hitstart )
{
  int istart = last_hitend;
  if( last_hitend > 0 ) istart += 2*fPadRawRight;
  
  // don't calculate anything if overlapping with the last hit
  int delta = this_hitstart - istart;
  if( delta < fPedSamples )
    return -999;
  
  int pedstart  = this_hitstart - fPedSamples;
  return TMath::Mean( digi.begin() + pedstart, digi.begin() + this_hitstart );
}

// calculate pedestal RMS
float test::HitNormCheck::PedestalRMS( const raw::RawDigit::ADCvector_t &digi,
				       int last_hitend, int this_hitstart )
{
  int istart = last_hitend;
  if( last_hitend > 0 ) istart += 2*fPadRawRight;
  
  // don't calculate anything if overlapping with the last hit
  int delta = this_hitstart - istart;
  if( delta < fPedSamples )
    return -999;
  
  int pedstart  = this_hitstart - fPedSamples;
  return TMath::RMS( digi.begin() + pedstart, digi.begin() + this_hitstart );
}

//
float test::HitNormCheck::PulseSumNoSignal( const raw::RawDigit::ADCvector_t &digi,
					    int last_hitend, int this_hitstart, float pedestal )
{
  int istart = last_hitend;
  if( last_hitend > 0 ) istart += 2*fPadRawRight;
  
  // don't calculate anything if overlapping with the last hit
  int delta = this_hitstart - istart;
  if( delta < fPedSamples )
    return -9999;
  
  if( pedestal < 0 )
    return -9999;

  int pedstart  = this_hitstart - (int)(1.5 * fPadRawRight);
  int pedend    = pedstart + fPadRawRight;
  float sum = std::accumulate(digi.begin() + pedstart, digi.begin() + pedend, 0.);
  int nsa = pedend - pedstart;
  sum -= (pedestal * nsa);
  return sum;
}

//
// the recob::Wire assumes zero suppression of pedestal
void test::HitNormCheck::PedestalProperties(  const raw::RawDigit::ADCvector_t &digi, 
					      const std::vector<float> &wire, 
					      int this_hitstart,
					      hitinfo_t &hinfo )
{
  const string myname = "test::HitNormCheck::PedestalProperties: ";
  hinfo.rawped      = -999;
  hinfo.rawpedrms   = 0.;
  hinfo.rawsumnosig = 0.;
  
  if( wire.size() != digi.size() ){
    cerr<<myname<<"size mismatch"<<endl;
    return;
  }
  
  int pedend = -1;
  for( int i = this_hitstart; i>=0; --i)
    {
      if( wire[i] == 0 ) {
	pedend = i;
	break;
      }
    }
  if( pedend < 0 ) return;
  
  // too many samples to skip
  if( (this_hitstart - pedend) > fMaxPadRoi ) return;
  
  int pedstart = -1;
  for( int i= pedend; i>=0; --i )
    {
      if( wire[i] != 0 ) break;
      pedstart = i;
    }

  if( pedstart < 0 ) return;
  
  if( (pedend - pedstart) < fPedSamples ) return;
  
  hinfo.rawped    = TMath::Mean( digi.begin() + pedstart, digi.begin() + pedend );
  hinfo.rawpedrms = TMath::RMS( digi.begin() + pedstart, digi.begin() + pedend );

  // calculate sum in no signal window
  float sum = std::accumulate(digi.begin() + pedstart, digi.begin() + pedstart + fPadRawRight, 0.);
  sum -= (hinfo.rawped * fPadRawRight);
  hinfo.rawsumnosig = sum;
  
  return;
}




// calculate pulse sum in raw data window
float test::HitNormCheck::PulseSum( const raw::RawDigit::ADCvector_t &digi, 
				    int pstart, int pend, float pedestal )

{
  float sum = std::accumulate(digi.begin() + pstart, digi.begin() + pend, 0.);
  int nsa = pend - pstart;
  sum -= (pedestal * nsa);
  return sum;
}


template< typename Iterator >  
void test::HitNormCheck::DumpRoi( const std::string roi_type,
				  unsigned eve, unsigned chan, unsigned hit, 
				  int roi_start, float ped, 
				  Iterator begin, Iterator end )
{
  art::ServiceHandle<art::TFileService> fs;
  
  int nbins = std::distance(begin, end);
  string hname = makeHistoName( roi_type, eve, chan, hit );
  TH1F *h = fs->make<TH1F>( hname.c_str(), hname.c_str(), nbins, 
			    roi_start, roi_start + nbins );
  h->GetXaxis()->SetTitle("ticks"); 
  //h->SetDirectory( 0 );
  int ibin = 1;
  for( auto it = begin; it != end; ++it ){
    h->SetBinContent( ibin, *it - ped );
    ibin++;
  }
  
  //fout->cd();
  //h->Write();
  //h->Delete();
}


DEFINE_ART_MODULE(test::HitNormCheck)
