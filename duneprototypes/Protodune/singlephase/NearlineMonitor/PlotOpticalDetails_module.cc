////////////////////////////////////////////////////////////////////////
//
//  PlotOpticalDetails_module.cc
//
//  Author: Leigh Whitehead (leigh.howard.whitehead@cern.ch)
//
//  Module to provide basic photon detector information for protoDUNE
//  Nearline Monitoring
//
////////////////////////////////////////////////////////////////////////

// ROOT includes
#include "TH1D.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// C++ includes
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace nlana {
 
  class PlotOpticalDetails : public art::EDAnalyzer{
  public:
 
    // Standard constructor and destructor for an ART module.
    PlotOpticalDetails(const fhicl::ParameterSet&);
    virtual ~PlotOpticalDetails();

    // This method is called once, at the start of the job. In this
    // example, it will define the histogram we'll write.
    void beginJob();

    // The analyzer routine, called once per event. 
    void analyze (const art::Event&); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fOpFlashModuleLabel;       // Input tag for OpFlash collection
   
    TH1D *hFlashTime;
    TH1D *hNFlash;
    TH1D *hNHitPerFlash;
  };

} 

//-----------------------------------------------------------------------
// Constructor
nlana::PlotOpticalDetails::PlotOpticalDetails(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
{

  // Indicate that the Input Module comes from .fcl
  fOpFlashModuleLabel = pset.get<std::string>("OpFlashLabel");

  art::ServiceHandle< art::TFileService > tfs;

  // Make a few plots
  hFlashTime    = tfs->make<TH1D>("hFlashTime","",50,-5000,5000);
  hNFlash          = tfs->make<TH1D>("hNFlash","",50,0,250);
  hNHitPerFlash    = tfs->make<TH1D>("hNHitPerFlash","",50,0,200);
}

//-----------------------------------------------------------------------
// Destructor
nlana::PlotOpticalDetails::~PlotOpticalDetails() 
{}
 
//-----------------------------------------------------------------------
void nlana::PlotOpticalDetails::beginJob()
{}

//-----------------------------------------------------------------------
void nlana::PlotOpticalDetails::analyze(const art::Event& evt) 
{
  // Get flashes from event
  auto FlashHandle = evt.getHandle< std::vector< recob::OpFlash > >(fOpFlashModuleLabel);

  // Get assosciations between flashes and hits
  art::FindManyP< recob::OpHit > Assns(FlashHandle, evt, fOpFlashModuleLabel);

  hNFlash->Fill(FlashHandle->size());

  // For every OpFlash in the vector
  for(unsigned int i = 0; i < FlashHandle->size(); ++i)
  {

    // Get OpFlash and associated hits
    art::Ptr< recob::OpFlash > TheFlashPtr(FlashHandle, i);
    recob::OpFlash TheFlash = *TheFlashPtr;

    // Find the hits associated to this flash
    std::vector< art::Ptr<recob::OpHit> > matchedHits = Assns.at(i);
    
    hFlashTime->Fill(TheFlash.Time());

    hNHitPerFlash->Fill(matchedHits.size());

  } // End loop over flashes

}

namespace nlana {
  DEFINE_ART_MODULE(PlotOpticalDetails)
}

