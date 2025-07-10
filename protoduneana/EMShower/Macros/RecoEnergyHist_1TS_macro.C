#include "AnalysisModule/helper.h"


void RecoEnergyHist_1TS_macro() {
    gROOT->SetBatch(kTRUE);
    std::string energyStr = "1";
    std::string filePath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/electron_"+energyStr+"GeV/analysisOutput.root";
    TTree *tree = help::getTree(filePath.c_str(), "tree");

    const int NMaxParticles=500;
    
    double fShowerEnergy[NMaxParticles];
    double fTrackLength[NMaxParticles];

    tree->SetBranchAddress("ShowerEnergy", &fShowerEnergy);
    tree->SetBranchAddress("TrackLength", &fTrackLength);

    Long64_t nentries = tree->GetEntries();
    double Energy[nentries];
    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        double bestShowerEnergy = 0;
        double bestTrackEnergy = 0;
        for (int j=0; j<NMaxParticles; j++) {
            if (fShowerEnergy[j] > bestShowerEnergy) {
                bestShowerEnergy = fShowerEnergy[j];
            }
            if (fTrackLength[j] > bestTrackEnergy) {
                bestTrackEnergy = fTrackLength[j];
            }
        }
        // Total energy is the sum of the best shower and track energies
        Energy[i] = bestShowerEnergy + (bestTrackEnergy * 2 / 1000.);
    }

    // ========== Create canvas ==========
    TCanvas* c1 = new TCanvas("c1", "Energy", canvaSizeX, canvaSizeY);
    // ========== Create h ==========
    double xmin = *std::minmax_element(Energy, Energy+nentries).first;
    double xmax = *std::minmax_element(Energy, Energy+nentries).second;
    TH1D* h = new TH1D("h", ";Energy (GeV);Number of events", 100, xmin, xmax);
    for (Long64_t i=0; i<nentries; i++) {
        h->Fill(Energy[i]);
    }
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    // h->SetStats(0);
    h->Draw();

    // ========== Draw the canvas ==========
    if (constants::saveData) {
        c1->SaveAs("/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/RecoEnergy_1GeV_1T&S.png");
    }
    
    // ========== Clean up ==========
    delete h;
    delete c1;
}