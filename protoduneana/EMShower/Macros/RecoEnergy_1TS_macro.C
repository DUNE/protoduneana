#include "AnalysisModule/SpecialMacroClass.h"


void RecoEnergy_1TS_macro() {
    gROOT->SetBatch(kTRUE);
    SpecialMacro tree("/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/electron_1GeV/analysisOutput.root", 1);

    double Energy[tree.nentries];
    for (Long64_t i=0; i<tree.nentries; i++) {
        tree.tree->GetEntry(i);
        double bestShowerEnergy = 0;
        double bestTrackEnergy = 0;
        for (int j=0; j<tree.NMaxParticles; j++) {
            if (tree.fShowerEnergy[j] > bestShowerEnergy) {
                bestShowerEnergy = tree.fShowerEnergy[j];
            }
            if (tree.fTrackLength[j] > bestTrackEnergy) {
                bestTrackEnergy = tree.fTrackLength[j];
            }
        }
        // Total energy is the sum of the best shower and track energies
        Energy[i] = bestShowerEnergy + (bestTrackEnergy * 2 / 1000.);
    }

    // === Create canvas ===
    TCanvas* c1 = new TCanvas("c1", "Energy", canvaSizeX, canvaSizeY);
    // === Create histogram ===
    double xmin = *std::minmax_element(Energy, Energy+tree.nentries).first;
    double xmax = *std::minmax_element(Energy, Energy+tree.nentries).second;
    TH1D* h1 = new TH1D("h1", ";Energy (GeV);Number of events", 100, xmin, xmax);
    for (Long64_t i=0; i<tree.nentries; i++) {
        h1->Fill(Energy[i]);
    }
    h1->SetLineColor(kBlue);
    h1->SetLineWidth(2);
    // h1->SetStats(0);
    h1->Draw();

    // === Draw the canvas ===
    if (saveData) {
        c1->SaveAs("/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Special/RecoEnergy_1GeV_1T&S.png");
    }
    
    // === Clean up ===
    delete h1;
    delete c1;
}