#include "AnalysisModule/helper.h"
#include "AnalysisModule/SpecialMacroClass.h"


void EnergyVSLength_macro() {
    gROOT->SetBatch(kTRUE);
    std::string filePath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/electron_1GeV/analysisOutput.root";
    SpecialMacro tree(filePath, 1.0);

    double TotalEnergy[tree.nentries];
    double TotalTrackLength[tree.nentries];
    for (Long64_t i=0; i<tree.nentries; i++) {
        tree.tree->GetEntry(i);
        TotalEnergy[i] = *(tree.fEnergy);
        TotalTrackLength[i] = *(tree.fTotalTrackLength);
    }

    // === Filter data ===
    std::vector<double> filteredEnergy;
    std::vector<double> filteredTrackLength;
    for (Long64_t i=0; i<tree.nentries; i++) {
        if (TotalTrackLength[i] != 0) {
            filteredEnergy.push_back(TotalEnergy[i]);
            filteredTrackLength.push_back(TotalTrackLength[i]);
        }
    }

    // === Create TH2D ===
    int nBinsX = 100; // Number of bins for track length
    int nBinsY = 100; // Number of bins for energy
    double xMin = 0;  // Minimum track length
    double xMax = 300; // Maximum track length (adjust as needed)
    double yMin = 0;  // Minimum energy
    double yMax = 1.6; // Maximum energy (adjust as needed)

    TH2D* h = new TH2D("h", ";Track Length (cm);Energy (GeV)", nBinsX, xMin, xMax, nBinsY, yMin, yMax);

    // Fill the histogram
    for (size_t i = 0; i < filteredTrackLength.size(); i++) {
        h->Fill(filteredTrackLength[i], filteredEnergy[i]);
    }

    // === Create canvas ===
    TCanvas* c1 = new TCanvas("c1", "Energy vs Track Length", 1920, 1080);

    // Draw the histogram
    h->Draw("COLZ");
    h->SetStats(0);

    // === Save the plot ===
    if (saveData) {
        std::string savePath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/test.png";
        c1->SaveAs(savePath.c_str());
    }

    // === Clean up ===
    delete c1;
    delete h;
}