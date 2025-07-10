#include "AnalysisModule/SpecialMacroClass.h"

void EnergyIC_macro() {
    // gROOT->SetBatch(kTRUE);
    std::cout << "Construct tree" << std::endl;
    SpecialMacro tree("/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/electron_1GeV_BGCosmics/1/analysisOutput.root", 1);

    double height = 200.;
    double r1 = 30.;
    double r2 = 60.;

    double EnergyIC[tree.nentries];
    std::cout << "Compute energyIC" << std::endl;
    tree.MCEnergyIC(EnergyIC, height, r1, r2);

    // =========== Create canvas ===========
    TCanvas* c1 = new TCanvas("c1", "Energy Completion", 1920, 1080);

    // =========== Create histogram ===========
    double minVal = *std::minmax_element(EnergyIC, EnergyIC+tree.nentries).first;
    double maxVal = *std::minmax_element(EnergyIC, EnergyIC+tree.nentries).second;
    TH1F* h1 = new TH1F("h1", ";Energy (GeV);Number of events", 100, minVal, maxVal);
    for (int i=0; i<tree.nentries; i++) {
        h1->Fill(EnergyIC[i]);
    }
    h1->SetLineColor(kBlue);
    h1->SetLineWidth(2);
    h1->Draw();
    h1->SetStats(0);
    double mean = h1->GetMean();
    double rms = h1->GetRMS();
    std::cout << "Mean: " << mean << std::endl;
    std::cout << "RMS: " << rms << std::endl;

    // =========== Save the plot ===========
    if (saveData) {
        std::string savePath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Special/RecoEnergyIC_1GeV_BGCosmics.png";
        c1->SaveAs(savePath.c_str());
    }

    // =========== Clean up ===========
    delete c1;
    delete h1;
}