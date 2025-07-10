#include "AnalysisModule/helper.h"
#include "AnalysisModule/logicUtils.h"
#include "AnalysisModule/constants.h"
#include <TGraphErrors.h>

void RecoRslCplVSEnergy_macro() {
    gROOT->SetBatch(true);
    std::string fileCompletionPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Completion/RecoFinal.txt";
    std::string fileResolutionPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Resolution/RecoFinal.txt";
    std::string fileResErrPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Resolution/RecoFInalErr.txt";
    std::string fileCplErrPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Completion/RecoFinalErr.txt";

    std::ifstream mcCompletionFile(fileCompletionPath);
    if (!mcCompletionFile) {
        std::cerr << "Erreur d'ouverture du fichier " << fileCompletionPath << std::endl;
        return;
    }
    double Completion[nEnergies];
    for (int i=0; i<nEnergies; i++) {
        mcCompletionFile >> Completion[i];
    }
    mcCompletionFile.close();

    std::ifstream mcResolutionFile(fileResolutionPath);
    if (!mcResolutionFile) {
        std::cerr << "Erreur d'ouverture du fichier " << fileResolutionPath << std::endl;
        return;
    }
    double Resolution[nEnergies];
    for (int i=0; i<nEnergies; i++) {
        mcResolutionFile >> Resolution[i];
    }
    mcResolutionFile.close();

    std::ifstream mcResErrFile(fileResErrPath);
    if (!mcResErrFile) {
        std::cerr << "Erreur d'ouverture du fichier " << fileResErrPath << std::endl;
        return;
    }
    double ResolutionErr[nEnergies];
    for (int i=0; i<nEnergies; i++) {
        mcResErrFile >> ResolutionErr[i];
    }
    mcResErrFile.close();

    std::ifstream mcCplErrFile(fileCplErrPath);
    if (!mcCplErrFile) {
        std::cerr << "Erreur d'ouverture du fichier " << fileCplErrPath << std::endl;
        return;
    }
    double CompletionErr[nEnergies];
    for (int i=0; i<nEnergies; i++) {
        mcCplErrFile >> CompletionErr[i];
    }
    mcCplErrFile.close();

    // ========== Draw Completion errorgraph ==========
    TCanvas *c1 = new TCanvas("c1", "Energy Distribution", canvaSizeX, canvaSizeY);
    TGraphErrors* g1 = new TGraphErrors(nEnergies, ElectronEnergy, Completion,
        nullptr, CompletionErr);

    g1->SetTitle(";Energy (GeV);Completion (%)");
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kBlue);
    g1->SetLineColor(kBlue);
    g1->SetLineWidth(2);
    g1->Draw("ALP");

    // ========== Save histogram ==========
    if (constants::saveData) {
        std::string savePath1 = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/RecoCompletionVSEnergy.png";
        c1->SaveAs(savePath1.c_str());
    }

    // ========== Clean up ==========
    delete c1;
    delete g1;

    // ========== Draw Resolution histogram ==========
    TCanvas *c2 = new TCanvas("c2", "Energy Distribution", canvaSizeX, canvaSizeY);
    TGraphErrors* g2 = new TGraphErrors(nEnergies, ElectronEnergy, Resolution,
        nullptr, ResolutionErr);
    g2->SetTitle(";Energy (GeV);Resolution (%)");
    g2->SetMarkerStyle(20);
    g2->SetMarkerColor(kBlue);
    g2->SetLineColor(kBlue);
    g2->SetLineWidth(2);
    g2->Draw("ALP");

    // ========== Save histogram ==========
    if (constants::saveData) {
        std::string savePath2 = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/RecoResolutionVSEnergy.png";
        c2->SaveAs(savePath2.c_str());
    }

    // ========== Clean up ==========
    delete c2;
    delete g2;
}