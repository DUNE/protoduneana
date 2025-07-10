#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "AnalysisModule/constants.h"
#include "AnalysisModule/draw.h"
#include "AnalysisModule/helper.h"

void PlotResolution_macro() {
    gROOT->SetBatch(kTRUE);
    double RecoMinResolutionFinal[nRmax1][nVolume];
    // double MCMinResolutionFinal[nRmax1][nVolume];
    double Volume[nVolume];

    for (int i=0; i<nEnergies; i++) {
        if (ElectronEnergy[i]==0.5 || ElectronEnergy[i]==0.3) {
            help::linspace(10e3, 1e6, nVolume, Volume);
        } else {
            help::linspace(21e3, 6.3e6, nVolume, Volume);
        }
        std::string energyStr = help::doubleToString(ElectronEnergy[i]);
        
        // =========== Open and read MC resolution file ===========
        // std::string mcFileName = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Resolution/MCResolution_" + energyStr + "GeV.txt";
        // std::ifstream mcFile(mcFileName);
        // if (!mcFile) {
        //     std::cerr << "Erreur d'ouverture du fichier " << mcFileName << std::endl;
        //     continue;
        // }
        // for (int r=0; r<nRmax1; r++) {
        //     for (int v=0; v<nVolume; v++) {
        //         mcFile >> MCMinResolutionFinal[r][v];
        //     }
        // }
        // mcFile.close();

        // =========== Open and read Reco resolution file ===========
        std::string recoFileName = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Resolution/RecoResolution_"+energyStr+"GeV_BGCosmics.txt";
        std::ifstream recoFile(recoFileName);
        if (!recoFile) {
            std::cerr << "Erreur d'ouverture du fichier " << recoFileName << std::endl;
            continue;
        }
        for (int r=0; r<nRmax1; r++) {
            for (int v=0; v<nVolume; v++) {
                recoFile >> RecoMinResolutionFinal[r][v];
            }
        }
        recoFile.close();

        // =========== Draw final resolution plot ===========
        // draw::FinalMinResolution("MC", energyStr, Rmax1, MCMinResolutionFinal, Volume);
        draw::FinalMinResolution("Reco", energyStr, Rmax1, RecoMinResolutionFinal, Volume);
    }
}