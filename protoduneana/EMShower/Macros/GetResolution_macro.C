#include "AnalysisModule/MCResolutionClass.h"
#include "AnalysisModule/RecoResolutionClass.h"

void GetResolution_macro() {
    gROOT->SetBatch(kTRUE);
    for (int e=0; e<nEnergies; e++) {
        std::string energyStr = help::doubleToString(ElectronEnergy[e]);
        std::string filePath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/electron_"+energyStr+"GeV/analysisOutput.root";
        RecoResolution recoResolution(filePath, ElectronEnergy[e]);
        MCResolution mcResolution(filePath, ElectronEnergy[e]);

        for (int j=0; j<nR1; j++) {
            std::string r1Str = help::doubleToString(R1[j]);
            for (int k=0; k<nR2; k++) {
                for (int l=0; l<nHeight; l++) {
                    recoResolution.Reset();
                    mcResolution.Reset();
                    std::cout << "Computing for R1 = " << R1[j] << " cm, R2 = " << R2[k] << " cm, Height = " << bHeights[l] << " cm" << std::endl;
                    recoResolution.ComputeResolution(j, k, l);
                    recoResolution.ComputeResolutionError(j, k, l);
                    mcResolution.ComputeResolution(j, k, l);
                    mcResolution.ComputeResolutionError(j, k, l);
                }
            }
            std::string recoMapPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Resolution"+recoResolution.dataType+"/"+energyStr+"GeV/"+recoResolution.dataType+"ResolutionIC_"+r1Str+"_map";
            recoResolution.Draw2DMap(j, recoMapPath);
            std::string recoMaxPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Resolution"+recoResolution.dataType+"/"+energyStr+"GeV/"+recoResolution.dataType+"ResolutionIC_"+r1Str+"_hist";
            recoResolution.DrawMin(j, recoMaxPath);
            std::string mcMapPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Resolution"+mcResolution.dataType+"/"+energyStr+"GeV/"+mcResolution.dataType+"ResolutionIC_"+r1Str+"_map";
            mcResolution.Draw2DMap(j, mcMapPath);
            std::string mcMaxPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Resolution"+mcResolution.dataType+"/"+energyStr+"GeV/"+mcResolution.dataType+"ResolutionIC_"+r1Str+"_hist";
            mcResolution.DrawMin(j, mcMaxPath);
        }
        std::string recoFinalPlot = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Resolution"+recoResolution.dataType+"/"+recoResolution.dataType+"FinalResolution_"+energyStr+"GeV";
        std::string recoFinalData = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Resolution/"+recoResolution.dataType+"Resolution_"+energyStr+"GeV.txt";
        recoResolution.DrawMinFinal(recoFinalPlot, recoFinalData);
        std::string mcFinalPlot = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Resolution"+mcResolution.dataType+"/"+mcResolution.dataType+"FinalResolution_"+energyStr+"GeV";
        std::string mcFinalData = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Resolution/"+mcResolution.dataType+"Resolution_"+energyStr+"GeV.txt";
        mcResolution.DrawMinFinal(mcFinalPlot, mcFinalData);
    }
}

