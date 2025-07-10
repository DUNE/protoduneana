#include "AnalysisModule/MCCompletionClass.h"
#include "AnalysisModule/RecoCompletionClass.h"

void GetCompletion_macro() {
    gROOT->SetBatch(kTRUE);
    for (int e=0; e<nEnergies; e++) {
        std::string energyStr = help::doubleToString(ElectronEnergy[e]);
        std::string filePath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/electron_"+energyStr+"GeV/analysisOutput.root";
        RecoCompletion recoCompletion(filePath, ElectronEnergy[e]);
        MCCompletion mcCompletion(filePath, ElectronEnergy[e]);

        for (int j=0; j<nR1; j++) {
            std::string r1Str = help::doubleToString(R1[j]);
            for (int k=0; k<nR2; k++) {
                for (int l=0; l<nHeight; l++) {
                    recoCompletion.Reset();
                    mcCompletion.Reset();
                    std::cout << "Computing for R1 = " << R1[j] << " cm, R2 = " << R2[k] << " cm, Height = " << bHeights[l] << " cm" << std::endl;
                    recoCompletion.ComputeCompletion(j, k, l);
                    recoCompletion.ComputeCompletionError(j, k, l);
                    mcCompletion.ComputeCompletion(j, k, l);
                    mcCompletion.ComputeCompletionError(j, k, l);
                }
            }
            std::string recoMapPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Containment"+recoCompletion.dataType+"/"+energyStr+"GeV/"+recoCompletion.dataType+"EnergyICRatio_"+r1Str+"_map";
            recoCompletion.Draw2DMap(j, recoMapPath);
            std::string recoMaxPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Containment"+recoCompletion.dataType+"/"+energyStr+"GeV/"+recoCompletion.dataType+"EnergyICRatio_"+r1Str+"_hist";
            recoCompletion.DrawMax(j, recoMaxPath);
            std::string mcMapPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Containment"+mcCompletion.dataType+"/"+energyStr+"GeV/"+mcCompletion.dataType+"EnergyICRatio_"+r1Str+"_map";
            mcCompletion.Draw2DMap(j, mcMapPath);
            std::string mcMaxPath = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Containment"+mcCompletion.dataType+"/"+energyStr+"GeV/"+mcCompletion.dataType+"EnergyICRatio_"+r1Str+"_hist";
            mcCompletion.DrawMax(j, mcMaxPath);
        }
        std::string recoFinalPlot = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Containment"+recoCompletion.dataType+"/"+recoCompletion.dataType+"FinalCompletion_"+energyStr+"GeV";
        std::string recoFinalData = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Completion/"+recoCompletion.dataType+"Completion_"+energyStr+"GeV.txt";
        recoCompletion.DrawMaxFinal(recoFinalPlot, recoFinalData);
        std::string mcFinalPlot = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/Figures/Containment"+mcCompletion.dataType+"/"+mcCompletion.dataType+"FinalCompletion_"+energyStr+"GeV";
        std::string mcFinalData = "/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/SavedData/Completion/"+mcCompletion.dataType+"Completion_"+energyStr+"GeV.txt";
        mcCompletion.DrawMaxFinal(mcFinalPlot, mcFinalData);
    }
}
