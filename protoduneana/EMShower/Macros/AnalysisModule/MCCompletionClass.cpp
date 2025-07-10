#include "MCCompletionClass.h"

// === Constructor ===
MCCompletion::MCCompletion(std::string filePath, float electronEnergy) : MC(filePath, electronEnergy) {
    Completion = new double[nR1][nR2][nHeight];
    CompletionErr = new double[nR1][nR2][nHeight];
    CompletionMax = new double[nR1][nVolume];
    ReductedCompletionError = new double[nR1][nVolume];
    VolumeFiltered = new std::vector<double>[nR1];
    CompletionFiltered = new std::vector<double>[nR1];
    CompletionErrorFiltered = new std::vector<double>[nR1];
}

// === Destructor ===
MCCompletion::~MCCompletion() {
    // === Delete allocated memory ===
    delete[] Completion;
    delete[] CompletionErr;
    delete[] CompletionMax;
    delete[] ReductedCompletionError;
    delete[] VolumeFiltered;
    delete[] CompletionFiltered;
    delete[] CompletionErrorFiltered;
}

// === Methods ===
void MCCompletion::ComputeCompletion(int iR1, int iR2, int iHeight) {
    // This method will contain the logic to compute the completion for a given geometry based on the MC data
    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        for (unsigned int j=0; j<fNMCParticles; ++j) {
            if (logic_utils::particleInCone(fMCParticleStartPositionX[j], fMCParticleStartPositionY[j],
                fMCParticleStartPositionZ[j], bHeights[iHeight], R1[iR1], R2[iR2])) {
                EnergyIC[i] += fMCParticleEnergy[j];
            }
        }
        Completion[iR1][iR2][iHeight] += (EnergyIC[i])/(electronEnergy*nentries)*100.0;
    }
    std::cout << dataType << "Completion[" << iR1 << "][" << iR2 << "][" << iHeight << "] = " 
              << Completion[iR1][iR2][iHeight] << " %" << std::endl;
}

void MCCompletion::ComputeCompletionError(int iR1, int iR2, int iHeight) {
    // This method will contain the logic to compute the completion error for a given geometry based on the MC data
    double xmin = *std::minmax_element(EnergyIC, EnergyIC+nentries).first;
    double xmax = *std::minmax_element(EnergyIC, EnergyIC+nentries).second;
    TH1D* h = new TH1D("h", "", 100, xmin, xmax);
    for (Long64_t i=0; i<nentries; i++) {
        h->Fill(EnergyIC[i]);
    }
    double meanErr = h->GetMeanError();
    CompletionErr[iR1][iR2][iHeight] = (meanErr/electronEnergy)*100.;
    delete h;
}

void MCCompletion::Draw2DMap(int iR1, std::string savePath) {
    // This method will contain the logic to draw a 2D map of the completion for a specific R1
    // === Create canvas ===
    TCanvas* c1 = new TCanvas("c1", "", canvaSizeX, canvaSizeY);
    c1->SetRightMargin(0.15); // Make room for color scale
    // === Create histogram ===
    std::string histTitle = ";Height (cm);R2 (cm)";
    TH2D* h1 = nullptr;
    if (electronEnergy==0.3 || electronEnergy==0.5) {
        h1 = new TH2D("h1", histTitle.c_str(), nHeight, sHeights[0], sHeights[nHeight-1]+1, nR2, R2[0], R2[nR2-1]+1);
    } else {
        h1 = new TH2D("h1", histTitle.c_str(), nHeight, bHeights[0], bHeights[nHeight-1]+1, nR2, R2[0], R2[nR2-1]+1);
    }
    
    // === Draw 2D map ===
    draw::ColorMap(Completion[iR1], h1);

    // === Drawing iso-volume curves ===
    double Height[nR2];
    std::vector<TGraph*> Graphs;
    draw::IsoVolume(R1[iR1], R2, Volume, Height, Graphs);

    // === Compute completion max along iso-volume lines ===
    for (int i=0; i<nVolume; i++) {
        compute_utils::IsoVolume(R1[iR1], R2, Volume[i], Height);
        int binX, binY;
        CompletionMax[iR1][i] = compute_utils::FindHistMaxAlongLine(h1, Height, R2, binX, binY);
        std::cout << "CompletionMax[" << iR1 << "][" << i << "] = " 
                  << CompletionMax[iR1][i] << " %" << std::endl;
        if (binX==-1 || binY==-1) {
            ReductedCompletionError[iR1][i] = 0;
            std::cout << "Invalid bin indices: binX = " << binX << ", binY = " << binY << std::endl;
            continue;
        } else ReductedCompletionError[iR1][i] = CompletionErr[iR1][binY-1][binX-1];
    }

    // === Save the plot ===
    if (saveData) {
        c1->SaveAs((savePath+".png").c_str());
        c1->SaveAs((savePath+".root").c_str());
    }

    // === Clean up ===
    delete c1;
    delete h1;
    for (auto graph : Graphs) {
        if (graph) delete graph;
    }
}

void MCCompletion::DrawMax(int iR1, std::string savePath) {
    // This method will contain the logic to draw the max completion as a function of volume for each R1
    // === Create canvas ===
    TCanvas* c1 = new TCanvas("c1", "", canvaSizeX, canvaSizeY);

    // === Filter data ===
    for (int i=0; i<nVolume; i++) {
        if (CompletionMax[iR1][i] > 10) {
            VolumeFiltered[iR1].push_back(Volume[i]);
            CompletionFiltered[iR1].push_back(CompletionMax[iR1][i]);
            CompletionErrorFiltered[iR1].push_back(ReductedCompletionError[iR1][i]);
        }
    }

    // === Draw max completion ===
    TGraphErrors* graph = nullptr;
    draw::Extremum(CompletionFiltered[iR1], CompletionErrorFiltered[iR1], VolumeFiltered[iR1], graph);

    // === Save the plot ===
    if (saveData) {
        c1->SaveAs((savePath+".png").c_str());
        c1->SaveAs((savePath+".root").c_str());
    }

    // === Clean up ===
    delete c1;
    delete graph;
}

void MCCompletion::DrawMaxFinal(std::string savePathPlot, std::string savePathData) {
    // This method will contain the logic to draw the max completion as a function of volume for all R1
    // === Create canvas ===
    TCanvas* c1 = new TCanvas("c1", "", canvaSizeX, canvaSizeY);
    
    // === Draw max completion for all R1 ===
    // Note: cutR1 is used to skip the first R1 (which is usually 0 cm)
    // This is because the first R1 is not representative of the completion
    // and can skew the results.
    int cutR1=0;
    TMultiGraph* multiGraph = nullptr;
    TLegend* legend = nullptr;
    draw::FinalExtremum(CompletionFiltered, VolumeFiltered, cutR1, multiGraph, legend);

    // === Save the plot and data ===
    if (saveData) {
        c1->SaveAs((savePathPlot+".png").c_str());
        c1->SaveAs((savePathPlot+".root").c_str());

        help::saveData(savePathData, CompletionMax);
    }

    // === Clean up ===
    delete c1;
    delete multiGraph;
    delete legend;
}