#include "RecoResolutionClass.h"

// === Constructor ===
RecoResolution::RecoResolution(std::string filePath, float electronEnergy) : Reco(filePath, electronEnergy) {
    Resolution = new double[nR1][nR2][nHeight];
    ResolutionErr = new double[nR1][nR2][nHeight];
    ResolutionMin = new double[nR1][nVolume];
    ReductedResolutionError = new double[nR1][nVolume];
    VolumeFiltered = new std::vector<double>[nR1];
    ResolutionFiltered = new std::vector<double>[nR1];
    ResolutionErrorFiltered = new std::vector<double>[nR1];
}

// === Destructor ===
RecoResolution::~RecoResolution() {
    // === Delete allocated memory ===
    delete[] Resolution;
    delete[] ResolutionErr;
    delete[] ResolutionMin;
    delete[] ReductedResolutionError;
    delete[] VolumeFiltered;
    delete[] ResolutionFiltered;
    delete[] ResolutionErrorFiltered;
}

// === Methods ===
void RecoResolution::ComputeResolution(int iR1, int iR2, int iHeight) {
    // This method will contain the logic to compute the Resolution for a given geometry based on the Reco data
    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        double trackLengthIC = 0.0;
        double showerEnergyIC = 0.0;
        for (unsigned int j=0; j<fNParticles; ++j) {
            // Track : start and end in the cone
            if (logic_utils::particleInCone(fTrackStartX[j], fTrackStartY[j], fTrackStartZ[j], Heights[iHeight], R1[iR1], R2[iR2]) && logic_utils::particleInCone(fTrackEndX[j], fTrackEndY[j], fTrackEndZ[j], Heights[iHeight], R1[iR1], R2[iR2]) && IsTrackBeforeShower[i][j]) {
                trackLengthIC += fTrackLength[j];
            }
            // Shower : start in the cone and ElectronEnergy>0
            if (logic_utils::particleInCone(fShowerStartX[j], fShowerStartY[j], fShowerStartZ[j], Heights[iHeight], R1[iR1], R2[iR2]) && fShowerEnergy[j]>0) {
                showerEnergyIC += fShowerEnergy[j];
            }
        }
        EnergyIC[i] = showerEnergyIC+(trackLengthIC*2/1000.);
    }
    Resolution[iR1][iR2][iHeight] += compute_utils::Resolution(EnergyIC, nentries);
    std::cout << dataType << "Resolution[" << iR1 << "][" << iR2 << "][" << iHeight << "] = " 
              << Resolution[iR1][iR2][iHeight] << " %" << std::endl;
}

void RecoResolution::ComputeResolutionError(int iR1, int iR2, int iHeight) {
    // This method will contain the logic to compute the Resolution error for a given geometry based on the Reco data
    double xmin = *std::minmax_element(EnergyIC, EnergyIC+nentries).first;
    double xmax = *std::minmax_element(EnergyIC, EnergyIC+nentries).second;
    TH1D* h = new TH1D("h", "", 100, xmin, xmax);
    for (Long64_t i=0; i<nentries; i++) {
        h->Fill(EnergyIC[i]);
    }
    double mean = h->GetMean();
    double meanErr = h->GetMeanError();
    double stddev = h->GetStdDev();
    double stddevErr = h->GetStdDevError();
    ResolutionErr[iR1][iR2][iHeight] = sqrt(std::pow(stddevErr/mean, 2)+std::pow(stddev*meanErr/std::pow(mean, 2), 2))*100;
    delete h;
}

void RecoResolution::Draw2DMap(int iR1, std::string savePath) {
    // This method will contain the logic to draw a 2D map of the Resolution for a specific R1
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
    draw::ColorMap(Resolution[iR1], h1);

    // === Drawing iso-volume curves ===
    double Height[nR2];
    std::vector<TGraph*> Graphs;
    draw::IsoVolume(R1[iR1], R2, Volume, Height, Graphs);

    // === Compute Resolution max along iso-volume lines ===
    for (int i=0; i<nVolume; i++) {
        compute_utils::IsoVolume(R1[iR1], R2, Volume[i], Height);
        int binX, binY;
        ResolutionMin[iR1][i] = compute_utils::FindHistMinAlongLine(h1, Height, R2, binX, binY);
        std::cout << "ResolutionMin[" << iR1 << "][" << i << "] = " 
                  << ResolutionMin[iR1][i] << " %" << std::endl;
        if (binX==-1 || binY==-1) {
            ReductedResolutionError[iR1][i] = 0;
            std::cout << "Invalid bin indices: binX = " << binX << ", binY = " << binY << std::endl;
            continue;
        } else ReductedResolutionError[iR1][i] = ResolutionErr[iR1][binY-1][binX-1];
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

void RecoResolution::DrawMin(int iR1, std::string savePath) {
    // This method will contain the logic to draw the max Resolution as a function of volume for each R1
    // === Create canvas ===
    TCanvas* c1 = new TCanvas("c1", "", canvaSizeX, canvaSizeY);

    // === Filter data ===
    for (int i=0; i<nVolume; i++) {
        if (ResolutionMin[iR1][i] < 100) {
            VolumeFiltered[iR1].push_back(Volume[i]);
            ResolutionFiltered[iR1].push_back(ResolutionMin[iR1][i]);
            ResolutionErrorFiltered[iR1].push_back(ReductedResolutionError[iR1][i]);
        }
    }
    
    // === Draw min Resolution ===
    TGraphErrors* graph = nullptr;
    draw::Extremum(ResolutionFiltered[iR1], ResolutionErrorFiltered[iR1], VolumeFiltered[iR1], graph);

    // === Save the plot ===
    if (saveData) {
        c1->SaveAs((savePath+".png").c_str());
        c1->SaveAs((savePath+".root").c_str());
    }

    // === Clean up ===
    delete c1;
    delete graph;
}

void RecoResolution::DrawMinFinal(std::string savePathPlot, std::string savePathData) {
    // This method will contain the logic to draw the max Resolution as a function of volume for all R1
    // === Create canvas ===
    TCanvas* c1 = new TCanvas("c1", "", canvaSizeX, canvaSizeY);
    
    // === Draw max Resolution for all R1 ===
    // Note: cutR1 is used to skip the first R1 (which is usually 0 cm)
    // This is because the first R1 is not representative of the Resolution
    // and can skew the results.
    int cutR1=0;
    TMultiGraph* multiGraph = nullptr;
    TLegend* legend = nullptr;
    draw::FinalExtremum(ResolutionFiltered, VolumeFiltered, cutR1, multiGraph, legend);

    // === Save the plot and data ===
    if (saveData) {
        c1->SaveAs((savePathPlot+".png").c_str());
        c1->SaveAs((savePathPlot+".root").c_str());

        help::saveData(savePathData, ResolutionMin);
    }

    // === Clean up ===
    delete c1;
    delete multiGraph;
    delete legend;
}