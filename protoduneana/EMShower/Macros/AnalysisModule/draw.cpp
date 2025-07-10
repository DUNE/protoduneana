#include "draw.h"

namespace draw {

    void ColorMap(double Data[nR2][nHeight], TH2D* h1) {
        // Draw maps of resolution and completion
        // === Find min and max values ===
        double minVal = compute_utils::findMin2D(Data, nR2, nHeight);
        double maxVal = compute_utils::findMax2D(Data, nR2, nHeight);
        h1->SetMinimum(minVal);
        h1->SetMaximum(maxVal);
        h1->SetStats(0);

        // === Fill histogram ===
        for (int i=0; i<nR2; i++) {
            for (int j=0; j<nHeight; j++) {
                if (std::isnan(Data[i][j])) {
                    std::cout << "NaN value found at index i = " << i << ", j = " << j << std::endl;
                    h1->SetBinContent(j+1, i+1, 100); // Cap the Data at 100%
                } else {
                    // std::cout << "Data[" << i << "][" << j << "] = " << Data[i][j] << std::endl;
                    h1->SetBinContent(j+1, i+1, Data[i][j]); // Note: bin counting starts at 1
                }
            }
        }

        // === Draw histogram ===
        h1->Draw("colz");
        gStyle->SetPalette(kBird);
    }

    void Extremum(std::vector<double> Data, std::vector<double> DataErr, std::vector<double> Volume,
        TGraphErrors* graph) {
        graph = new TGraphErrors(Volume.size(), Volume.data(), Data.data(), nullptr,
        DataErr.data());
        std::string title = ";Volume (cm3);Completion (%)";
        graph->SetTitle(title.c_str());
        graph->SetLineColor(kRed);
        graph->SetLineWidth(2);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.0);
        graph->Draw("ALP");
    }
    
    void FinalExtremum(std::vector<double> Data[nR1], std::vector<double> Volume[nR1], int cutR1,
            TMultiGraph* multiGraph, TLegend* legend) {
        int Colors[10] = {kOrange, kRed, kPink+1, kMagenta, kViolet+2, kBlue, kAzure+2, kCyan, kTeal+1, kGreen};
        
        // === Draw all TGraphs on the same canvas ===
        multiGraph = new TMultiGraph();
        for (int i=cutR1; i<nR1; i++) {
            TGraph* graph = new TGraph(Volume[i].size(), Volume[i].data(), Data[i].data());
            graph->SetLineColor(Colors[i]);
            graph->SetLineWidth(2);
            graph->SetMarkerSize(1.5);
            graph->SetMarkerStyle(20);
            multiGraph->Add(graph);
        }
        
        multiGraph->SetTitle(";Volume (cm3);Completion (%)");
        multiGraph->Draw("AL");
        
        // === Add a legend ===
        legend = new TLegend(0.7, 0.1, 0.9, 0.5);
        for (int i=cutR1; i<nR1; ++i) {
            std::string r1Str = help::doubleToString(R1[i]);
            legend->AddEntry(multiGraph->GetListOfGraphs()->At(i-cutR1), ("r1 = "+r1Str+" cm").c_str(), "l");
        }
        legend->Draw();
    }
    
    void IsoVolume(const double R1, const double R2[nR2], const double Volume[nVolume], double Height[nR2],
        std::vector<TGraph*> Graphs) {
        for (int i=0; i<nVolume; i++) {
            compute_utils::IsoVolume(R1, R2, Volume[i], Height);
            // Create and draw the curve
            TGraph* graph = new TGraph(nR2, Height, R2);
            graph->SetLineColor(kRed); // Red color for the curves
            graph->SetLineWidth(2); // Line thickness
            graph->SetMarkerStyle(20); // Marker style for points
            graph->SetMarkerSize(1.5); // Marker size
            graph->SetMarkerColor(kRed); // Marker color
            graph->Draw("PL SAME"); // Draw on the same canvas
            Graphs.push_back(graph);
        }
    }
} // namespace draw