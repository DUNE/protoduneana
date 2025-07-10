#ifndef DRAW_H
#define DRAW_H

#include <iostream>
#include <string>
#include <vector>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include "computeUtils.h"
#include "helper.h"
#include "constants.h"

namespace draw {

    void ColorMap(double Data[nR2][nHeight], TH2D* h1);

    void Extremum(std::vector<double> DataMax, std::vector<double> DataErr, std::vector<double> Volume,
        TGraphErrors* graph);
    
    void FinalExtremum(std::vector<double> Data[nR1], std::vector<double> Volume[nR1], int cutR1,
        TMultiGraph* multiGraph, TLegend* legend);
    
    void IsoVolume(const double R1, const double R2[nR2], const double Volume[nVolume], double Height[nR2],
        std::vector<TGraph*> Graphs);
} // namespace draw

#endif