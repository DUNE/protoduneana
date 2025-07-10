#ifndef COMPUTE_UTILS_H
#define COMPUTE_UTILS_H

#include <string>
#include <iostream>
#include <TH1D.h>
#include <TH2D.h>
#include "constants.h"

namespace compute_utils {
    double Resolution(const double* Energy, int n);

    void IsoVolume(const double r1, const double R2[nR2], const double volume, double Heights[nR2]);

    double FindHistMaxAlongLine(const TH2D* h, const double* X, const double* Y, int& maxBinX, int& maxBinY);

    double FindHistMinAlongLine(const TH2D* h, const double* X, const double* Y, int& minBinX, int& minBinY);

    double VolumeTruncatedCone(double r1, double r2, double height);

    double findMin2D(double array[][nHeight], int nRows, int nCols);

    double findMax2D(double array[][nHeight], int nRows, int nCols);
}

#endif