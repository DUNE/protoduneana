#include "computeUtils.h"

namespace compute_utils {
    double Resolution(const double* Energy, int n) {
        double xmin = *std::minmax_element(Energy, Energy+n).first;
        double xmax = *std::minmax_element(Energy, Energy+n).second;
        TH1D *energyHisto = new TH1D("energyHisto", "", 100, xmin, xmax);
        
        for (int i=0; i<n; i++) {
            energyHisto->Fill(Energy[i]);
        }
    
        double mean = energyHisto->GetMean();
        double sigma = energyHisto->GetRMS();
        delete energyHisto;
        double resolution = (sigma/mean)*100;
        return resolution;
    }

    void IsoVolume(const double r1, const double R2[nR2], const double volume, double Heights[nR2]) {
        // Compute Heights knowing the volume, R1 and R2
        for (int i=0; i<nR2; ++i) {
            Heights[i] = 3*volume / (M_PI*(R2[i]*R2[i] + r1*r1 + r1*R2[i]));
        }
    }

    double FindHistMaxAlongLine(const TH2D* h, const double X[nR2], const double Y[nR2], int& maxBinX, int& maxBinY) {
        double MaxContent = std::numeric_limits<double>::min();
        maxBinX = -1;
        maxBinY = -1;
        for (int i=0; i<nR2; i++) {
            // Find the bins corresponding to the values of X, and Y
            int binX = h->GetXaxis()->FindBin(X[i]);
            int binY = h->GetYaxis()->FindBin(Y[i]);

            // Get the content of the histogram at that bin
            double content = h->GetBinContent(binX, binY);

            // Find the maximum of the content
            if (content>MaxContent) {
                MaxContent = content;
                maxBinX = binX;
                maxBinY = binY;
            }
        }
        return MaxContent;
    }

    double FindHistMinAlongLine(const TH2D* h, const double* X, const double* Y, int& minBinX, int& minBinY) {
        double MinContent = std::numeric_limits<double>::max();
        minBinX = -1;
        minBinY = -1;

        for (int i=0; i<nHeight; i++) {
            // Trouver les bins correspondant aux valeurs de X et Y
            int binX = h->GetXaxis()->FindBin(X[i]);
            int binY = h->GetYaxis()->FindBin(Y[i]);

            // Obtenir le contenu de l'histogramme à ce bin
            double content = h->GetBinContent(binX, binY);

            // Trouver le minimum du contenu
            if (content<MinContent && content!=0) {
                MinContent = content;
                minBinX = binX;
                minBinY = binY;
            }
        }
        return MinContent;
    }

    double VolumeTruncatedCone(double r1, double r2, double height) {
        return (M_PI*height/3.0) * (r1*r1 + r1*r2 + r2*r2);
    }

    double findMin2D(double array[][nHeight], int nRows, int nCols) {
        double minVal = std::numeric_limits<double>::max(); // Initialize to the largest possible value
        for (int i=0; i<nRows; ++i) {
            for (int j=0; j<nCols; ++j) {
                if (array[i][j]<minVal) {
                    minVal = array[i][j];
                }
            }
        }
        return minVal;
    }

    double findMax2D(double array[][nHeight], int nRows, int nCols) {
        double maxVal = std::numeric_limits<double>::lowest(); // Initialize to the smallest possible value
        for (int i=0; i<nRows; ++i) {
            for (int j=0; j<nCols; ++j) {
                if (array[i][j]>maxVal) {
                    maxVal = array[i][j];
                }
            }
        }
        return maxVal;
    }
}