#ifndef MC_RESOLUTION_CLASS_H
#define MC_RESOLUTION_CLASS_H

#include "MCRecoClass.h"

class MCResolution : public MC {
public:
double (*Resolution)[nR2][nHeight];
double (*ResolutionErr)[nR2][nHeight];
double (*ResolutionMin)[nVolume];
double (*ReductedResolutionError)[nVolume];
std::vector<double>* VolumeFiltered;
std::vector<double>* ResolutionFiltered;
std::vector<double>* ResolutionErrorFiltered;

    // === Constructor ===
    MCResolution(std::string filePath, float electronEnergy);

    // === Destructor ===
    ~MCResolution() override;

    // === Methods ===
    void ComputeResolution(int iR1, int iR2, int iHeight);
    
    void ComputeResolutionError(int iR1, int iR2, int iHeight);

    void Draw2DMap(int iR1, std::string savePath);

    void DrawMin(int iR1, std::string savePath);

    void DrawMinFinal(std::string savePathPlot, std::string savePathData);
};

#endif // MC_RESOLUTION_CLASS_H