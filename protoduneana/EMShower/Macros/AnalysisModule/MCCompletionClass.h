#ifndef MC_COMPLETION_CLASS_H
#define MC_COMPLETION_CLASS_H

#include "MCRecoClass.h"

class MCCompletion : public MC {
    public:
    double (*Completion)[nR2][nHeight];
    double (*CompletionErr)[nR2][nHeight];
    double (*CompletionMax)[nVolume];
    double (*ReductedCompletionError)[nVolume];
    std::vector<double>* VolumeFiltered;
    std::vector<double>* CompletionFiltered;
    std::vector<double>* CompletionErrorFiltered;

    // === Constructor ===
    MCCompletion(std::string filePath, float electronEnergy);

    // === Destructor ===
    ~MCCompletion() override;

    // === Methods ===
    void ComputeCompletion(int iR1, int iR2, int iHeight);
    
    void ComputeCompletionError(int iR1, int iR2, int iHeight);

    void Draw2DMap(int iR1, std::string savePath);

    void DrawMax(int iR1, std::string savePath);

    void DrawMaxFinal(std::string savePathPlot, std::string savePathData);
};

#endif // MC_COMPLETION_CLASS_H