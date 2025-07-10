#ifndef RECO_COMPLETION_CLASS_H
#define RECO_COMPLETION_CLASS_H

#include "MCRecoClass.h"

class RecoCompletion : public Reco {
public:
double (*Completion)[nR2][nHeight];
double (*CompletionErr)[nR2][nHeight];
double (*CompletionMax)[nVolume];
double (*ReductedCompletionError)[nVolume];
std::vector<double>* VolumeFiltered;
std::vector<double>* CompletionFiltered;
std::vector<double>* CompletionErrorFiltered;

    // === Constructor ===
    RecoCompletion(std::string filePath, float electronEnergy);

    // === Destructor ===
    ~RecoCompletion() override;

    // === Methods ===
    void ComputeCompletion(int iR1, int iR2, int iHeight);
    
    void ComputeCompletionError(int iR1, int iR2, int iHeight);

    void Draw2DMap(int iR1, std::string savePath);

    void DrawMax(int iR1, std::string savePath);

    void DrawMaxFinal(std::string savePathPlot, std::string savePathData);
};

#endif // RECO_COMPLETION_CLASS_H