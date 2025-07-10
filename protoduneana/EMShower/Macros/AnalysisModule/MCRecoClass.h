#ifndef MC_RECO_CLASS_H
#define MC_RECO_CLASS_H

#include <string>

#include "TTree.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TLegend.h"

#include "treeClass.h"
#include "draw.h"

class MC : public Tree {
    public:
        unsigned int fNMCParticles;
        double* fMCParticleEnergy;
        double* fMCParticleStartPositionX;
        double* fMCParticleStartPositionY;
        double* fMCParticleStartPositionZ;
        double* fMCParticleEndPositionX;
        double* fMCParticleEndPositionY;
        double* fMCParticleEndPositionZ;
        double fMCTotalEnergy;

        // === Constructor ===
        MC(std::string filePath, float electronEnergy);

        // === Destructor ===
        ~MC() override;
};

class Reco : public Tree {
    public:
        unsigned int fNParticles;
        double* fTrackLength;
        double* fTrackStartX;
        double* fTrackStartY;
        double* fTrackStartZ;
        double* fTrackEndX;
        double* fTrackEndY;
        double* fTrackEndZ;
        double* fShowerEnergy;
        double* fShowerStartX;
        double* fShowerStartY;
        double* fShowerStartZ;
        double* Heights;
        std::vector<bool>* IsTrackBeforeShower;

        // === Constructor ===
        Reco(std::string filePath, float electronEnergy);

        // === Destructor ===
        ~Reco() override;
};

#endif // MC_RECO_CLASS_H