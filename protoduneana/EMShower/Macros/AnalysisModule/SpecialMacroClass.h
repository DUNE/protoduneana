#ifndef SPECIAL_MACRO_CLASS_H
#define SPECIAL_MACRO_CLASS_H

#include "treeClass.h"

class SpecialMacro : public Tree {
    public:
        // === MC ===
        unsigned int fNMCParticles;
        double fMCTotalEnergy;
        double* fMCParticleEnergy;
        double* fMCParticleStartPositionX;
        double* fMCParticleStartPositionY;
        double* fMCParticleStartPositionZ;
        double* fMCParticleEndPositionX;
        double* fMCParticleEndPositionY;
        double* fMCParticleEndPositionZ;
        unsigned int* fMCParticlePdgCode;

        // === Reco ===
        unsigned int fNParticles;
        double fEnergy;
        double* fTrackLength;
        double fTotalTrackLength;
        double* fTrackStartX;
        double* fTrackStartY;
        double* fTrackStartZ;
        double* fTrackEndX;
        double* fTrackEndY;
        double* fTrackEndZ;
        double* fShowerEnergy;
        double fTotalShowerEnergy;
        double* fShowerStartX;
        double* fShowerStartY;
        double* fShowerStartZ;
        double* fShowerEndX;
        double* fShowerEndY;
        double* fShowerEndZ;
        double* fShowerLength;
        std::vector<bool>* IsTrackBeforeShower;

        // === Background info ===
        // Ar39
        unsigned int fNAr39Particles;
        double* fAr39Energy;
        double* fAr39StartX;
        double* fAr39StartY;
        double* fAr39StartZ;
        unsigned int* fAr39PdgCode;
        // Ar42
        unsigned int fNAr42Particles;
        double* fAr42Energy;
        double* fAr42StartX;
        double* fAr42StartY;
        double* fAr42StartZ;
        unsigned int* fAr42PdgCode;
        // Kr85
        unsigned int fNKr85Particles;
        double* fKr85Energy;
        double* fKr85StartX;
        double* fKr85StartY;
        double* fKr85StartZ;
        unsigned int* fKr85PdgCode;
        // Cosmic
        unsigned int fNCosmicParticles;
        double* fCosmicEnergy;
        double* fCosmicStartX;
        double* fCosmicStartY;
        double* fCosmicStartZ;
        unsigned int* fCosmicPdgCode;

        // === Constructor ===
        SpecialMacro(std::string filePath, float electronEnergy);
        // === Avoid duplicate ===
        SpecialMacro(const SpecialMacro&) = delete;
        SpecialMacro& operator=(const SpecialMacro&) = delete;
        // === Destructor ===
        ~SpecialMacro() override;

        // === Methods ===
        void RecoEnergyIC(double* EnergyIC, double height, double r1, double r2);

        void MCEnergyIC(double* EnergyIC, double height, double r1, double r2, std::vector<int> pdgCode);
};

#endif // SPECIAL_MACRO_CLASS_H