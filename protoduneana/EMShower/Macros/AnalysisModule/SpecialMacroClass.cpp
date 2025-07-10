#include "SpecialMacroClass.h"


// === SpecialMacro Constructor ===
SpecialMacro::SpecialMacro(std::string filePath, float electronEnergy) : Tree(filePath, electronEnergy) {
    // === Allocate memory for MC ===
    std::cout << "Allocate memory MC" << std::endl;
    fMCParticleEnergy = new double[NMaxMCParticles];
    fMCParticleStartPositionX = new double[NMaxMCParticles];
    fMCParticleStartPositionY = new double[NMaxMCParticles];
    fMCParticleStartPositionZ = new double[NMaxMCParticles];
    fMCParticleEndPositionX = new double[NMaxMCParticles];
    fMCParticleEndPositionY = new double[NMaxMCParticles];
    fMCParticleEndPositionZ = new double[NMaxMCParticles];

    // === Set branch addresses for MC ===
    std::cout << "SetBranchAddress MC" << std::endl;
    tree->SetBranchAddress("nMCParticles", &fNMCParticles);
    tree->SetBranchAddress("MCParticleEnergy", fMCParticleEnergy);
    tree->SetBranchAddress("MCTotalEnergy", &fMCTotalEnergy);
    tree->SetBranchAddress("MCParticleStartPositionX", fMCParticleStartPositionX);
    tree->SetBranchAddress("MCParticleStartPositionY", fMCParticleStartPositionY);
    tree->SetBranchAddress("MCParticleStartPositionZ", fMCParticleStartPositionZ);
    tree->SetBranchAddress("MCParticleEndPositionX", fMCParticleEndPositionX);
    tree->SetBranchAddress("MCParticleEndPositionY", fMCParticleEndPositionY);
    tree->SetBranchAddress("MCParticleEndPositionZ", fMCParticleEndPositionZ);

    // === Allocate memory for Reco ===
    std::cout << "Allocate memory Reco" << std::endl;
    fTrackLength = new double[NMaxParticles];
    fTrackStartX = new double[NMaxParticles];
    fTrackStartY = new double[NMaxParticles];
    fTrackStartZ = new double[NMaxParticles];
    fTrackEndX = new double[NMaxParticles];
    fTrackEndY = new double[NMaxParticles];
    fTrackEndZ = new double[NMaxParticles];
    fShowerEnergy = new double[NMaxParticles];
    fShowerStartX = new double[NMaxParticles];
    fShowerStartY = new double[NMaxParticles];
    fShowerStartZ = new double[NMaxParticles];
    fShowerEndX = new double[NMaxParticles];
    fShowerEndY = new double[NMaxParticles];
    fShowerEndZ = new double[NMaxParticles];
    fShowerLength = new double[NMaxParticles];

    // === Set branch addresses for Reco ===
    std::cout << "SetBranchAddress Reco" << std::endl;
    tree->SetBranchAddress("nParticles", &fNParticles);
    tree->SetBranchAddress("Energy", &fEnergy);
    tree->SetBranchAddress("TrackLength", fTrackLength);
    tree->SetBranchAddress("TotalTrackLength", &fTotalTrackLength);
    tree->SetBranchAddress("TrackStartX", fTrackStartX);
    tree->SetBranchAddress("TrackStartY", fTrackStartY);
    tree->SetBranchAddress("TrackStartZ", fTrackStartZ);
    tree->SetBranchAddress("TrackEndX", fTrackEndX);
    tree->SetBranchAddress("TrackEndY", fTrackEndY);
    tree->SetBranchAddress("TrackEndZ", fTrackEndZ);
    tree->SetBranchAddress("ShowerEnergy", fShowerEnergy);
    tree->SetBranchAddress("TotalShowerEnergy", &fTotalShowerEnergy);
    tree->SetBranchAddress("ShowerStartX", fShowerStartX);
    tree->SetBranchAddress("ShowerStartY", fShowerStartY);
    tree->SetBranchAddress("ShowerStartZ", fShowerStartZ);
    tree->SetBranchAddress("ShowerEndX", fShowerEndX);
    tree->SetBranchAddress("ShowerEndY", fShowerEndY);
    tree->SetBranchAddress("ShowerEndZ", fShowerEndZ);
    tree->SetBranchAddress("ShowerLength", fShowerLength);

    // === Check if track start before shower ===
    std::cout << "IsTrackBeforeShower" << std::endl;
    IsTrackBeforeShower = new std::vector<bool>[nentries];
    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        IsTrackBeforeShower[i] = logic_utils::TrackBeforeShower(fNParticles, fShowerStartX, fShowerStartY, fShowerStartZ,
            fTrackStartX, fTrackStartY, fTrackStartZ);
    }

    // === Allocate memory for BG ===
    std::cout << "Allocate memory BG" << std::endl;
    // Ar39
    fAr39Energy = new double[NMaxBGParticles];
    fAr39StartX = new double[NMaxBGParticles];
    fAr39StartY = new double[NMaxBGParticles];
    fAr39StartZ = new double[NMaxBGParticles];
    fAr39PdgCode = new unsigned int[NMaxBGParticles];
    // Ar42
    fAr42Energy = new double[NMaxBGParticles];
    fAr42StartX = new double[NMaxBGParticles];
    fAr42StartY = new double[NMaxBGParticles];
    fAr42StartZ = new double[NMaxBGParticles];
    fAr42PdgCode = new unsigned int[NMaxBGParticles];
    // Kr85
    fKr85Energy = new double[NMaxBGParticles];
    fKr85StartX = new double[NMaxBGParticles];
    fKr85StartY = new double[NMaxBGParticles];
    fKr85StartZ = new double[NMaxBGParticles];
    fKr85PdgCode = new unsigned int[NMaxBGParticles];
    // Cosmic
    fCosmicEnergy = new double[NMaxBGParticles];
    fCosmicStartX = new double[NMaxBGParticles];
    fCosmicStartY = new double[NMaxBGParticles];
    fCosmicStartZ = new double[NMaxBGParticles];
    fCosmicPdgCode = new unsigned int[NMaxBGParticles];

    // === Set branch addresses for BG ===
    std::cout << "SetBranchAddress BG" << std::endl;
    tree->SetBranchAddress("nAr39Particles", &fNAr39Particles);
    tree->SetBranchAddress("nAr42Particles", &fNAr42Particles);
    tree->SetBranchAddress("nKr85Particles", &fNKr85Particles);
    tree->SetBranchAddress("nCosmicParticles", &fNCosmicParticles);
    // Ar39
    tree->SetBranchAddress("Ar39Energy", fAr39Energy);
    tree->SetBranchAddress("Ar39StartX", fAr39StartX);
    tree->SetBranchAddress("Ar39StartY", fAr39StartY);
    tree->SetBranchAddress("Ar39StartZ", fAr39StartZ);
    tree->SetBranchAddress("Ar39PdgCode", fAr39PdgCode);
    // Ar42
    tree->SetBranchAddress("Ar42Energy", fAr42Energy);
    tree->SetBranchAddress("Ar42StartX", fAr42StartX);
    tree->SetBranchAddress("Ar42StartY", fAr42StartY);
    tree->SetBranchAddress("Ar42StartZ", fAr42StartZ);
    tree->SetBranchAddress("Ar42PdgCode", fAr42PdgCode);
    // Kr85
    tree->SetBranchAddress("Kr85Energy", fKr85Energy);
    tree->SetBranchAddress("Kr85StartX", fKr85StartX);
    tree->SetBranchAddress("Kr85StartY", fKr85StartY);
    tree->SetBranchAddress("Kr85StartZ", fKr85StartZ);
    tree->SetBranchAddress("Kr85PdgCode", fKr85PdgCode);
    // Cosmic
    tree->SetBranchAddress("CosmicEnergy", fCosmicEnergy);
    tree->SetBranchAddress("CosmicStartX", fCosmicStartX);
    tree->SetBranchAddress("CosmicStartY", fCosmicStartY);
    tree->SetBranchAddress("CosmicStartZ", fCosmicStartZ);
    tree->SetBranchAddress("CosmicPdgCode", fCosmicPdgCode);
}

// === SpecialMacro Destructor ===
SpecialMacro::~SpecialMacro() {
    std::cout << "SpecialMacro destructor" << std::endl;
    // === Delete allocated memory for MC ===
    // delete[] fMCParticleEnergy;
    // delete[] fMCParticleStartPositionX;
    // delete[] fMCParticleStartPositionY;
    // delete[] fMCParticleStartPositionZ;
    // delete[] fMCParticleEndPositionX;
    // delete[] fMCParticleEndPositionY;
    // delete[] fMCParticleEndPositionZ;

    // === Delete allocated memory for Reco ===
    // delete[] fTrackLength;
    // delete[] fTrackStartX;
    // delete[] fTrackStartY;
    // delete[] fTrackStartZ;
    // delete[] fTrackEndX;
    // delete[] fTrackEndY;
    // delete[] fTrackEndZ;
    // delete[] fShowerEnergy;
    // delete[] fShowerStartX;
    // delete[] fShowerStartY;
    // delete[] fShowerStartZ;
    // delete[] fShowerEndX;
    // delete[] fShowerEndY;
    // delete[] fShowerEndZ;
    // delete[] fShowerLength;

    // === Delete allocated memory for BG ===
    // delete[] fAr39Energy;
    // delete[] fAr39StartX;
    // delete[] fAr39StartY;
    // delete[] fAr39StartZ;
    // delete[] fAr39PdgCode;

    // delete[] fAr42Energy;
    // delete[] fAr42StartX;
    // delete[] fAr42StartY;
    // delete[] fAr42StartZ;
    // delete[] fAr42PdgCode;

    // delete[] fKr85Energy;
    // delete[] fKr85StartX;
    // delete[] fKr85StartY;
    // delete[] fKr85StartZ;
    // delete[] fKr85PdgCode;

    // delete[] fCosmicEnergy;
    // delete[] fCosmicStartX;
    // delete[] fCosmicStartY;
    // delete[] fCosmicStartZ;
    // delete[] fCosmicPdgCode;
}

void SpecialMacro::RecoEnergyIC(double* EnergyIC, double height, double r1, double r2) {
    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        double trackLengthIC = 0;
        double showerEnergyIC = 0;
        for (unsigned int j=0; j<fNParticles; ++j) {
            // Track : start and end in the cone
            if (logic_utils::particleInCone(fTrackStartX[j], fTrackStartY[j], fTrackStartZ[j], height, r1, r2)
            && logic_utils::particleInCone(fTrackEndX[j], fTrackEndY[j], fTrackEndZ[j], height, r1, r2)) {
                trackLengthIC += fTrackLength[j];
            }
            // Shower : start in the cone
            if (logic_utils::particleInCone(fShowerStartX[j], fShowerStartY[j], fShowerStartZ[j], height, r1, r2)
            && fShowerEnergy[j]>0) {
                showerEnergyIC += fShowerEnergy[j];
            }
        }
        EnergyIC[i] = showerEnergyIC+(trackLengthIC*2/1000.);
        std::cout << "RecoEnergyIC[" << i << "] : " << EnergyIC[i] << std::endl;
    }
}

void SpecialMacro::MCEnergyIC(double* EnergyIC, double height, double r1, double r2) {
    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        for (unsigned int j=0; j<fNMCParticles; j++) {
            if (logic_utils::particleInCone(fMCParticleStartPositionX[j],fMCParticleStartPositionY[j],
                fMCParticleStartPositionZ[j], height, r1, r2)) {
                EnergyIC[i] += fMCParticleEnergy[j];
            }
        }
        std::cout << "MCEnergyIC[" << i << "] : " << EnergyIC[i] << std::endl;
    }
}
