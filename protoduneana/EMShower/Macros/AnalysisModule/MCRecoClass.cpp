#include "MCRecoClass.h"

// === MC constructor ===
MC::MC(std::string filePath, float electronEnergy) : Tree(filePath, electronEnergy) {
    dataType = "MC";
    // === Allocate memory ===
    fMCParticleEnergy = new double[NMaxMCParticles];
    fMCParticleStartPositionX = new double[NMaxMCParticles];
    fMCParticleStartPositionY = new double[NMaxMCParticles];
    fMCParticleStartPositionZ = new double[NMaxMCParticles];
    fMCParticleEndPositionX = new double[NMaxMCParticles];
    fMCParticleEndPositionY = new double[NMaxMCParticles];
    fMCParticleEndPositionZ = new double[NMaxMCParticles];

    // === Set branch addresses ===
    tree->SetBranchAddress("nMCParticles", &fNMCParticles);
    tree->SetBranchAddress("MCParticleEnergy", fMCParticleEnergy);
    tree->SetBranchAddress("MCTotalEnergy", &fMCTotalEnergy);
    tree->SetBranchAddress("MCParticleStartPositionX", fMCParticleStartPositionX);
    tree->SetBranchAddress("MCParticleStartPositionY", fMCParticleStartPositionY);
    tree->SetBranchAddress("MCParticleStartPositionZ", fMCParticleStartPositionZ);
    tree->SetBranchAddress("MCParticleEndPositionX", fMCParticleEndPositionX);
    tree->SetBranchAddress("MCParticleEndPositionY", fMCParticleEndPositionY);
    tree->SetBranchAddress("MCParticleEndPositionZ", fMCParticleEndPositionZ);
};

// === MC destructor ===
MC::~MC() {
    // === Delete allocated memory ===
    delete[] fMCParticleEnergy;
    delete[] fMCParticleStartPositionX;
    delete[] fMCParticleStartPositionY;
    delete[] fMCParticleStartPositionZ;
    delete[] fMCParticleEndPositionX;
    delete[] fMCParticleEndPositionY;
    delete[] fMCParticleEndPositionZ;
};

// === Reco constructor ===
Reco::Reco(std::string filePath, float electronEnergy) : Tree(filePath, electronEnergy) {
    dataType = "Reco";
    // === Allocate memory ===
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
    
    // === Set branch addresses ===
    tree->SetBranchAddress("nParticles", &fNParticles);
    tree->SetBranchAddress("TrackLength", fTrackLength);
    tree->SetBranchAddress("TrackStartX", fTrackStartX);
    tree->SetBranchAddress("TrackStartY", fTrackStartY);
    tree->SetBranchAddress("TrackStartZ", fTrackStartZ);
    tree->SetBranchAddress("TrackEndX", fTrackEndX);
    tree->SetBranchAddress("TrackEndY", fTrackEndY);
    tree->SetBranchAddress("TrackEndZ", fTrackEndZ);
    tree->SetBranchAddress("ShowerEnergy", fShowerEnergy);
    tree->SetBranchAddress("ShowerStartX", fShowerStartX);
    tree->SetBranchAddress("ShowerStartY", fShowerStartY);
    tree->SetBranchAddress("ShowerStartZ", fShowerStartZ);
    
    Heights = new double[nHeight];
    if (electronEnergy==0.3 || electronEnergy==0.5) for (int i=0; i<nHeight; i++) Heights[i] = sHeights[i];
    else for (int i=0; i<nHeight; i++) Heights[i] = bHeights[i];

    // === Check if track start before shower ===
    IsTrackBeforeShower = new std::vector<bool>[nentries];
    for (Long64_t i=0; i<nentries; i++) {
        tree->GetEntry(i);
        IsTrackBeforeShower[i] = logic_utils::TrackBeforeShower(fNParticles, fShowerStartX, fShowerStartY, fShowerStartZ,
            fTrackStartX, fTrackStartY, fTrackStartZ);
    }
};

// === Reco destructor ===
Reco::~Reco() {
    // === Delete allocated memory ===
    delete[] fTrackLength;
    delete[] fTrackStartX;
    delete[] fTrackStartY;
    delete[] fTrackStartZ;
    delete[] fTrackEndX;
    delete[] fTrackEndY;
    delete[] fTrackEndZ;
    delete[] fShowerEnergy;
    delete[] fShowerStartX;
    delete[] fShowerStartY;
    delete[] fShowerStartZ;
    delete[] Heights;
};