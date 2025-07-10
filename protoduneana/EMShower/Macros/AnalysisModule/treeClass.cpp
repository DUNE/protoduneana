#include "treeClass.h"

// === Tree constructor ===
Tree::Tree(std::string filePath, float electronEnergy) {
    tree = help::getTree(filePath.c_str(), "tree");
    nentries = tree->GetEntries();
    this->electronEnergy = electronEnergy;
    Volume = new double[nVolume];
    if (electronEnergy==0.3 || electronEnergy==0.5) {
        double minVol = compute_utils::VolumeTruncatedCone(sHeights[0], R1[0], R2[0]);
        double maxVol = compute_utils::VolumeTruncatedCone(sHeights[nHeight-1], R1[nR1-1], R2[nR2-1]);
        help::linspace(minVol, maxVol, nVolume, Volume);
    }
    else {
        double minVol = compute_utils::VolumeTruncatedCone(bHeights[0], R1[0], R2[0]);
        double maxVol = compute_utils::VolumeTruncatedCone(bHeights[nHeight-1], R1[nR1-1], R2[nR2-1]);
        help::linspace(minVol, maxVol, nVolume, Volume);
    }
    EnergyIC = new double[nentries];
};

// === Tree Destructor ===
Tree::~Tree() {
    std::cout << "Tree destructor" << std::endl;
    // delete tree;
    delete[] Volume;
    delete[] EnergyIC;
};

// === Reset method ===
void Tree::Reset() {
    std::fill(EnergyIC, EnergyIC + nentries, 0.0);
}