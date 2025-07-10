#ifndef TREE_CLASS_H
#define TREE_CLASS_H

#include <string>

#include "TTree.h"

#include "constants.h"
#include "helper.h"
#include "computeUtils.h"
#include "logicUtils.h"

class Tree {
    public:
        TTree *tree;
        Long64_t nentries;
        float electronEnergy;
        int NMaxParticles = 1000;
        int NMaxMCParticles = 900000;
        int NMaxBGParticles = 7000;
        double* Volume;
        double* EnergyIC;
        std::string dataType;

        // === Constructor ===
        Tree(std::string filePath, float electronEnergy);

        // === Destructor ===
        virtual ~Tree();

        // === Methods ===
        void Reset();

};

#endif