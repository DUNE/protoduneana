#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include "constants.h"

namespace help {
    std::string doubleToString(double value);

    TTree* getTree(const char* fileName, const char* treeName);

    void linspace(double start, double end, int num, double* array);

    bool ends_with(const std::string& str, const std::string& suffix);

    void saveData(const std::string& fileName, double data[nR1][nVolume]);

    bool isClose(double x1, double y1, double z1,
                 double x2, double y2, double z2,
                 double tolerance);
}

#endif