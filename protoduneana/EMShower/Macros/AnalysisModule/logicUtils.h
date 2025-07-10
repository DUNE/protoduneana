#ifndef LOGIC_UTILS_H
#define LOGIC_UTILS_H

#include <vector>
#include <cmath>
#include <iostream>
#include <string>

#include "constants.h"

namespace logic_utils {
    bool particleInCone(double particleX, double particleY, double particleZ, double height, double r1, double r2);

    std::vector<bool> TrackBeforeShower(unsigned int fNParticles, double* fShowerStartX, double* fShowerStartY, double* fShowerStartZ,
    double* fTrackStartX, double* fTrackStartY, double* fTrackStartZ);
}

#endif