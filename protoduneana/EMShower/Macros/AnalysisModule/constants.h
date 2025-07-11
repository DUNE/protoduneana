#ifndef CONSTANTS_H
#define CONSTANTS_H

constexpr int nHeight = 3;
constexpr int nR2 = 2;
constexpr int nR1 = 1;
constexpr int nVolume = 10;
constexpr int nEnergies = 4;

constexpr double ElectronEnergy[nEnergies] = {0.3, 0.5, 0.75, 1.0};
constexpr double R1[nR1] = {40.0};
constexpr double R2[nR2] = {70.0, 100.0};
constexpr double sHeights[nHeight] = {20.0};
constexpr double bHeights[nHeight] = {200.0, 300.0, 400.0};

constexpr int canvaSizeX = 1920;
constexpr int canvaSizeY = 1080;

constexpr bool saveData = true;

constexpr double beamX = 94.8; // beam output at x=94.8
constexpr double beamY = 142.6; // beam output at y=142.6
constexpr double beamZ = 0.7; // beam output at z=0.7

#endif