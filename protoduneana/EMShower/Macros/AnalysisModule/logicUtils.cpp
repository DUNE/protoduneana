#include "logicUtils.h"

namespace logic_utils {
    bool particleInCone(double particleX, double particleY, double particleZ, double height, double r1, double r2) {
        // Cone parameters
        const double angleX0Z = -10.825 * M_PI / 180.0; // 169.175-180
        const double angleY0Z = -44.486 * M_PI / 180.0; // 135.514-180

        // Unit vector along cone axis
        double axisX = sin(angleX0Z);
        double axisY = sin(angleY0Z);
        double axisZ = cos(angleX0Z)*cos(angleY0Z);

        // Normalize the direction vector
        double norm = sqrt(axisX*axisX + axisY*axisY + axisZ*axisZ);
        axisX /= norm;
        axisY /= norm;
        axisZ /= norm;

        // Cone apex coordinates
        const double coneX = beamX-30*axisX;
        const double coneY = beamY-30*axisY;
        const double coneZ = beamZ-30*axisZ;

        // Vector from cone apex to point
        double dx = particleX-coneX;
        double dy = particleY-coneY;
        double dz = particleZ-coneZ;

        // Project point onto cone axis
        double projection = dx*axisX + dy*axisY + dz*axisZ;

        // Check if projection is within the cone height
        if (projection<0 || projection>height) return false;

        // Calculate allowed radius at this height
        double fraction = projection/height;
        double allowedRadius = r1+(r2-r1)*fraction;

        // Calculate perpendicular distance from point to axis
        double perpX = dx-projection*axisX;
        double perpY = dy-projection*axisY;
        double perpZ = dz-projection*axisZ;
        double perpDistance = sqrt(perpX*perpX + perpY*perpY + perpZ*perpZ);

        return perpDistance <= allowedRadius;
    }

    std::vector<bool> TrackBeforeShower(unsigned int fNParticles, double* fShowerStartX, double* fShowerStartY, double* fShowerStartZ,
        double* fTrackStartX, double* fTrackStartY, double* fTrackStartZ) {
        std::vector<bool> isTrackBeforeShower(fNParticles, false);
        for (unsigned int i=0; i<fNParticles; ++i) {
            if (fShowerStartX[i]==-999 && fShowerStartY[i]==-999 && fShowerStartZ[i]==-999) continue;
            if (fTrackStartX[i]==0 && fTrackStartY[i]==0 && fTrackStartZ[i]==0) continue;
            for (unsigned int j=0; j<fNParticles; j++) {
                if (fTrackStartZ[i]<fShowerStartZ[j]) {
                    isTrackBeforeShower[i] = true;
                    continue;
                }
            }        
        }
        return isTrackBeforeShower;
    }
}

