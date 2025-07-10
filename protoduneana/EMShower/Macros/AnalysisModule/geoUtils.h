#ifndef GEO_UTILS_H
#define GEO_UTILS_H

#include <TEveManager.h>
#include <TEvePointSet.h>
#include <TEveArrow.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoCone.h>
#include <TGeoMatrix.h>
#include <TEveLine.h>

namespace geo_utils {
    void createAxis(double length, double dx, double dy, double dz, double x, double y, double z, Color_t color,
        const char* name);

    void createCone(double xRef, double yRef, double zRef, double angle_x0z, double angle_y0z,
        double height, double rmax1, double rmax2);

    TEvePointSet* createPointSet(const char* name, int markerStyle, double markerSize, Color_t color);

    void createLine(double startX, double startY, double startZ, double endX, double endY, double endZ,
    int lineWidth, Color_t color);
}

#endif