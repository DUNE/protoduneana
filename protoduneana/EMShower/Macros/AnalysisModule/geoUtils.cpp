#include "geoUtils.h"

namespace geo_utils {
    void createAxis(double length, double dx, double dy, double dz, double x, double y, double z,
        Color_t color, const char* name) {
        TEveArrow* Axis = new TEveArrow(length*dx, length*dy, length*dz, x, y, z);
        Axis->SetMainColor(color);
        Axis->SetTubeR(0.002);
        Axis->SetPickable(kTRUE);
        Axis->SetName(name);
        gEve->AddElement(Axis);
    }

    void createCone(double refX, double refY, double refZ,
        double angle_x0z, double angle_y0z,
        double height, double r1, double r2) {
        // Créer la géométrie
        TGeoManager *geom = new TGeoManager("world", "World with cone");
        TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
        TGeoMedium *vacuum = new TGeoMedium("Vacuum", 1, matVacuum);

        TGeoVolume *top = geom->MakeBox("Top", vacuum, 500, 500, 500);
        geom->SetTopVolume(top);

        // Définir le cône
        double rmin1 = 0;
        double rmin2 = 0;
        TGeoCone *solidCone = new TGeoCone(height, rmin1, r1, rmin2, r2);
        TGeoVolume *coneSolid = new TGeoVolume("coneSolid", solidCone, vacuum);
        coneSolid->SetLineColor(kWhite);
        coneSolid->SetTransparency(50);            

        // Direction en radians
        double theta_xz = (angle_x0z-180) * TMath::DegToRad();
        double theta_yz = (angle_y0z-180) * TMath::DegToRad();

        // Calcul direction normalisée
        double dx = sin(theta_xz);
        double dy = sin(theta_yz);
        double dz = cos(theta_xz) * cos(theta_yz);

        // Normaliser
        double norm = sqrt(dx*dx + dy*dy + dz*dz);
        dx /= norm;
        dy /= norm;
        dz /= norm;

        // Axes selon la direction du cône
        // geo_utils::createAxis(1000.0, dx, dy, dz, refX, refY, refZ, kYellow, "Cone direction");

        // Create rotation first
        TGeoRotation *rot = new TGeoRotation();
        rot->RotateX(angle_y0z);
        rot->RotateY(angle_x0z);

        // Calculate the offset needed to place the summit at the reference point
        // The cone's summit is at height/2 in local coordinates
        double localSummitZ = height-30;

        // Calculate the translation that puts the summit at the reference point
        // We need to apply the inverse rotation to our desired direction vector
        double translation_x = refX + (localSummitZ * dx);
        double translation_y = refY + (localSummitZ * dy);
        double translation_z = refZ + (localSummitZ * dz);

        // Create the translation
        TGeoTranslation *trans = new TGeoTranslation(translation_x, translation_y, translation_z);

        // Combine rotation and translation
        TGeoCombiTrans *combi = new TGeoCombiTrans(*trans, *rot);

        // Ajouter le cône avec transformation
        top->AddNode(coneSolid, 1, combi);

        geom->CloseGeometry();
    }

    TEvePointSet* createPointSet(const char* name, int markerStyle, double markerSize, Color_t color) {
        TEvePointSet *refPoint = new TEvePointSet(name);
        refPoint->SetMarkerStyle(markerStyle); // style point
        refPoint->SetMarkerSize(markerSize);
        refPoint->SetMarkerColor(color);
        return refPoint;
    }

    void createLine(double startX, double startY, double startZ, double endX, double endY, double endZ,
        int lineWidth, Color_t color) {
        // Draw a line from track start to track end
        TEveLine* trackLine = new TEveLine();
        trackLine->SetLineWidth(lineWidth);
        trackLine->SetMainColor(color); // Choose your color
        trackLine->SetNextPoint(startX, startY, startZ);
        trackLine->SetNextPoint(endX, endY, endZ);
        gEve->AddElement(trackLine);
    }
}