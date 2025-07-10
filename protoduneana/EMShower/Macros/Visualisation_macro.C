
#include "AnalysisModule/SpecialMacroClass.h"
#include "AnalysisModule/geoUtils.h"

void Visualisation_macro() {
    SpecialMacro tree("/exp/dune/app/users/chalamet/development/srcs/duneana/duneana/AnaTest/electron_1GeV_BGCosmics/1/analysisOutput.root", 1);

    // === Initialiser TEve ===
    TEveManager::Create();

    // === Créer les axes X, Y, Z ===
    double axisLength = 1000.0;
    geo_utils::createAxis(axisLength, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, kRed, "X axis");
    geo_utils::createAxis(axisLength, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, kGreen, "Y axis");
    geo_utils::createAxis(axisLength, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, kBlue, "Z axis");
    
    
    // === Créer un point de référence ===
    float refX = 94.8; // beam output at x=94.8
    float refY = 142.6; // beam output at y=142.6
    float refZ = 0.7; // beam output at z=0.7
    TEvePointSet *refPoint = geo_utils::createPointSet("Reference Point", kCircle, 2.0, kRed);
    refPoint->SetNextPoint(refX, refY, refZ); // Position du point de référence
    
    // === Paramètres du cône ===
    double height = 200.0;
    double r1 = 30.0;
    double r2 = 60.0;
    double angle_x0z = 169.1;
    double angle_y0z = 135.514;
    
    // === Créer la géométrie avec un cône ===
    geo_utils::createCone(beamX, beamY, beamZ, angle_x0z, angle_y0z, height/2, r1, r2);
    
    // === Populate window ===
    TEvePointSet *TrackStart = geo_utils::createPointSet("track start", kCircle, 1.0, kCyan);
    TEvePointSet *Shower = geo_utils::createPointSet("shower", kCircle, 1.0, kYellow);
    for (Long64_t i=0; i<tree.nentries; ++i) {
        if (i != 1) continue; // Skip to the desired event
        tree.tree->GetEntry(i);
        for (unsigned int j=0; j<tree.fNParticles; j++) {
            TrackStart->SetNextPoint(tree.fTrackStartX[j], tree.fTrackStartY[j], tree.fTrackStartZ[j]);
            Shower->SetNextPoint(tree.fShowerStartX[j], tree.fShowerStartY[j], tree.fShowerStartZ[j]);
            if (tree.fTrackLength[j]>50) {
                geo_utils::createLine(tree.fTrackStartX[j], tree.fTrackStartY[j], tree.fTrackStartZ[j],
                    tree.fTrackEndX[j], tree.fTrackEndY[j], tree.fTrackEndZ[j], 2, kCyan);
            }
        }
    }
    
    // === Ajouter les éléments au viewer TEve ===
    gEve->AddElement(TrackStart);
    gEve->AddElement(Shower);
    gEve->AddElement(refPoint);
    
    // === Lier la géométrie au viewer TEve ===
    if (gGeoManager) {
        TEveGeoTopNode *topNode = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
        gEve->AddGlobalElement(topNode);
    }

    // === Affichage ===
    gEve->Redraw3D(kTRUE);
}
