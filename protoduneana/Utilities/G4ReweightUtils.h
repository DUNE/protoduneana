#ifndef G4REWEIGHTUTILS_h
#define G4REWEIGHTUTILS_h

#include <vector>

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "geant4reweight/src/ReweightBase/G4ReweightTraj.hh"
#include "geant4reweight/src/ReweightBase/G4MultiReweighter.hh"

namespace protoana {
 namespace G4ReweightUtils {

  bool IsSkippable(int pdg);

  //Function to create reweightable object
  bool CreateRWTraj(const simb::MCParticle & part,
                    const sim::ParticleList & plist,
                    const geo::GeometryCore & geom, int event,
                    G4ReweightTraj * theTraj); 

  std::vector<G4ReweightTraj *> CreateNRWTrajs(
      const simb::MCParticle & part,
      const sim::ParticleList & plist,
      const geo::GeometryCore & geom, int event,
      std::string material_name  = "LAr",
      bool fVerbose=false);

  std::vector<std::vector<G4ReweightTraj *>> BuildHierarchy(
      int ID, int PDG,
      const sim::ParticleList & plist,
      const geo::GeometryCore & geom, int event,
      std::string material_name = "LAr", bool skip_first=false,
      bool fVerbose=false);

  double GetNTrajWeightFromSetPars(
    const std::vector<G4ReweightTraj *> & trajs, G4MultiReweighter & rw);

  std::pair<double, double> GetNTrajPMSigmaWeights(
    const std::vector<G4ReweightTraj *> & trajs, G4MultiReweighter & rw,
    size_t iPar);
 }
}

#endif
