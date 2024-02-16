#include "G4ReweightUtils.h"
#include "geant4reweight/src/ReweightBase/G4ReweightStep.hh"

bool protoana::G4ReweightUtils::CreateRWTraj(
    const simb::MCParticle & part, const sim::ParticleList & plist,
    const geo::GeometryCore & geom, int event,
    G4ReweightTraj * theTraj) {

  //Loop over daughters
  for (int i = 0; i < part.NumberDaughters(); ++i) {
    int d_index = part.Daughter(i);
    auto d_part = plist[d_index];
    
    int d_PDG = d_part->PdgCode();
    int d_ID = d_part->TrackId();

    theTraj->AddChild(new G4ReweightTraj(d_ID, d_PDG, part.TrackId(),
                      event, {0,0}));
  }

  //Create process map
  auto procs = part.Trajectory().TrajectoryProcesses();
  std::map<size_t, std::string> proc_map;
  for (auto it = procs.begin(); it != procs.end(); ++it) {
    proc_map[it->first] = part.Trajectory().KeyToProcess(it->second);
  }

  std::vector<double> traj_X, traj_Y, traj_Z;
  std::vector<double> traj_PX, traj_PY, traj_PZ;
  std::vector<size_t> elastic_indices;

  bool found_last = false;
  for (size_t i = 0; i < part.NumberTrajectoryPoints(); ++i) {
    double x = part.Position(i).X();
    double y = part.Position(i).Y();
    double z = part.Position(i).Z();
    
    geo::Point_t test_point{x, y, z};
    const TGeoMaterial * test_material = geom.Material(test_point);

    if (!strcmp(test_material->GetName(), "LAr")) {
      traj_X.push_back(x);
      traj_Y.push_back(y);
      traj_Z.push_back(z);

      traj_PX.push_back(part.Px(i));
      traj_PY.push_back(part.Py(i));
      traj_PZ.push_back(part.Pz(i));

      auto itProc = proc_map.find(i);
      if (itProc != proc_map.end() && itProc->second == "hadElastic") {
        elastic_indices.push_back(i);
      }
    }

    if (i == part.NumberTrajectoryPoints() - 1)
      found_last = true;
  }

  double mass = 0.;

  switch (abs(part.PdgCode())) {
    case 211: {
      mass = 139.57;
      break;
    }
    case 2212: {
      mass = 938.28;
      break;
    }
    default: {
      return false;
      break;
    }
  }

  for (size_t i = 1; i < traj_X.size(); ++i) {
    std::string proc = "default";
    if (found_last && i == traj_X.size() - 1) {
      proc = part.EndProcess();
    }
    else if (std::find(elastic_indices.begin(), elastic_indices.end(), i) !=
             elastic_indices.end()){
      proc = "hadElastic";
    }

    double dX = traj_X[i] - traj_X[i-1];
    double dY = traj_Y[i] - traj_Y[i-1];
    double dZ = traj_Z[i] - traj_Z[i-1];

    double len = sqrt(dX*dX + dY*dY + dZ*dZ);

    double preStepP[3] = {traj_PX[i-1]*1.e3, 
                          traj_PY[i-1]*1.e3, 
                          traj_PZ[i-1]*1.e3};

    double postStepP[3] = {traj_PX[i]*1.e3, 
                           traj_PY[i]*1.e3, 
                           traj_PZ[i]*1.e3};
    if (i == 1) {
      double p_squared = preStepP[0]*preStepP[0] + preStepP[1]*preStepP[1] +
                         preStepP[2]*preStepP[2];
      theTraj->SetEnergy(sqrt(p_squared + mass*mass));
    }

    G4ReweightStep * step = new G4ReweightStep(part.TrackId(), part.PdgCode(),
                                               0, event, preStepP, postStepP,
                                               len, proc);
    theTraj->AddStep(step);
  }

  return true;
}

std::vector<G4ReweightTraj *> protoana::G4ReweightUtils::CreateNRWTrajs(
    const simb::MCParticle & part, const sim::ParticleList & plist,
    const geo::GeometryCore &geom, int event,
    std::string material_name, bool fVerbose) {
  std::vector<G4ReweightTraj *> results;


  //Create process map
  auto procs = part.Trajectory().TrajectoryProcesses();
  std::map<size_t, std::string> proc_map;
  for (auto it = procs.begin(); it != procs.end(); ++it) {
    proc_map[it->first] = part.Trajectory().KeyToProcess(it->second);
  }

  std::vector<double> traj_X, traj_Y, traj_Z;
  //std::vector<double> traj_PX, traj_PY, traj_PZ;
  //std::vector<size_t> elastic_indices;

  std::vector<std::pair<size_t, size_t>> ranges;

  //bool found_last = false;
  bool found_material = false;
  size_t start = 0, end = 0;
  //G4ReweightTraj theTraj(part.TrackId(), part.PdgCode(), 0, event, {0,0});
  if (fVerbose) std::cout << "N traj pts: " <<
                             part.NumberTrajectoryPoints() << std::endl;
  for (size_t i = 0; i < part.NumberTrajectoryPoints(); ++i) {
    double x = part.Position(i).X();
    double y = part.Position(i).Y();
    double z = part.Position(i).Z();

    geo::Point_t test_point{x, y, z};
    const TGeoMaterial * test_material = geom.Material(test_point);
    if (!test_material) continue;
    //if (!strcmp(test_material->GetName(), material_name)) {
    if (test_material->GetName() == material_name) {
      if (fVerbose) {
        std::cout << i << " " << "LAr: " << test_material->GetDensity() << " " <<
                     test_material->GetA() << " " << test_material->GetZ() <<
                     " " << x << " " << y << " " << z << 
                     std::endl;
      }

      if (!found_material) {
        found_material = true;
        start = i;
      }

      //traj_PX.push_back(part.Px(i));
      //traj_PY.push_back(part.Py(i));
      //traj_PZ.push_back(part.Pz(i));

      //auto itProc = proc_map.find(i);
      //if (itProc != proc_map.end() && itProc->second == "hadElastic") {
      //  elastic_indices.push_back(i);
      //}
    }
    else {
      if (fVerbose) {
        std::cout << i << " " << test_material->GetName() << " " <<
                     test_material->GetDensity() << " " <<
                     test_material->GetA() << " " << test_material->GetZ() <<
                     " " << x << " " << y << " " << z << 
                     std::endl;
      }
      if (found_material) {
        found_material = false;
        end = i;
        ranges.push_back({start, end});
      }
    }

    //if (i == part.NumberTrajectoryPoints() - 1)
    //  found_last = true;
  }
  if (found_material) {
    //size_t np = part.NumberTrajectoryPoints();
    ranges.push_back({start, part.NumberTrajectoryPoints() - 1});
    //double x = part.Position(np - 1).X();
    //double y = part.Position(np - 1).Y();
    //double z = part.Position(np - 1).Z();
  }

  double mass = part.Mass()*1.e3;

  /*switch (abs(part.PdgCode())) {
    case 211: {
      mass = 139.57;
      break;
    }
    case 2212: {
      mass = 938.28;
      break;
    }
    case 2112: {
      mass = 939.57;
      break;
    }
    case  321: {

    }
    default: {
      return results;
      break;
    }
  }*/

  for (size_t i = 0; i < ranges.size(); ++i) {
    //std::cout << ranges[i].first << " " << ranges[i].second << std::endl;
    G4ReweightTraj * theTraj = new G4ReweightTraj(i, part.PdgCode(), 0, event, {0,0});
    
    for (size_t j = ranges[i].first; j < ranges[i].second; ++j) {
      double dx = part.Position(j+1).X() - part.Position(j).X();
      double dy = part.Position(j+1).Y() - part.Position(j).Y();
      double dz = part.Position(j+1).Z() - part.Position(j).Z();

      double len = sqrt(dx*dx + dy*dy + dz*dz);
  
      double preStepP[3] = {part.Px(j)*1.e3,
                            part.Py(j)*1.e3,
                            part.Pz(j)*1.e3};
  
      double postStepP[3] = {part.Px(j + 1)*1.e3,
                             part.Py(j + 1)*1.e3,
                             part.Pz(j + 1)*1.e3};
      if (j == ranges[i].first) {
        double p_squared = preStepP[0]*preStepP[0] + preStepP[1]*preStepP[1] +
                           preStepP[2]*preStepP[2];
        theTraj->SetEnergy(sqrt(p_squared + mass*mass));
      }

      std::string proc = "default";
      auto itProc = proc_map.find(j);
      if (itProc != proc_map.end() &&
          j != (part.NumberTrajectoryPoints() - 2)) {
        proc = itProc->second;

        /*if (proc == "Unknown" ) {
          proc = "CoulombScat";
        }*/
      }
      //- 2 because the last element is the end of the last step
      else if (j == (part.NumberTrajectoryPoints() - 2)) {
        proc = part.EndProcess();
      }
      //std::cout << j << " Proc: " << proc << std::endl;
      G4ReweightStep * step = new G4ReweightStep(i, part.PdgCode(),
                                                 0, event, preStepP, postStepP,
                                                 len, proc);
      theTraj->AddStep(step);
    }

    results.push_back(theTraj);
  }

  if (results.size()) {
    //Loop over daughters
    for (int i = 0; i < part.NumberDaughters(); ++i) {
      int d_index = part.Daughter(i);
      auto d_part = plist[d_index];

      int d_PDG = d_part->PdgCode();
      int d_ID = d_part->TrackId();

      if (IsSkippable(d_PDG)) continue;

      results.back()->AddChild(new G4ReweightTraj(d_ID, d_PDG,
                               results.size() - 1, event, {0,0}));

      auto * child = results.back()->GetChildren().back();
      const auto & pos0 = d_part->Position(0);
      const auto & pos1 = d_part->Position(1);
      double d_len = sqrt(
        std::pow((pos0.X() - pos1.X()), 2) +
        std::pow((pos0.Y() - pos1.Y()), 2) +
        std::pow((pos0.Z() - pos1.Z()), 2)
      );
      const auto & p0_lv = d_part->Momentum(0);
      const auto & p1_lv = d_part->Momentum(1);
      double d_p0[3] = {p0_lv.X()*1.e3, p0_lv.Y()*1.e3, p0_lv.Z()*1.e3};
      double d_p1[3] = {p1_lv.X()*1.e3, p1_lv.Y()*1.e3, p1_lv.Z()*1.e3};
      if (abs(d_PDG) == 211 || d_PDG == 111)
        /*std::cout << "Adding step to child " << i << " " << d_PDG << " " <<
                     sqrt(d_p0[0]*d_p0[0] + d_p0[1]*d_p0[1] + d_p0[2]*d_p0[2]) <<
                     " " << d_len << std::endl;*/
      child->AddStep(
          new G4ReweightStep(i, d_PDG, 0, event, d_p0, d_p1, d_len, "default"));
    }
  }
  return results;
}

double protoana::G4ReweightUtils::GetNTrajWeightFromSetPars(
    const std::vector<G4ReweightTraj *> & trajs, G4MultiReweighter & rw) {
  double weight = 1.;
  for (size_t i = 0; i < trajs.size(); ++i) {
    if (trajs[i]->GetNSteps() > 0)
      weight *= rw.GetWeightFromSetParameters(*trajs[i]); 
  }
  return weight;
}

std::pair<double, double> protoana::G4ReweightUtils::GetNTrajPMSigmaWeights(
    const std::vector<G4ReweightTraj *> & trajs, G4MultiReweighter & rw,
    size_t iPar) {
  std::pair<double, double> results = {1., 1.};
  for (size_t i = 0; i < trajs.size(); ++i) {
    if (trajs[i]->GetNSteps() > 0) {
      std::pair<double, double> temp_weight
          = rw.GetPlusMinusSigmaParWeight(*trajs[i], iPar);
      results.first *= temp_weight.first;
      results.second *= temp_weight.second;
    }
  }
  return results;
}

std::vector<std::vector<G4ReweightTraj *>>
    protoana::G4ReweightUtils::BuildHierarchy(
        int ID, int PDG, const sim::ParticleList & plist,
        const geo::GeometryCore & geom, int event,
        std::string material_name, bool skip_first, bool verbose) {

  std::deque<int> to_create = {ID};
  std::vector<std::vector<G4ReweightTraj *>> full_created;

  if (skip_first && verbose) {
    std::cout << "Skipping first" << std::endl;
  }

  while (to_create.size()) {
    auto part = plist[to_create[0]];
    for (int i = 0; i < part->NumberDaughters(); ++i) {
      int daughter_ID = part->Daughter(i);
      auto d_part = plist[daughter_ID];
      if ((d_part->PdgCode() == 2212) || (d_part->PdgCode() == 2112) ||
          (abs(d_part->PdgCode()) == 211) || (abs(d_part->PdgCode()) == 321)) {
        to_create.push_back(daughter_ID);
        if (verbose)
          std::cout << "Adding daughter " << to_create.back() << std::endl;
      }
    }
  
  
    if (skip_first && verbose) {
      std::cout << "Skipping " << ID << std::endl;
    }

    if (!skip_first || (skip_first && ID != part->TrackId())) {
      if (skip_first && verbose) {
        std::cout << "Not skipping " << part->TrackId() << std::endl;
      }
      std::vector<G4ReweightTraj *> temp_trajs =
          CreateNRWTrajs(*part, plist, geom,
                         event, material_name, verbose);
    
      if (temp_trajs.size()) {
        auto last_traj = temp_trajs.back();
        if (verbose) 
          std::cout << "created " << last_traj->GetTrackID() << " " <<
                       last_traj->GetPDG() << std::endl;
  
        if (temp_trajs[0]->GetPDG() == PDG) {
          full_created.push_back(temp_trajs);
        }
      }
    }
    to_create.pop_front();
  }

  return full_created;
}

bool protoana::G4ReweightUtils::IsSkippable(int pdg) {
  return ((abs(pdg) == 11) || (pdg == 22) || (pdg > 1000000000));
}
