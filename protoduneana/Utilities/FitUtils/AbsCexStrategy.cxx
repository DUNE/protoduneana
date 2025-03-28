#include "AbsCexStrategy.h"
#include "ThinSliceStrategyFactory.h"
#include "cetlib_except/exception.h"

#include "ThinSliceConsts.h"
DECLARE_THINSLICESTRATEGY_FACTORY_NS(protoana::AbsCexStrategy, protoana, AbsCexStrategy)



protoana::AbsCexStrategy::AbsCexStrategy(const fhicl::ParameterSet & pset)
    : ThinSliceStrategy(pset),
      fERecoSelections(pset.get<std::vector<int>>("ERecoSelections", {})),
      fEndZSelections(pset.get<std::vector<int>>("EndZSelections", {})),
      fOneBinSelections(pset.get<std::vector<int>>("OneBinSelections", {})),
      fPitch(pset.get<double>("WirePitch")),
      fSliceCut(pset.get<int>("SliceCut")),
      fTrajZStart(pset.get<double>("TrajZStart"))
    {

}


void protoana::AbsCexStrategy::BuildDists(
    ThinSliceDistHolder & holder, 
    const fhicl::ParameterSet & pset,
    std::string label) {
  builder.Build(holder, pset, label);
}

double protoana::AbsCexStrategy::GetEventWeight(
    const ThinSliceEvent & event) const {
  // Implement the logic to calculate the event weight
  // based on the event data and any other necessary parameters
  return 1.0; // Placeholder value, replace with actual calculation
}

int protoana::AbsCexStrategy::GetBeamBin(
    const std::vector<double> & beam_energy_bins,
    const ThinSliceEvent & event,
    bool restrict_P) const {
  
        double momentum = event.GetBeamInstP();

        int bin = -1;
        for (size_t j = 1; j < beam_energy_bins.size(); ++j) {
        if ((beam_energy_bins[j-1] <= 1.e3*momentum) &&
            (1.e3*momentum < beam_energy_bins[j])) {
            bin = j - 1;
            break;
        }
        }
        if (bin == -1 && !restrict_P) {
            std::string message = "Could not find beam energy bin for " +
                                    std::to_string(momentum);
            throw std::runtime_error(message);
        }
        return bin;
  }

void protoana::AbsCexStrategy::FillHistsFromEvent(
    const ThinSliceEvent & event, ThinSliceDistHolder & dists,
    size_t beam_bin,
    std::mutex & mutex,
    double weight) const {
        double val = 0.;
    
        int new_selection = event.GetSelectionID();
    
        // //Change the selection here
        // if (!fFakeDataActive) {
        //   if (fVaryMCCalibration && fVaryMCCalSelection) { 
        //     new_selection = (
        //       (fMCCalibrationFactor > 1.) ?
        //       event.GetCalUpSelectionID() :
        //       event.GetCalDownSelectionID()
        //     );
        //   }
        //   else if (fMCBackUpSelection) {
        //     new_selection = event.GetBackUpSelectionID(); 
        //   }
        //   else if (fMCBackDownSelection) {
        //     new_selection = event.GetBackDownSelectionID(); 
        //   }
        //   else if (fMCFrontUpSelection) {
        //     new_selection = event.GetFrontUpSelectionID(); 
        //   }
        //   else if (fMCFrontDownSelection) {
        //     new_selection = event.GetFrontDownSelectionID(); 
        //   }
        // }
    
        // std::cout << event.GetSampleID() << " " << new_selection << " " <<
        //              beam_bin << " " << GetSignalBin(event, dists) << std::endl;
        auto selected_hist
            = dists.GetSelectionHist(
                event.GetSampleID(), new_selection,
                beam_bin, GetSignalBin(event, dists)
            );

        double reco_beam_endZ = event.GetRecoEndZ();

        const std::vector<double> & reco_beam_incidentEnergies
        = event.GetRecoIncidentEnergies();
        double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();

        if (std::find(fEndZSelections.begin(), fEndZSelections.end(), new_selection)
            != fEndZSelections.end()) {
          if (selected_hist->FindBin(reco_beam_endZ) == 0) {
            val = selected_hist->GetBinCenter(1);
          }
          else if (selected_hist->FindBin(reco_beam_endZ) >
                   selected_hist->GetNbinsX()) {
            val = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
          }
          else {
            val = reco_beam_endZ;
          }
        }
        else if (
            std::find(fOneBinSelections.begin(), fOneBinSelections.end(), new_selection)
            != fOneBinSelections.end()) {
          val = .5;
        }
        else if (reco_beam_incidentEnergies.size()) {
    
          bool fDoEnergyFix = true;
          double fEnergyFix = 10.;

          double energy = 0.;
          double beam_energy_delta = 0.;
        //   if (fFakeResolution < 0) {
            energy = {reco_beam_interactingEnergy + beam_energy_delta};
            if (fDoEnergyFix) {
              for (size_t k = 1; k < reco_beam_incidentEnergies.size(); ++k) {
                double deltaE = ((reco_beam_incidentEnergies)[k-1] -
                                 (reco_beam_incidentEnergies)[k]);
                // if (fVaryCalibrationFakeData && fFillFakeDataInMain) {
                //   energy += deltaE;
                //   if (deltaE*fDataCalibrationFactor < fEnergyFix) {
                //     energy -= deltaE*fDataCalibrationFactor;
                //   }
                // }
                // else if (fVaryMCCalibration && !fFakeDataActive) {
                //   energy += deltaE;
                //   if (deltaE*fMCCalibrationFactor < fEnergyFix) {
                //     energy -= deltaE*fMCCalibrationFactor;
                //   }
                // }
                // else {
                  if (deltaE > fEnergyFix) {
                    energy += deltaE; 
                  }
                // }
              }
            }
        //   }
        //   else {
        //     if (fUseStoredEnergies) {
        //       energy[0] = fStoredEnergies[worker_id][energy_index]; 
        //       ++energy_index;
        //     }
        //     else {
        //       energy[0] = fRNG.Gaus(end_energy, fFakeResolution);
        //       fStoredEnergies[worker_id].push_back(energy[0]);
        //     }
        //   }
    
          if (selected_hist->FindBin(energy) == 0) {
              val = selected_hist->GetBinCenter(1);
          }
          else if (selected_hist->FindBin(energy) > selected_hist->GetNbinsX()) {
            val = selected_hist->GetBinCenter(selected_hist->GetNbinsX());
          }
          else {
            val = energy;
          }
        }
        else {
          val = selected_hist->GetBinCenter(1);
        }

        std::shared_ptr<TH1> interaction_hist = nullptr;
        bool fill_interaction = dists.IsInteraction(event.GetSampleID());
        if (fill_interaction) {
            interaction_hist = dists.GetInteractionHist(
                event.GetSampleID(),
                beam_bin
            );
        }

        std::vector<double> inc_energies;
        bool fill_incident = dists.IsIncident(event.GetSampleID());
        if (fill_incident) inc_energies = MakeTrueIncidentEnergies(event.GetTrueTrajZ(), event.GetTrueTrajKE());

        std::lock_guard<std::mutex> lock(mutex);
        selected_hist->Fill(val, weight);
        if (fill_interaction) {
            interaction_hist->Fill(GetTrueEndEnergy(event), weight);
        }
        if (fill_incident) {
            //Loop over the measurement incident histograms and fill with 
            //any contributions to it from this event 
            for (const auto & map : dists.GetIncidentHists()) {
                for (auto & [key, hist] : map) {
                    for (const auto & inc_energy : inc_energies) {
                        hist->Fill(inc_energy, weight);
                    }        
                }
            }
        }

}


//Put in a util class?
std::vector<double> protoana::AbsCexStrategy::MakeTrueIncidentEnergies(
    const std::vector<double> & true_beam_traj_Z,
    const std::vector<double> & true_beam_traj_KE) const {

    std::vector<double> results;
    double next_slice_z = fTrajZStart;
    int next_slice_num = 0;
    for (size_t j = 1; j < true_beam_traj_Z.size() - 1; ++j) {
        double z = true_beam_traj_Z[j];
        double ke = true_beam_traj_KE[j];

        if (z < fTrajZStart) {
            continue;
        }

        if (z >= next_slice_z) {
            double temp_z = true_beam_traj_Z[j-1];
            double temp_e = true_beam_traj_KE[j-1];
            
            while (next_slice_z < z && next_slice_num < fSliceCut) {
                double sub_z = next_slice_z - temp_z;
                double delta_e = true_beam_traj_KE[j-1] - ke;
                double delta_z = z - true_beam_traj_Z[j-1];
                temp_e -= (sub_z/delta_z)*delta_e;
                results.push_back(temp_e);
                temp_z = next_slice_z;
                next_slice_z += fPitch;
                ++next_slice_num;
            }
        }
    }

    return results;
}

//PUt in a util class
double protoana::AbsCexStrategy::GetTrueEndEnergy(const ThinSliceEvent & event) const {
    double true_beam_endP = event.GetTrueEndP();
    return sqrt(true_beam_endP*true_beam_endP*1.e6 + ChargedPiMass*ChargedPiMass) - ChargedPiMass;
}

int protoana::AbsCexStrategy::GetSignalBin(const ThinSliceEvent & event, const ThinSliceDistHolder & dists) const {
    int sample_ID = event.GetSampleID();

    const auto & signal_ids = dists.GetSignalIDs();

    //Non signal events
    if (std::find(signal_ids.begin(), signal_ids.end(), sample_ID) == signal_ids.end()) {
        return 0;
    }

    const auto & signal_ranges = dists.GetSignalRanges();

    //Single bin signal events
    if (signal_ranges.find(sample_ID) == signal_ranges.end()) {
        return 0;
    }

    //Multi bin signal events
    const auto & ranges = signal_ranges.at(sample_ID);
    double end_energy = GetTrueEndEnergy(event); //sqrt(true_beam_endP*true_beam_endP*1.e6 + ChargedPiMass*ChargedPiMass) - ChargedPiMass;

    //Underflow
    if (end_energy < ranges[0].first) {
        return 0;
    }
    //Overflow
    if (end_energy > ranges.back().second) {
        return ranges.size()+1;
    }
    for (size_t i = 0; i < ranges.size(); ++i) {
        if (end_energy > ranges[i].first && end_energy < ranges[i].second) {
            return i+1;
        }
    }

    //If we reached here then something went wrong
    throw cet::exception("AbsCexStrategy")
        << "Could not find signal bin for event with end energy: "
        << end_energy << " and sample ID: " << sample_ID << "\n";
}


void protoana::AbsCexStrategy::CompareDataMC(const ThinSliceDistHolder & holder, const ThinSliceDataSet & dataset, TFile & fout) const {
    fout.cd();
    for (const auto & map : holder.GetInteractionHists()) {
        for (const auto & [true_id, hist] : map) {
            hist->Write();
        }
    }
    for (const auto & map : holder.GetIncidentHists()) {
        for (const auto & [true_id, hist] : map) {
            hist->Write();
        }
    }
    for (const auto & [true_id, hist] : holder.GetXSecHists()) {
        hist->Write();
        holder.GetTotalIncidentHists().at(true_id)->Write();
    }
    for (const auto & map : holder.GetSelectionHists()) {
        for (const auto & [ids, hists] : map) {
            for (const auto & hist : hists) {
                hist->Write();
            }
        }
    }
}

void protoana::AbsCexStrategy::CalcXSecs(ThinSliceDistHolder & holder, double scale) const {
    builder.CalcXSecs(holder, scale);
}
