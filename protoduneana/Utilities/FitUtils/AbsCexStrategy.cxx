#include "AbsCexStrategy.h"
#include "ThinSliceStrategyFactory.h"
#include "cetlib_except/exception.h"

#include "ThinSliceConsts.h"

#include "THStack.h"
#include "TLegend.h"
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
  fBuilder.Build(holder, pset, label);
}

double protoana::AbsCexStrategy::GetEventWeight(
    const ThinSliceEvent & event) const {
  // Implement the logic to calculate the event weight
  // based on the event data and any other necessary parameters
  return 1.0; // Placeholder value, replace with actual calculation
}

void protoana::AbsCexStrategy::FillDataHistsFromEvent(
    ThinSliceDataSet & data_set, const ThinSliceEvent & event) const {


  int selection_ID = event.GetSelectionID();
  auto & selected_hists = data_set.GetSelectionHists();
  auto * selected_hist = selected_hists[selection_ID];

  double reco_beam_endZ = event.GetRecoEndZ();
  const std::vector<double> & reco_energies
      = event.GetRecoIncidentEnergies();
  double reco_beam_interactingEnergy = event.GetRecoInteractingEnergy();

  double val = 0.;

  bool fDoEnergyFix = true;
  double fEnergyFix = 10.;

  if (std::find(fEndZSelections.begin(), fEndZSelections.end(), selection_ID)
      != fEndZSelections.end()) {
    if (selected_hist->FindBin(reco_beam_endZ) == 0) {
      val = selected_hist->GetBinCenter(1);
    }
    else if (selected_hist->FindBin(reco_beam_endZ) > selected_hist->GetNbinsX()) {
      val = selected_hist->GetBinCenter(
          selected_hist->GetNbinsX());
    }
    else {
      val = reco_beam_endZ;
    }
  }
  else if (
    std::find(fOneBinSelections.begin(), fOneBinSelections.end(),
              selection_ID)
    != fOneBinSelections.end()) {
    val = .5;
  }
  else if (reco_energies.size()) {
    if (selected_hists.find(selection_ID) != selected_hists.end()) {
      double energy = reco_beam_interactingEnergy/* + deltaE_scale*/;

      if (fDoEnergyFix) {
        for (size_t k = 1; k < reco_energies.size(); ++k) {
          double deltaE = (reco_energies[k-1] - reco_energies[k]);
          // if (fVaryDataCalibration) {
          //   energy += deltaE;
          //   if (deltaE*fDataCalibrationFactor < fEnergyFix) {
          //     energy -= deltaE*fDataCalibrationFactor;
          //   }
          // }
          // else {
            if (deltaE > fEnergyFix) {
              energy += deltaE; 
            }
          // }
        }
      }
      if (selected_hist->FindBin(energy) == 0) {
        val = selected_hist->GetBinCenter(1);
      }
      else if (selected_hist->FindBin(energy) >
                selected_hist->GetNbinsX()) {
        val = selected_hist->GetBinCenter(
            selected_hist->GetNbinsX());
      }
      else {
        val = energy;
      }
    }
  }
  else {
    val = selected_hist->GetBinCenter(1);
  }
  selected_hist->Fill(val);
}

//TODO -- remove
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

        const std::vector<double> & reco_energies
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
        else if (reco_energies.size()) {
    
          bool fDoEnergyFix = true;
          double fEnergyFix = 10.;

          double energy = 0.;
          double beam_energy_delta = 0.;
        //   if (fFakeResolution < 0) {
            energy = {reco_beam_interactingEnergy + beam_energy_delta};
            if (fDoEnergyFix) {
              for (size_t k = 1; k < reco_energies.size(); ++k) {
                double deltaE = ((reco_energies)[k-1] -
                                 (reco_energies)[k]);
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
        bool fill_incident = (fFillIncident && dists.IsIncident(event.GetSampleID()));
        if (fill_incident) inc_energies = MakeTrueIncidentEnergies(event.GetTrueTrajZ(), event.GetTrueTrajKE());

        std::lock_guard<std::mutex> lock(mutex);
        selected_hist->Fill(val, weight);
        if (fill_interaction) {
            interaction_hist->Fill(GetTrueEndEnergy(event), weight);
        }
        if (fill_incident) {
            //Loop over the measurement incident histograms and fill with 
            //any contributions to it from this event 
            // for (const auto & map : dists.GetIncidentHists()) {
                auto & map = dists.GetIncidentHists()[beam_bin];
                for (auto & [key, hist] : map) {
                    for (const auto & inc_energy : inc_energies) {
                        hist->Fill(inc_energy, weight);
                    }
                }
            // }
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

void protoana::AbsCexStrategy::CompareSelections(
    const ThinSliceDistHolder & holder,
    ThinSliceDataSet & data_set, TFile & output_file,
    std::vector<std::pair<int, int>> plot_style,
    bool plot_rebinned,
    bool post_fit,
    TDirectory * plot_dir) const {

  plot_dir->cd();
  //output_file.cd();
  std::map<int, TH1*> data_hists
      = (plot_rebinned ?
         data_set.GetRebinnedSelectionHists() :
         data_set.GetSelectionHists());
  for (const auto & [selection_ID, data_hist] : data_hists) {  
    // int selection_ID = it->first;
    // TH1D * data_hist = (TH1D*)it->second;
    data_hist->SetLineColor(kBlack);
    data_hist->SetMarkerColor(kBlack);
    data_hist->SetMarkerStyle(20);

    std::string canvas_name_no_data = "c";
    canvas_name_no_data += (post_fit ? "PostFit" : "Nominal") +
                   data_set.GetSelectionName(selection_ID) +
                   "NoData";
    TCanvas cSelectionNoData(canvas_name_no_data.c_str(), "");
    cSelectionNoData.SetTicks();

    std::string canvas_name = "c";
    canvas_name += (post_fit ? "PostFit" : "Nominal") +
                   data_set.GetSelectionName(selection_ID);
    TCanvas cSelection(canvas_name.c_str(), "");
    cSelection.SetTicks();

    const auto & mc_selection_hists = holder.GetSelectionHists();
    std::map<int, std::vector<TH1D *>> temp_hists;
    //Loop over beam bins
    for (size_t i = 0; i < mc_selection_hists.size(); ++i) {
      const auto & mc_selection_map = mc_selection_hists[i];
      
      //Loop over true, selection ID pairs
      for (const auto & [ids, hist_vec] : mc_selection_map) {
        int sample_ID = ids.first;
        int mc_selection_ID = ids.second;
        if (mc_selection_ID != selection_ID) continue;
        if (i == 0) {
          temp_hists[sample_ID] = std::vector<TH1D *>();
          for (const auto & hist : hist_vec) {
            temp_hists[sample_ID].push_back(
              (TH1D*)hist->Clone()
            );
            // temp_hists[sample_ID].push_back((TH1D*)(
            // plot_rebinned ?
            // hist.GetRebinnedSelectionHist(selection_ID) :
            // hist.GetSelectionHist(selection_ID))->Clone());
          }
        }
        else {
          for (size_t j = 0; j < hist_vec.size(); ++j) {
            temp_hists[sample_ID][j]->Add(
              hist_vec[j].get()
            );
            // temp_hists[sample_ID].push_back((TH1D*)(
            // plot_rebinned ?
            // hist.GetRebinnedSelectionHist(selection_ID) :
            // hist.GetSelectionHist(selection_ID))->Clone());
          }
        }
      }
    }
    //Build the mc stack
    std::string stack_name = (post_fit ? "PostFit" : "Nominal") +
                             data_set.GetSelectionName(selection_ID) +
                             "Stack";
    THStack mc_stack(stack_name.c_str(), "");
    TLegend leg;
    std::vector<TH1D*> temp_vec;
    size_t iColor = 0;
    //need to add second loop with temp hists
    for (auto it2 = temp_hists.begin(); it2 != temp_hists.end(); ++it2) {
      for (size_t i = 0; i < it2->second.size(); ++i) {
        TH1D * sel_hist = it2->second.at(i);
        std::pair<int, int> color_fill = GetColorAndStyle(iColor, plot_style);
        sel_hist->SetFillColor(color_fill.first);
        sel_hist->SetFillStyle(color_fill.second);
        sel_hist->SetLineColor(kBlack);
        mc_stack.Add(sel_hist);
        temp_vec.push_back(sel_hist);
        ++iColor;
      }
    }
    mc_stack.Write();

    for (auto it = temp_vec.rbegin(); it != temp_vec.rend(); ++it) {
      leg.AddEntry(*it);
    }

    cSelectionNoData.cd();
    mc_stack.Draw("hist");
    std::string title_no_data = ";";
    title_no_data += data_hist->GetXaxis()->GetTitle();
    mc_stack.SetTitle(title_no_data.c_str());
    mc_stack.GetHistogram()->SetTitleSize(.04, "X");
    // mc_stack.Draw("hist");
    // if (it == data_hists.begin())
    //   leg.Write("leg_no_data");
    cSelectionNoData.Write();
    leg.Draw("same");

    leg.AddEntry(data_hist, "Data");
    
    double chi2 = CalcChi2(holder, data_set);
    if (chi2 < 0. && chi2 > -1.e7) {
      chi2 = 0.;
    }
    TString chi2_str;
    chi2_str.Form("#chi^{2} = %.2f", chi2);
    leg.AddEntry((TObject*)0x0, chi2_str, "");

    cSelection.cd();
    std::string title = data_set.GetSelectionName(selection_ID);
    title += ";";
    title += data_hist->GetXaxis()->GetTitle();
    mc_stack.SetTitle(title.c_str());

    mc_stack.Draw("hist");
    double max_mc = mc_stack.GetHistogram()->GetMaximum();
    int max_data_bin = data_hist->GetMaximumBin();
    double max_data = data_hist->GetBinContent(max_data_bin) +
                      data_hist->GetBinError(max_data_bin);
    mc_stack.SetMaximum((max_data > max_mc ? max_data : max_mc));
    mc_stack.Draw("hist");
    data_hist->Draw("e1 same");
    leg.Draw("same");

    cSelection.Write();

    //Get the full incident hist from stack
    TList * l = (TList*)mc_stack.GetHists();
    TH1D * hMC = (TH1D*)l->At(0)->Clone();
    for (int i = 1; i < l->GetSize(); ++i) {
      hMC->Add((TH1D*)l->At(i));
    }

    std::string ratio_name = data_set.GetSelectionName(selection_ID) + "Ratio" +
                             (post_fit ? "PostFit" : "Nominal");
    TH1D * hRatio
        = (TH1D*)data_hist->Clone(ratio_name.c_str());
    hRatio->Divide(hMC);
    hRatio->Write(); 
    std::string total_name = (post_fit ? "PostFit" : "Nominal") +
                             data_set.GetSelectionName(selection_ID) +
                             "Total";
    hMC->Write(total_name.c_str());

    canvas_name += "Ratio";
    TCanvas cRatio(canvas_name.c_str(), "");
    cRatio.SetTicks();
    TPad p1((canvas_name + "pad1").c_str(), "", 0, 0.3, 1., 1.);
    p1.SetBottomMargin(0.1);
    p1.Draw();
    p1.cd();
    mc_stack.Draw("hist");
    mc_stack.GetHistogram()->SetTitle("Abs;;");
    for (int i = 1; i < mc_stack.GetHistogram()->GetNbinsX(); ++i) {
      hRatio->GetXaxis()->SetBinLabel(
          i, mc_stack.GetHistogram()->GetXaxis()->GetBinLabel(i));
      mc_stack.GetHistogram()->GetXaxis()->SetBinLabel(i, "");
    }
    mc_stack.Draw("hist");
    data_hist->Draw("e1 same");

    cRatio.cd();
    TPad p2((canvas_name + "pad2").c_str(), "", 0, 0, 1., 0.3);
    p2.SetTopMargin(0.1);
    p2.SetBottomMargin(.2);
    p2.Draw();
    p2.cd();
    hRatio->Sumw2();
    std::string r_title = "";
    hRatio->GetYaxis()->SetTitle("Data/MC");
    hRatio->SetTitleSize(.04, "XY");
    hRatio->Draw("ep");

    cRatio.Write();


    //Differences
    //Get the full incident hist from stack
    // TH1D * hMC2 = (TH1D*)l->At(0)->Clone();
    // for (int i = 1; i < l->GetSize(); ++i) {
    //   hMC2->Add((TH1D*)l->At(i));
    // }

    // std::string diff_name = data_set.GetSelectionName(selection_ID) + "Diff" +
    //                          (post_fit ? "PostFit" : "Nominal");
    // TH1D * hDiff
    //     = (TH1D*)data_hist->Clone(diff_name.c_str());
    // hMC2->Scale(-1.);
    // hDiff->Add(hMC2);
    // hMC2->Scale(-1.);
    // hDiff->Divide(hMC2);
    // hDiff->Write(); 

    // canvas_name += "Diff";
    // TCanvas cDiff(canvas_name.c_str(), "");
    // cDiff.SetTicks();
    // cDiff.cd();
    // hDiff->GetYaxis()->SetTitle("r");
    // hDiff->SetTitleSize(.04, "XY");
    // hDiff->Draw("ep");

    // cDiff.Write();

    ///Write in NoStacks here

    // THStack full_mc_stack((stack_name + "Full").c_str(), "");
    // //TLegend leg;
    // size_t iColorFull = 0;
    // //need to add second loop with temp hists
    // for (auto it2 = temp_hists.begin(); it2 != temp_hists.end(); ++it2) {
    //   TH1D * sel_hist = it2->second.at(0);
    //   std::pair<int, int> color_fill = GetColorAndStyle(iColorFull, plot_style);
    //   sel_hist->SetFillColor(color_fill.first);
    //   sel_hist->SetFillStyle(color_fill.second);
    //   sel_hist->SetLineColor(kBlack);


    //   for (size_t i = 1; i < it2->second.size(); ++i) {
    //     sel_hist->Add(it2->second.at(i));
    //   }
    //   temp_vec.push_back(sel_hist);
    //   full_mc_stack.Add(sel_hist);
    //   ++iColorFull;
    // }
    // full_mc_stack.Write();
  }
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


void protoana::AbsCexStrategy::CompareDataMC(
  const ThinSliceDistHolder & holder, 
  const ThinSliceDataSet & dataset, 
  TFile & fout) const {
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

    for (const auto & [key, hist] : holder.GetTotalSelectionHists()) {
        hist->Write();
    }
}

void protoana::AbsCexStrategy::CalcXSecs(ThinSliceDistHolder & holder, double scale) const {
    fBuilder.CalcXSecs(holder, scale);
}

double protoana::AbsCexStrategy::CalcChi2(const ThinSliceDistHolder & holder, ThinSliceDataSet & dataset) const {
    bool fDoBarlowBeeston = true;
    std::vector<int> to_skip = {5, 6};
    return fBuilder.CalcChi2(holder, dataset, fDoBarlowBeeston, to_skip);
}