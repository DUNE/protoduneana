#include "ThinSliceDataSet.h"
#include "ThinSliceSample.h"

protoana::ThinSliceDataSet::ThinSliceDataSet(
    const std::vector<double> & incident_bins,
    const std::vector<fhicl::ParameterSet> & selections,
    const std::vector<double> & beam_bins) {
  fIncidentHist = TH1D("Data_incident_hist",
                           "Data;Reconstructed KE (MeV)",
                           incident_bins.size() - 1,
                           &incident_bins[0]);
  for (auto it = selections.begin(); it != selections.end(); ++it) {
    fSelectionNames[it->get<int>("ID")] = it->get<std::string>("Name");
    std::string sel_name = "Data_selected_" + it->get<std::string>("Name") +
                           "_hist";
    std::vector<std::vector<double>> selected_bins =
        it->get<std::vector<std::vector<double>>>("RecoBins");

    std::vector<std::string> titles =
        it->get<std::vector<std::string>>("AxisTitles");
    TString title = "Data";
    for (auto & t : titles) {
      title += ";" + t; 
    }

    if (selected_bins.size() == 1) {
      fSelectionHists[it->get<int>("ID")] = new TH1D(
          sel_name.c_str(), title/*.c_str()"Data;Reconstructed KE (MeV)"*/,
          selected_bins[0].size() - 1, &selected_bins[0][0]);
      sel_name += "beam_bin_";
      fBeamBinSelectionHists[it->get<int>("ID")] = new TH2D(
          sel_name.c_str(), title,
          selected_bins[0].size() - 1, &selected_bins[0][0],
          beam_bins.size() - 1, &beam_bins[0]
      );
    }
    else if (selected_bins.size() == 2) {
      fSelectionHists[it->get<int>("ID")] = new TH2D(
          sel_name.c_str(), title/*.c_str()"Data;Reconstructed KE (MeV)"*/,
          selected_bins[0].size() - 1, &selected_bins[0][0],
          selected_bins[1].size() - 1, &selected_bins[1][0]);
    }
    else if (selected_bins.size() == 3) {
      fSelectionHists[it->get<int>("ID")] = new TH3D(
          sel_name.c_str(), title/*.c_str()"Data;Reconstructed KE (MeV)"*/,
          selected_bins[0].size() - 1, &selected_bins[0][0],
          selected_bins[1].size() - 1, &selected_bins[1][0],
          selected_bins[2].size() - 1, &selected_bins[2][0]);
    }
    /*else {
     * throw
     * }*/
  }


}

void protoana::ThinSliceDataSet::SetupExtraHists(
    const std::vector<fhicl::ParameterSet> & extra_hists) {
  for (const auto & hist_set : extra_hists) {
    std::string category = hist_set.get<std::string>("Category");
    std::string name = hist_set.get<std::string>("Name");
    std::string title = hist_set.get<std::string>("Title");

    bool fixed_bins = hist_set.get<bool>("DoFixedBins");

    //NEED TO MAKE SURE THIS IS ENSURED
    auto bins = hist_set.get<std::vector<double>>("Bins");
    auto binning = hist_set.get<std::vector<double>>("Binning");

    if (fixed_bins) {
      fExtraHists[category] = new TH1D(
        name.c_str(), title.c_str(), binning[0], binning[1], binning[2]
      );
    }
    else {
      fExtraHists[category] = new TH1D(
        name.c_str(), title.c_str(), bins.size()-1, &bins[0]
      );
    }
    fExtraHists[category]->SetDirectory(0);
  }
}

void protoana::ThinSliceDataSet::MakeRebinnedHists() {
  if (!fMadeRebinned) {
    std::string inc_name = fIncidentHist.GetName();
    inc_name += "Rebinned";
    fIncidentHistRebinned = TH1D(inc_name.c_str(), fIncidentHist.GetTitle(),
                                 fIncidentHist.GetNbinsX(), 0, fIncidentHist.GetNbinsX());
    for (int i = 1; i <= fIncidentHist.GetNbinsX(); ++i) {
      fIncidentHistRebinned.SetBinContent(i, fIncidentHist.GetBinContent(i));

      double low_edge = fIncidentHist.GetXaxis()->GetBinLowEdge(i);
      double up_edge = fIncidentHist.GetXaxis()->GetBinUpEdge(i);
      std::string bin_label = (low_edge < 0. ? "< 0." :
                               (protoana::PreciseToString(low_edge, 0) + " - " +
                                protoana::PreciseToString(up_edge, 0)));
      fIncidentHistRebinned.GetXaxis()->SetBinLabel(i, bin_label.c_str());
    }

    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      TH1 * sel_hist = (TH1 *)it->second;
      std::string name = sel_hist->GetName();
      name += "Rebinned";
      
      size_t nAxes = 1;
      if (sel_hist->GetNbinsY() > 1) ++nAxes;
      if (sel_hist->GetNbinsZ() > 1) ++nAxes;

      if (nAxes == 1) {
        TString title = sel_hist->GetTitle();
        title += ";";
        title += sel_hist->GetXaxis()->GetTitle();
        fSelectionHistsRebinned[it->first] = new TH1D(
            name.c_str(), title/*.c_str()sel_hist->GetTitle()*/,
            sel_hist->GetNbinsX(), 0, sel_hist->GetNbinsX());
        Rebin1D(sel_hist, fSelectionHistsRebinned[it->first]);
      }
      else if (nAxes == 2) {
        std::string title = sel_hist->GetTitle();
        title += ";";
        title += sel_hist->GetXaxis()->GetTitle();
        title += ";";
        title += sel_hist->GetYaxis()->GetTitle();

        fSelectionHistsRebinned[it->first] = new TH2D(
            name.c_str(), title.c_str()/*sel_hist->GetTitle()*/,
            sel_hist->GetNbinsX(), 0, sel_hist->GetNbinsX(),
            sel_hist->GetNbinsY(), 0, sel_hist->GetNbinsY());
        Rebin2D(sel_hist, fSelectionHistsRebinned[it->first]);
      }
      else if (nAxes == 3) {
        std::string title = sel_hist->GetTitle();
        title += ";";
        title += sel_hist->GetXaxis()->GetTitle();
        title += ";";
        title += sel_hist->GetYaxis()->GetTitle();
        title += ";";
        title += sel_hist->GetZaxis()->GetTitle();

        fSelectionHistsRebinned[it->first] = new TH3D(
            name.c_str(), title.c_str()/*sel_hist->GetTitle()*/,
            sel_hist->GetNbinsX(), 0, sel_hist->GetNbinsX(),
            sel_hist->GetNbinsY(), 0, sel_hist->GetNbinsY(),
            sel_hist->GetNbinsZ(), 0, sel_hist->GetNbinsZ());
        Rebin3D(sel_hist, fSelectionHistsRebinned[it->first]);
      }
    }

    fMadeRebinned = true;
  }
}

void protoana::ThinSliceDataSet::Refill1DRebinned() {

  if (!fMadeRebinned) MakeRebinnedHists();

  for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
    for (int i = 1; i <= it->second->GetNbinsX(); ++i) {
      fSelectionHistsRebinned[it->first]->SetBinContent(
          i, it->second->GetBinContent(i));
    }
  }
}

void protoana::ThinSliceDataSet::Rebin1D(TH1 * sel_hist, TH1 * rebinned) {
  for (int i = 1; i <= sel_hist->GetNbinsX(); ++i) {
    double low_x = sel_hist->GetXaxis()->GetBinLowEdge(i);
    double up_x = sel_hist->GetXaxis()->GetBinUpEdge(i);
    std::string bin_label = (low_x < 0. ? "< 0." :
                             (protoana::PreciseToString(low_x, 0) + " - " +
                              protoana::PreciseToString(up_x, 0)));
    rebinned->GetXaxis()->SetBinLabel(i, bin_label.c_str());

    rebinned->SetBinContent(i, sel_hist->GetBinContent(i));
  }
}

void protoana::ThinSliceDataSet::Rebin2D(TH1 * sel_hist, TH1 * rebinned) {
  for (int i = 1; i <= sel_hist->GetNbinsX(); ++i) {
    double low_x = sel_hist->GetXaxis()->GetBinLowEdge(i);
    double up_x = sel_hist->GetXaxis()->GetBinUpEdge(i);
    std::string bin_label = (low_x < 0. ? "< 0." :
                             (protoana::PreciseToString(low_x, 0) + " - " +
                              protoana::PreciseToString(up_x, 0)));
    rebinned->GetXaxis()->SetBinLabel(
        i, bin_label.c_str());
    for (int j = 1; j <= sel_hist->GetNbinsY(); ++j) {
      double low_y = sel_hist->GetYaxis()->GetBinLowEdge(j);
      double up_y = sel_hist->GetYaxis()->GetBinUpEdge(j);
      std::string y_label = (low_y < 0. ? "< 0." :
                             (protoana::PreciseToString(low_y, 0) + " - " +
                              protoana::PreciseToString(up_y, 0)));
      rebinned->GetYaxis()->SetBinLabel(j, bin_label.c_str());
      rebinned->SetBinContent(i, j, sel_hist->GetBinContent(i, j));
    }
  }
}

void protoana::ThinSliceDataSet::Rebin3D(TH1 * sel_hist, TH1 * rebinned) {
  for (int i = 1; i <= sel_hist->GetNbinsX(); ++i) {
    double low_x = sel_hist->GetXaxis()->GetBinLowEdge(i);
    double up_x = sel_hist->GetXaxis()->GetBinUpEdge(i);
    std::string bin_label = (low_x < 0. ? "< 0." :
                             (protoana::PreciseToString(low_x, 0) + " - " +
                              protoana::PreciseToString(up_x, 0)));
    rebinned->GetXaxis()->SetBinLabel(i, bin_label.c_str());
    for (int j = 1; j <= sel_hist->GetNbinsY(); ++j) {
      double low_y = sel_hist->GetYaxis()->GetBinLowEdge(j);
      double up_y = sel_hist->GetYaxis()->GetBinUpEdge(j);
      std::string y_label = (low_y < 0. ? "< 0." :
                             (protoana::PreciseToString(low_y, 0) + " - " +
                              protoana::PreciseToString(up_y, 0)));
      rebinned->GetYaxis()->SetBinLabel(j, bin_label.c_str());

      for (int k = 1; k <= sel_hist->GetNbinsY(); ++k) {
        double low_z = sel_hist->GetYaxis()->GetBinLowEdge(k);
        double up_z = sel_hist->GetYaxis()->GetBinUpEdge(k);
        std::string y_label = (low_z < 0. ? "< 0." :
                               (protoana::PreciseToString(low_z, 0) + " - " +
                                protoana::PreciseToString(up_z, 0)));
        rebinned->GetZaxis()->SetBinLabel(k, bin_label.c_str());

        rebinned->SetBinContent(i, j, k, sel_hist->GetBinContent(i, j, k));
      }
    }
  }
}

void protoana::ThinSliceDataSet::GenerateStatFluctuation(
    std::vector<double> & beam_fluxes, bool poisson) {

  std::cout << "Beam Fluxes" << std::endl;
  for (auto & bf : beam_fluxes) {
    std::cout << bf << std::endl;
    bf = 0.;
  }
  //bool retry = true;
  //while (retry) {
  //if (!poisson) {
    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      it->second->Reset();
      //fSelectionHistsRebinned[it->first]->Reset();
    }

    //int total = 0;
    std::cout << "Seed: " << fRNG.GetSeed() << std::endl;
    /*
    for (int i = 0; i < fTotal; ++i) {
      double r = fRNG.Uniform();
      if (i < 10) std::cout << i << " " << r << std::endl;
      std::pair<int, int> bin;
      for (size_t j = 0; j < fCumulatives.size(); ++j) {
        //std::cout << fCumulatives[j].second << " " <<  r <<
        //             " "  << fCumulatives[j].second - r << std::endl;
        if ((fCumulatives[j].second - r) > 0.) {
          bin = fCumulatives[j].first;
        }
        else {
          break;
        }
      }
      //std::cout << "Found bin: " << bin.first << " " << bin.second << std::endl;
      fSelectionHists[bin.first]->AddBinContent(bin.second); 
      ++total;
      //fSelectionHistsRebinned[bin.first]->AddBinContent(bin.second); 
    }

    double bin_total = 0.;
    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      for (int i = 1; i <= it->second->GetNbinsX(); ++i) {
        bin_total += it->second->GetBinContent(i);
      }
    }
    std::cout << total << " " << bin_total << std::endl;*/

  //}
  /*else {
    std::cout << "Poisson fluc" << std::endl;
    double total = 0.;
    for (auto & it : fSelectionHists) {
      for (int i = 1; i <= it.second->GetNbinsX(); ++i) {
        std::cout << it.first << " " << it.second << " " <<
                     it.second->GetBinContent(i) << " ";
        it.second->SetBinContent(
            i, fRNG.PoissonD(it.second->GetBinContent(i)));
        std::cout << it.second->GetBinContent(i) << std::endl;
        total += it.second->GetBinContent(i);
      }
    }
    std::cout << "total: " << total << std::endl;
  }*/

    //bool good = true;
    /*
    double new_total = 0.;
    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      for (int i = 1; i <= it->second->GetNbinsX(); ++i) {
        if (it->second->GetBinContent(i) < 1.) {
          std::cout << "DataSet fluctuation: Found bin with 0 content. Adding" << std::endl;
          it->second->AddBinContent(i);
          //fSelectionHistsRebinned[it->first]->AddBinContent(i);
          //good = false;
        }
        new_total += it->second->GetBinContent(i);
      }
    }

    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      it->second->Scale(new_total/fTotal);
      //fSelectionHistsRebinned[it->first]->Scale(new_total/fTotal);
    }
    */

    for (auto & it : fBeamBinSelectionHists) {
      auto * sel_hist = fSelectionHists[it.first];
      for (int i = 1; i <= it.second->GetNbinsX(); ++i) {
        for (int j = 1; j <= it.second->GetNbinsY(); ++j) {
          double val = it.second->GetBinContent(i, j);
          double r_val = fRNG.PoissonD(val);
          sel_hist->AddBinContent(i, r_val);


          beam_fluxes[j-1] += r_val;
        }
      }
    }
    std::cout << "New Beam Fluxes" << std::endl;
    for (auto & bf : beam_fluxes) {
      std::cout << bf << std::endl;
    }

    Refill1DRebinned();

    //retry = !good;
  //}

 // std::cout << "Bin vals: ";
 // std::cout << bin.second << " ";
 // std::cout << std::endl;
}

//TODO chunk this out and put in the driver. Then add parts for the extra hists
void protoana::ThinSliceDataSet::FillHistsFromSamples(
    const std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    double & flux, std::vector<double> & fluxes_by_beam, bool fluctuate,
    const std::vector<int> & to_skip) {

  flux = 0.;
  //std::vector<double> temp_flux(/*fluxes_by_beam.size()*/);
  fluxes_by_beam.clear();
  for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
    it->second->Reset();
  }


  if (fluctuate) {
    //Get total
    double total = 0.;
    int a = 0;
    for (auto & sample : samples ) {
      for (size_t i = 0; i < sample.second.size(); ++i) {
        for (size_t j = 0; j < sample.second[i].size(); ++j) {
          const auto & hists = sample.second[i][j].GetSelectionHists();
          for (auto & hist : hists) {
            total += hist.second->Integral();
            std::cout << a << " " << i << " "  << j << " " << hist.second->Integral() << std::endl;
            //std::cout << temp_flux.size() << std::endl;
            //temp_flux[i] += hist.second->Integral();
          }
        }
      }
      ++a;
    }
    std::cout << "Total: " << total << std::endl;
    //for (auto & f : temp_flux) {std::cout << f << " ";} std::cout << std::endl;


    //First, iterate through all samples/hists. Vary the contents
    //Skip accordingly. Get varied total
    double varied_total = 0.;
    double nominal_skip_total = 0.;
    //std::map<int, double> nominal_skip_vals;
    //for (auto i : to_skip) nominal_skip_vals[i] = 0.;

    for (auto & sample : samples ) {
      auto & sample_vec = sample.second;
      for (size_t i = 0; i < sample_vec.size(); ++i) {
        for (size_t j = 0; j < sample_vec[i].size(); ++j) {
          const auto & hists = sample_vec[i][j].GetSelectionHists();
          for (auto & hist : hists) {
            if (std::find(to_skip.begin(), to_skip.end(), hist.first)
                != to_skip.end()) {
              nominal_skip_total += hist.second->Integral();
              //nominal_skip_vals[hist.first] += hist.second->Integral();
              continue;
            }
            for (int k = 1; k <= hist.second->GetNbinsX(); ++k) {
              double val = hist.second->GetBinContent(k);
              double rand = fRNG.PoissonD(val);
              hist.second->SetBinContent(k, rand);
            }
            varied_total += hist.second->Integral();
          }
        }
      }
    }

    double leftover = total - varied_total;
    double ratio = leftover/nominal_skip_total;
    std::cout << "Leftover: " << leftover << std::endl;
    std::cout << "Nominal skip total " << nominal_skip_total << std::endl;

    for (auto & sample : samples ) {
      auto & sample_vec = sample.second;
      for (size_t i = 0; i < sample_vec.size(); ++i) {
        for (size_t j = 0; j < sample_vec[i].size(); ++j) {
          const auto & hists = sample_vec[i][j].GetSelectionHists();
          for (auto k : to_skip) {
            if (hists.find(k) == hists.end()) continue;
            //std::cout << "Scaling " << k << " " << hists.at(k)->Integral() <<
            //             " " << ratio << std::endl;
            hists.at(k)->Scale(ratio);
          }
        }
      }
    }
  }

  for (auto & extra_hist : fExtraHists) {
    extra_hist.second->Reset();
  }

  int a = 0;
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (a == 0) fluxes_by_beam.push_back(0.);
      for (size_t j = 0; j < it->second[i].size(); ++j) {

        const auto & hists = it->second[i][j].GetSelectionHists();
        for (auto it2 = hists.begin(); it2 != hists.end(); ++it2) {

          fSelectionHists[it2->first]->Add(it2->second);
          for (int sel_bin = 1; sel_bin <= it2->second->GetNbinsX(); ++sel_bin) {
            double x
                = fBeamBinSelectionHists[it2->first]->GetXaxis()->GetBinCenter(
                    sel_bin);
            double y
                = fBeamBinSelectionHists[it2->first]->GetYaxis()->GetBinCenter(
                    i+1);
            fBeamBinSelectionHists[it2->first]->Fill(
                x, y, it2->second->GetBinContent(sel_bin));
          }
          flux += it2->second->Integral();
          fluxes_by_beam[i] += it2->second->Integral();
          //std::cout << "Adding " << it2->second->Integral() << " to " << i << std::endl;
        }

        for (auto & extra_hists : it->second[i][j].GetExtraHists()) {
          fExtraHists[extra_hists.first]->Add(extra_hists.second);
        }
      }
    }
    ++a;
  }
  std::cout << "Flux: " << flux << std::endl;
}

void protoana::ThinSliceDataSet::SetDirectory() {
  for (auto & hist : fSelectionHists) {
    hist.second->SetDirectory(0);
  }
  for (auto & hist : fSelectionHistsRebinned) {
    hist.second->SetDirectory(0);
  }

  for (auto & hist : fExtraHists) {
    hist.second->SetDirectory(0);
  }
}
