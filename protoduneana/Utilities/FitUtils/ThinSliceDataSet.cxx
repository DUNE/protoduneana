#include "ThinSliceDataSet.h"
#include "ThinSliceSample.h"

protoana::ThinSliceDataSet::ThinSliceDataSet(
    const std::vector<double> & incident_bins,
    const std::vector<fhicl::ParameterSet> & selections) {
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

void protoana::ThinSliceDataSet::GenerateStatFluctuation() {

  //bool retry = true;
  //while (retry) {
    for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
      it->second->Reset();
      //fSelectionHistsRebinned[it->first]->Reset();
    }

    for (int i = 0; i < fTotal; ++i) {
      double r = fRNG.Uniform();
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
      //fSelectionHistsRebinned[bin.first]->AddBinContent(bin.second); 
    }

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

    Refill1DRebinned();

    //retry = !good;
  //}

 // std::cout << "Bin vals: ";
 // std::cout << bin.second << " ";
 // std::cout << std::endl;
}

void protoana::ThinSliceDataSet::FillHistsFromSamples(
    const std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    double & flux, std::vector<double> & fluxes_by_beam) {

  flux = 0.;
  for (auto it = fSelectionHists.begin(); it != fSelectionHists.end(); ++it) {
    it->second->Reset();
  }

  int a = 0;
  for (auto it = samples.begin(); it != samples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (a == 0) fluxes_by_beam.push_back(0.);
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        const auto & hists = it->second[i][j].GetSelectionHists();
        for (auto it2 = hists.begin(); it2 != hists.end(); ++it2) {
          fSelectionHists[it2->first]->Add(it2->second);
          flux += it2->second->Integral();
          fluxes_by_beam[i] += it2->second->Integral();
        }
      }
    }
    ++a;
  }
}
