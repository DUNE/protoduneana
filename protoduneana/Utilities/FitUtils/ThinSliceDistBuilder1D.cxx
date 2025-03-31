#include "ThinSliceDistBuilder1D.h"
#include "cetlib_except/exception.h"
#include "TString.h"
#include "TH1D.h"


std::shared_ptr<TH1D> protoana::ThinSliceDistBuilder1D::MakeSelHist(
      TString true_bins_string,
      TString title,
      const fhicl::ParameterSet & sel,
      std::string sample_name, size_t beam_bin, std::string label,
      bool total) {
    std::string sel_name = sel.get<std::string>("Name");

    //Assume these have at least one
    //TODO -- make these single axes
    // std::string axis_title = sel.get<std::vector<std::string>>("AxisTitles")[0];
    // TString title;
    // title.Form("%s_")
    // (is_signal ?
    //     ("(" + protoana::PreciseToString(range.first) + " "
    //      + protoana::PreciseToString(range.second) + ")") :
    //      "") +
    // ";Reconstructed KE (MeV)";



    std::vector<double> reco_bins = sel.get<std::vector<std::vector<double>>>(
      "RecoBins")[0];

    TString hist_name;
    if (!total) {
      hist_name.Form(
        "sample_%s_%s_%s_selected_%s_hist_beam_%lu",
        sample_name.c_str(), label.c_str(), true_bins_string.Data(),
        sel_name.c_str(), beam_bin-1
      );
    }
    else {
      hist_name.Form(
        "total_%s_selected_%s_hist",
        label.c_str(), sel_name.c_str()
      );
    }

    // holder.fSelectionHists[i-1][{true_id, id}].push_back(
    std::cout << "Making " << hist_name << std::endl;
    return std::make_shared<TH1D>(
      hist_name, title,
      reco_bins.size()-1, &reco_bins[0]
    );
}

void protoana::ThinSliceDistBuilder1D::SelHistLoop(
  ThinSliceDistHolder & holder,
  const std::vector<fhicl::ParameterSet> & selections,
  std::string sample_name,
  TString true_bins_string, TString title, int true_id, size_t beam_bin, std::string label
) {
  for (const auto & sel : selections) {
    int sel_id = sel.get<int>("ID");
    std::string sel_name = sel.get<std::string>("Name");

    holder.fSelectionHists[beam_bin-1][{true_id, sel_id}].push_back(
      MakeSelHist(true_bins_string, title, sel, sample_name, beam_bin, label)
    );
  }
}

void protoana::ThinSliceDistBuilder1D::TotalSelHistLoop(
  ThinSliceDistHolder & holder,
  const std::vector<fhicl::ParameterSet> & selections, std::string label) {
  for (const auto & sel : selections) {
    int sel_id = sel.get<int>("ID");
    std::string sel_name = sel.get<std::string>("Name");

    holder.fTotalSelectionHists[sel_id] =
      MakeSelHist(TString(), TString(), sel, "", 0, label, true);
  }
}

void protoana::ThinSliceDistBuilder1D::BuildXSecs(
  ThinSliceDistHolder & holder,
  const fhicl::ParameterSet & pset,
  std::string label) {
    
    const auto & beam_energy_bins = 
    pset.get<std::vector<double>>("BeamEnergyBins");
    int nbins = beam_energy_bins.size()-1;
  
    holder.fInteractionHists = std::vector<TrueCatHist_map>(nbins);
    holder.fIncidentHists = std::vector<TrueCatHist_map>(nbins);
    
    const auto & measurement_sample_IDs = pset.get<std::vector<int>>("MeasurementSamples");
    const auto & incident_sample_IDs = pset.get<std::vector<int>>("IncidentSamples");
    holder.fIncidentIDs = incident_sample_IDs;
    holder.fMeasurementIDs = measurement_sample_IDs;

    const auto & samples = pset.get<std::vector<fhicl::ParameterSet>>("Samples");

    for (size_t i = 1; i < beam_energy_bins.size(); ++i) {
      for (const auto & sample : samples) {
        const auto & signal_bins = sample.get<std::vector<double>>("SignalBins", {});
        int true_id = sample.get<int>("ID");
        holder.fTrueCatIDs.push_back(true_id);

        if (std::find(measurement_sample_IDs.begin(), measurement_sample_IDs.end(), true_id) != measurement_sample_IDs.end()) {
          std::string sample_name = sample.get<std::string>("Name");
          std::vector<double> signal_bins = sample.get<std::vector<double>>("SignalBins", {});
          TString hist_name;
          //Making the incident histogram for this measured sample
          hist_name.Form(
            "sample_%s_%s_incident_hist_beam_%lu", sample_name.c_str(), label.c_str(), (i-1)
          );
          std::cout << "Making " << hist_name << std::endl;

          holder.fIncidentHists[i-1][true_id] =
            std::make_shared<TH1D>(
              hist_name, "",
              signal_bins.size()-1, &signal_bins[0]
          );

          //Only need to do 1 for each xsec
          if (i == 1) {
            //Making the xsec histogram for this measured sample
            hist_name.Form(
              "sample_%s_%s_xsec_hist", sample_name.c_str(), label.c_str()
            );
            std::cout << "Making " << hist_name << std::endl;

            holder.fXSecHists[true_id] =
              std::make_shared<TH1D>(
                hist_name, "",
                signal_bins.size()-1, &signal_bins[0]
            );

            hist_name.Form(
              "sample_%s_%s_total_incident_hist", sample_name.c_str(), label.c_str()
            );
            holder.fTotalIncidentHists[true_id] =
            std::make_shared<TH1D>(
              hist_name, "",
              signal_bins.size()-1, &signal_bins[0]
            );
          }


          //Making the interaction histogram for this measured sample
          hist_name.Form(
            "sample_%s_%s_interaction_hist_beam_%lu", sample_name.c_str(), label.c_str(), (i-1)
          );
          std::cout << "Making " << hist_name << std::endl;

          holder.fInteractionHists[i-1][true_id] =
            std::make_shared<TH1D>(
              hist_name, "",
              signal_bins.size()-1, &signal_bins[0]
          );
        }
    }
  }
}

void protoana::ThinSliceDistBuilder1D::BuildSels(
  ThinSliceDistHolder & holder,
  const fhicl::ParameterSet & pset,
  std::string label) {

    const auto & beam_energy_bins = 
    pset.get<std::vector<double>>("BeamEnergyBins");

    //Check that the beam energy bins are valid
    if (beam_energy_bins.size() < 2) {
      throw cet::exception("ThinSliceDistBuilder1D")
        << "Beam energy bins must have at least two entries";
    }
    int nbins = beam_energy_bins.size()-1;
    holder.fSelectionHists = std::vector<TrueCatSelID_map>(nbins);


    const auto & selections = pset.get<std::vector<fhicl::ParameterSet>>("Selections");
    const auto & samples = pset.get<std::vector<fhicl::ParameterSet>>("Samples");
    // const auto & measurement_samples = pset.get<int>("MeasurementSamples");

    for (size_t i = 1; i < beam_energy_bins.size(); ++i) {
      for (const auto & sample : samples) {

        std::string sample_name = sample.get<std::string>("Name");
        int true_id = sample.get<int>("ID");
        bool is_signal = sample.get<bool>("IsSignal");
        if (is_signal) {
          holder.fSignalIDs.push_back(true_id);
        }
        bool is_single_signal_bin = sample.get<bool>("SingleSignalBin", false);
        // bool has_bins = is_signal && !is_single_signal_bin;
        std::vector<double> signal_bins = sample.get<std::vector<double>>("SignalBins", {});

        //If it's signal and has multiple bins add underflow and overflow
        // then loop over the true bins
        TString title;
        if (is_signal && !is_single_signal_bin) {
          //Underflow
          TString true_bins_string("Underflow");
          title.Form("%s Underflow", sample_name.c_str());
          SelHistLoop(
            holder,
            selections, sample_name, true_bins_string,
            title,
            true_id, i,
            label
          );

          //Loop over signal bins, and make a selection histogram for each
          for (size_t j = 1; j < signal_bins.size(); ++j) {
            TString true_bins_string;
            true_bins_string.Form(
              "%.2f_%.2f",
              signal_bins[j-1], signal_bins[j]
            );
            title.Form("%s (%.2f %.2f)", sample_name.c_str(), signal_bins[j-1], signal_bins[j]);
            SelHistLoop(
              holder,
              selections, sample_name, true_bins_string,
              title,
              true_id, i,
              label
            );
          }
          //Overflow
          true_bins_string.Form("Overflow");
          title.Form("%s Underflow", sample_name.c_str());
          SelHistLoop(
            holder,
            selections, sample_name, true_bins_string,
            title,
            true_id, i, label
          );
        }
        else if (is_signal && is_single_signal_bin) {
          TString true_bins_string("Single");
          title.Form("%s", sample_name.c_str());
          SelHistLoop(
            holder,
            selections, sample_name, true_bins_string,
            title,
            true_id, i, label
          );
        }
        else {
          TString true_bins_string("");
          title.Form("%s", sample_name.c_str());
          SelHistLoop(
            holder,
            selections, sample_name, true_bins_string,
            title,
            true_id, i, label
          );
        }
      }
    }

    TotalSelHistLoop(holder, selections, label);
}

void protoana::ThinSliceDistBuilder1D::Build(
    ThinSliceDistHolder & holder,
    const fhicl::ParameterSet & pset,
    std::string label) {
  BuildSels(holder, pset, label);
  BuildXSecs(holder, pset, label);
}


  
void protoana::ThinSliceDistBuilder1D::CalcXSecs(
  ThinSliceDistHolder & holder,
  double scale) const {
    //Loop over measurements
    // for each, sum interaction and incident contribution from each beam bin
    for (auto & id : holder.fMeasurementIDs) {
      auto xsec_hist = holder.fXSecHists.at(id);
      xsec_hist->Reset();
      auto total_incident = holder.fTotalIncidentHists.at(id);
      total_incident->Reset();

      //Add each interaction to the xsec and sum the incidents
      for (auto & hists : holder.fInteractionHists) {
        xsec_hist->Add(hists.at(id).get());
      }
      for (auto & hists : holder.fIncidentHists) {
        total_incident->Add(hists.at(id).get());
      }

      //Now go through each bin and calculate the cross section per bin
      //Defined as -1.*log(1. - interaction/incident));
      for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
        double interactions = xsec_hist->GetBinContent(i);
        double incidents = total_incident->GetBinContent(i);
        xsec_hist->SetBinContent(i, -1.*log(1. - interactions/incidents));
      }
      xsec_hist->Scale(scale);
    }
}

double protoana::ThinSliceDistBuilder1D::CalcChi2(
  const ThinSliceDistHolder & holder, ThinSliceDataSet & dataset,
  bool do_barlow_beeston,
  std::vector<int> to_skip) const {
    double chi2 = 0.;

    //Iterate over selections
    for (const auto & [sel_ID, mc_hist] : holder.fTotalSelectionHists) {
      if (std::find(to_skip.begin(), to_skip.end(), sel_ID) == to_skip.end())
        continue;

      auto * data_hist = dataset.GetSelectionHist(sel_ID);
      for (int i = 1; i <= mc_hist->GetNbinsX(); ++i) {
        double data_val = data_hist->GetBinContent(i);
        /// Skip any bins with data == 0
        //
        //See PDG Stat Review:
        //https://pdg.lbl.gov/2018/reviews/rpp2018-rev-statistics.pdf
        //Page 6
        //
        if (data_val < 1.e-7) continue;

        double mc_val = mc_hist->GetBinContent(i);
        double mc_sumw2 = std::pow(mc_hist->GetBinError(i), 2);
        double sigma_squared = (mc_val > 1.e-7 ? mc_sumw2/(mc_val*mc_val) : 1.);

        double beta = .5*((1. - mc_val*sigma_squared) +
                   sqrt(std::pow((1. - mc_val*sigma_squared), 2) +
                        4.*data_val*sigma_squared));
        chi2 += 2.*((do_barlow_beeston ? beta : 1.)*mc_val - data_val +
        data_val*std::log(data_val/(mc_val*(do_barlow_beeston ? beta : 1.)))) +
        (do_barlow_beeston ? ((beta - 1.)*(beta - 1.)/sigma_squared) : 0.);
      }
    }
    return chi2;
}