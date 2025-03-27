#include "ThinSliceDistBuilder1D.h"
void protoana::ThinSliceDistBuilder1D::Build(
    ThinSliceDistHolder & holder,
    const fhicl::ParameterSet & pset) {
  return;

  const auto & beam_energy_bins = 
      pset.get<std::vector<double>>("BeamEnergyBins");
      holder.fSelectionHists = std::vector<TrueCatSelID_map>(beam_energy_bins.size()-1);
      holder.fInteractionHists = std::vector<TrueCatHist_map>(beam_energy_bins.size()-1);
      holder.fIncidentHists = std::vector<TrueCatHist_map>(beam_energy_bins.size()-1);
      holder.fXSecHists = std::vector<TrueCatHist_map>(beam_energy_bins.size()-1);
  
  const auto & selections = pset.get<std::vector<fhicl::ParameterSet>>("Selections");
  const auto & samples = pset.get<std::vector<fhicl::ParameterSet>>("Samples");
  // const auto & measurement_samples = pset.get<int>("MeasurementSamples");
  
  for (size_t i = 1; i < beam_energy_bins.size(); ++i) {
    for (const auto & sample : samples) {

      std::string sample_name = sample.get<std::string>("Name");
      int true_id = sample.get<int>("ID");
      bool is_signal = sample.get<bool>("IsSignal");
      bool is_single_signal_bin = sample.get<bool>("SingleSignalBin", false);
      bool has_bins = is_signal && !is_single_signal_bin;
      std::vector<double> signal_bins;
      if (has_bins) {
        signal_bins = sample.get<std::vector<double>>("SignalBins");
      }
      else {
        signal_bins = {0., 0.};
      }
        

      //Make one entry in the map for each true bin
      if (is_signal && is_single_signal_bin) {

      }
      //Just do a single one
      else {

      }


      for (size_t j = 1; j < signal_bins.size()-1; ++j) {
        for (const auto & sel : selections) {
          int id = sel.get<int>("ID");
          std::string sel_name = sel.get<std::string>("Name");

          //Assume these have at least one
          //TODO -- make these single axes
          std::string title = sel.get<std::vector<std::string>>("AxisTitles")[0];
          std::vector<double> bins = sel.get<std::vector<std::vector<double>>>(
            "RecoBins")[0];


          std::string hist_name = "asdf";
          std::string hist_title = "asdf";
          holder.fSelectionHists[i-1][{true_id, id}].push_back(
            new TH1D(
              sel_name.c_str(), title.c_str(),
              bins.size()-1, &bins[0])
            );
          
        }
      }
    }
  }

}