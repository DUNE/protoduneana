#include "ThinSliceDistHolder.h"

void protoana::ThinSliceDistHolder::Reset() {
    for (auto & hist_map : fSelectionHists) {
        for (auto & [key, vec] : hist_map) {
            for (auto & hist : vec) {
                hist->Reset();
            }
        }
    }

    for (auto & hist_map : fInteractionHists) {
        for (auto & [key, val] : hist_map) {
            val->Reset();
        }
    }
    for (auto & hist_map : fIncidentHists) {
        for (auto & [key, val] : hist_map) {
            val->Reset();
        }
    }
    for (auto & [key, val] : fXSecHists) {
        val->Reset();
        fTotalIncidentHists[key]->Reset();
    }

}

void protoana::ThinSliceDistHolder::SetIDs(const fhicl::ParameterSet & pset) {
    const auto & selections = pset.get<std::vector<fhicl::ParameterSet>>("Selections");
    const auto & samples = pset.get<std::vector<fhicl::ParameterSet>>("Samples");

    for (const auto & sample : samples) {
        int true_id = sample.get<int>("ID");
        fTrueCatIDs.push_back(true_id);
        bool is_signal = sample.get<bool>("IsSignal");
        if (is_signal) {
            fSignalIDs.push_back(true_id);
        }

        bool is_single_bin_signal = sample.get<bool>("SingleSignalBin", false);
        if (!is_single_bin_signal) {
            std::vector<double> signal_bins = sample.get<std::vector<double>>("SignalBins", {});
            for (size_t i = 1; i < signal_bins.size(); ++i) {
                fSignalRanges[true_id].emplace_back(signal_bins[i-1], signal_bins[i]);
            }
        }
    }

    for (const auto & sel : selections) {
        int sel_id = sel.get<int>("ID");
        fSelIDs.push_back(sel_id);
    }
}