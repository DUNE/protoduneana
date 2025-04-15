#ifndef PROTOANA_THINSLICESTRATEGY_H
#define PROTOANA_THINSLICESTRATEGY_H


#include "ThinSliceEvent.h"
#include "ThinSliceDistHolder.h"
#include "ThinSliceDistBuilder.h"
#include "TH1.h"
#include "TFile.h"
#include "ThinSliceDataSet.h"
#include <mutex>

namespace protoana {

class ThinSliceStrategy {
public:

    ThinSliceStrategy(const fhicl::ParameterSet & pset) {};
    virtual ~ThinSliceStrategy() = default;

    // Pure virtual method(s) to be implemented by derived classes
    virtual void FillHistsFromEvent(const ThinSliceEvent & event, ThinSliceDistHolder & dists, size_t beam_bin, std::mutex & mutex, double weight = 1.) const = 0;
    virtual void BuildDists(ThinSliceDistHolder & holder, const fhicl::ParameterSet & pset, std::string label = "") = 0;
    virtual double GetEventWeight(const ThinSliceEvent & event) const = 0;
    virtual int GetSignalBin(const ThinSliceEvent & event, const ThinSliceDistHolder & dists) const = 0;
    virtual int GetBeamBin(
        const std::vector<double> & beam_energy_bins,
        const ThinSliceEvent & event,
        bool restrict_P) const = 0;
    virtual void CompareDataMC(const ThinSliceDistHolder & holder, const ThinSliceDataSet & dataset, TFile & fout) const = 0;
    virtual void CalcXSecs(ThinSliceDistHolder & holder, double scale = 1.) const = 0;
    virtual double CalcChi2(const ThinSliceDistHolder & holder, ThinSliceDataSet & dataset) const = 0;
    virtual void CompareSelections(
        const ThinSliceDistHolder & holder,
        ThinSliceDataSet & data_set, TFile & output_file,
        std::vector<std::pair<int, int>> plot_style,
        bool plot_rebinned,
        bool post_fit,
        TDirectory * plot_dir) const = 0;
    virtual void FillDataHistsFromEvent(
        ThinSliceDataSet & data_set, const ThinSliceEvent & event) const = 0;

    void SetMCStatVar(bool statvar) {
        fMCStatVar = statvar;
    }
    // virtual void CalcTotalDists(ThinSliceDistHolder & holder) const = 0;

protected:

    // ThinSliceDistBuilder fDistBuilder;
    bool fMCStatVar = false;
};

} // namespace protoana

#endif // PROTOANA_THINSLICESTRATEGY_H