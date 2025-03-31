#ifndef PROTOANA_ABSCEXSTRATEGY_H
#define PROTOANA_ABSCEXSTRATEGY_H

#include "ThinSliceStrategy.h"

#include "ThinSliceDistBuilder1D.h"

namespace protoana {

class AbsCexStrategy : public ThinSliceStrategy {
public:
    AbsCexStrategy(const fhicl::ParameterSet & pset);
    ~AbsCexStrategy() override = default;

    // Implementation of the pure virtual method from ThinSliceStrategy
    virtual void FillHistsFromEvent(const ThinSliceEvent & event, ThinSliceDistHolder & dists, size_t beam_bin, std::mutex & mutex, double weight = 1.) const override;
    virtual void BuildDists(ThinSliceDistHolder & holder, const fhicl::ParameterSet & pset, std::string label = "") override;
    virtual double GetEventWeight(const ThinSliceEvent & event) const override;
    virtual int GetSignalBin(const ThinSliceEvent & event, const ThinSliceDistHolder & dists) const override;
    virtual int GetBeamBin(
        const std::vector<double> & beam_energy_bins,
        const ThinSliceEvent & event,
        bool restrict_P) const override;
    virtual void CompareDataMC(const ThinSliceDistHolder & holder, const ThinSliceDataSet & dataset, TFile & fout) const override;
    virtual void CalcXSecs(ThinSliceDistHolder & holder, double scale = 1.) const override;
    virtual double CalcChi2(const ThinSliceDistHolder & holder, ThinSliceDataSet & dataset) const override;

    virtual void CompareSelections(
        const ThinSliceDistHolder & holder,
        ThinSliceDataSet & data_set, TFile & output_file,
        std::vector<std::pair<int, int>> plot_style,
        bool plot_rebinned,
        bool post_fit,
        TDirectory * plot_dir) const override;
    // virtual void CalcTotalDists(ThinSliceDistHolder & holder) const override;

private:
    std::vector<int> fERecoSelections, fEndZSelections, fOneBinSelections;
    double GetTrueEndEnergy(const ThinSliceEvent & event) const;
    std::vector<double> MakeTrueIncidentEnergies(const std::vector<double> & true_beam_traj_Z,
        const std::vector<double> & true_beam_traj_KE) const;
    ThinSliceDistBuilder1D fBuilder;

    std::pair<int, int> GetColorAndStyle (
        size_t i, const std::vector<std::pair<int, int>> & plot_style) const {
          return {plot_style.at(i % plot_style.size()).first,
                  (i < plot_style.size() ? 1001: 3244)};
    };

    double fPitch;
    int fSliceCut;
    double fTrajZStart;

};

} // namespace protoana

#endif // PROTOANA_ABSCEXSTRATEGY_H
