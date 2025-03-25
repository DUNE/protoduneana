#ifndef THINSLICEDATASET_hh
#define THINSLICEDATASET_hh

#include <vector>
#include <map>

#include "fhiclcpp/ParameterSet.h"
#include "ThinSliceSample.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom3.h"

namespace protoana {
class ThinSliceDataSet {
 public:
  ThinSliceDataSet() {};

  ThinSliceDataSet(const std::vector<fhicl::ParameterSet> & selections,
                   const std::vector<double> & beam_bins);
  ~ThinSliceDataSet() {};

  const std::map<int, std::string> & GetSelectionNames() const {
    return fSelectionNames;
  };

  std::string & GetSelectionName(int id) {
    return fSelectionNames.at(id);
  }

  std::map<int, TH1 *> & GetSelectionHists() {
    return fSelectionHists;
  };

  std::map<int, TH2D *> & GetBeamBinSelectionHists() {
    return fBeamBinSelectionHists;
  };

  std::map<int, TH1 *> & GetRebinnedSelectionHists() {
    return fSelectionHistsRebinned;
  };

  TH1 * GetSelectionHist(int id) {
    return fSelectionHists.at(id);
  };


  TH1 * GetRebinnedSelectionHist(int id) {
    return fSelectionHistsRebinned.at(id);
  };

  void SetupExtraHists(const std::vector<fhicl::ParameterSet> & extra_hists);

  void FillExtraHist(const std::string & name, double val, double weight=1.) {
    if (fExtraHists.find(name) != fExtraHists.end()) {
      fExtraHists.at(name)->Fill(val, weight);
    }
  };

  std::map<std::string, TH1 *> & GetExtraHists() {
    return fExtraHists;
  };

  //TODO PROTECT
  TH1 * GetExtraHist(std::string name) {
    if (fExtraHists.find(name) != fExtraHists.end())
      return fExtraHists.at(name);
    return 0x0;
  };

  void FillSelectionHist(int id, double val) {
    if (fSelectionHists.find(id) != fSelectionHists.end()) {
      fSelectionHists.at(id)->Fill(val);
    }
  };

  template <size_t N> void FillSelectionHist(int id, const double (& vals)[N]) {
    if (N < 1 || N > 3) {
      std::string message = "Error: trying to fill hists with too many values";
      message += std::to_string(N);
      throw std::runtime_error(message);
    }

    if (fSelectionHists.find(id) != fSelectionHists.end()) {
      if (N == 1) {
        fSelectionHists.at(id)->Fill(vals[0]);
      }
      else if (N == 2) {
        ((TH2D*)fSelectionHists.at(id))->Fill(vals[0], vals[1]);
      }
      else if (N == 3) {
        ((TH3D*)fSelectionHists.at(id))->Fill(vals[0], vals[1], vals[2]);
      }
    }
  }

  void MakeRebinnedHists();
  void Refill1DRebinned();

  void GetCumulatives() {
    fTotal = 0;
    std::cout << "Making cumulatives" << std::endl;
    for (auto it = fSelectionHists.begin();
         it != fSelectionHists.end(); ++it) {
      for (int i = 1; i <= it->second->GetNbinsX(); ++i) {
        double val = it->second->GetBinContent(i);
        if (fCumulatives.size()) {
          val += fCumulatives.back().second;
        }
        fCumulatives.push_back({{it->first, i}, val});
        fTotal += it->second->GetBinContent(i);
        std::cout << fCumulatives.back().second << " " << fTotal << std::endl;
      }
    }

    for (auto it = fCumulatives.begin(); it != fCumulatives.end(); ++it) {
      it->second /= fTotal;
    }

    std::sort(fCumulatives.begin(), fCumulatives.end(),
              [](auto a, auto b){return (a.second > b.second);});
    std::cout << "N Cumulatives: " << fCumulatives.size() << std::endl;
    for (auto c : fCumulatives) {
      std::cout << c.second << std::endl;
    }
  };

  void GenerateStatFluctuation(std::vector<double> & beam_fluxes, bool poisson=false);

  void FillHistsFromSamples(
      const std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
      double & flux, std::vector<double> & fluxes_by_beam,
      bool fluctuate=false, const std::vector<int> & to_skip = {5,6});

  void SetDirectory();
 private:
  void Rebin1D(TH1 * sel_hist, TH1 * rebinned);
  void Rebin2D(TH1 * sel_hist, TH1 * rebinned);
  void Rebin3D(TH1 * sel_hist, TH1 * rebinned);
  std::map<int, TH2D *> fBeamBinSelectionHists;
  std::map<int, TH1 *> fSelectionHists;
  std::map<int, TH1 *> fSelectionHistsRebinned;
  std::map<std::string, TH1 *> fExtraHists;
  bool fMadeRebinned = false;
  std::map<int, std::string> fSelectionNames;
  std::vector<std::pair<std::pair<int, int>, double>> fCumulatives;
  TRandom3 fRNG = TRandom3(0);
  double fTotal;

};
}
#endif
