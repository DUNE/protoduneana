#include "SelectionDefinitions.h"

#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

std::vector<double> recalc_chi2_proton::operator()(
    std::vector<std::vector<double>> dedx,
    std::vector<std::vector<double>> range) {

  protoana::ProtoDUNETrackUtils utils;
  std::vector<double> results;
  for (size_t i = 0; i < dedx.size(); ++i) {
    auto scaled_dedx = dedx[i];
    for (auto & v : scaled_dedx) v *= fScale;
    auto temp_results = utils.Chi2PID(
        scaled_dedx, range[i], fTemplate);
    results.push_back(
      (temp_results.second > 0 ?
       (temp_results.first / temp_results.second) :
       -999.
      )
    );
  }
  return results;
}
