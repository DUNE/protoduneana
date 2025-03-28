#ifndef THINSLICEDISTBUILDER1D_hh
#define THINSLICEDISTBUILDER1D_hh

#include "TH1.h"

#include <map>

#include "fhiclcpp/ParameterSet.h"
#include "ThinSliceDistHolder.h"
#include "ThinSliceDistBuilder.h"
namespace protoana {

class ThinSliceDistBuilder1D : public ThinSliceDistBuilder {
public:
  ThinSliceDistBuilder1D() = default;
  virtual ~ThinSliceDistBuilder1D() override = default;

  //This will have to be overriden in the derived class
  //in order to set up the distributions in the holder
  virtual void Build(
    ThinSliceDistHolder & holder,
    const fhicl::ParameterSet & pset,
    std::string label = "") override;
  virtual void CalcXSecs(
      ThinSliceDistHolder & holder,
      double scale = 1.) const override;
private:
  std::shared_ptr<TH1D> MakeSelHist(
    TString true_bins_string, const fhicl::ParameterSet & sel,
    std::string sample_name, size_t beam_bin, std::string label);
  void SelHistLoop(
    ThinSliceDistHolder & holder,
    const std::vector<fhicl::ParameterSet> & selections,
    std::string sample_name,
    TString true_bins_string, 
    int true_id, size_t beam_bin, std::string label
  );

void BuildSels(
  ThinSliceDistHolder & holder,
  const fhicl::ParameterSet & pset,
  std::string label);
void BuildXSecs(
  ThinSliceDistHolder & holder,
  const fhicl::ParameterSet & pset,
  std::string label);
};

}
#endif
