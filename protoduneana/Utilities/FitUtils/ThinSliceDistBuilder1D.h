#ifndef THINSLICEDISTBUILDER1D_hh
#define THINSLICEDISTBUILDER1D_hh

#include "TH1.h"

#include <map>

#include "fhiclcpp/ParameterSet.h"
#include "ThinSliceDistHolder.h"
#include "ThinSliceDistBuilder.h"
namespace protoana {

class ThinSliceDistBuilder1D {
public:
  // ThinSliceDistBuilder1D() = default;
  // ~ThinSliceDistBuilder1D() override = default;

  //This will have to be overriden in the derived class
  //in order to set up the distributions in the holder
  void Build(
    ThinSliceDistHolder & holder,
    const fhicl::ParameterSet & pset);
private:

};

}
#endif
