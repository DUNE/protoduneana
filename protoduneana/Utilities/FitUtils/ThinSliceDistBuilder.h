#ifndef THINSLICEDISTBUILDER_hh
#define THINSLICEDISTBUILDER_hh

#include "TH1.h"

#include <map>

#include "fhiclcpp/ParameterSet.h"
#include "ThinSliceDistHolder.h"
namespace protoana {

class ThinSliceDistBuilder {
public:
  ThinSliceDistBuilder() = default;
  virtual ~ThinSliceDistBuilder() = default;

  //This will have to be overriden in the derived class
  //in order to set up the distributions in the holder
  virtual void Build(
    ThinSliceDistHolder & holder,
    const fhicl::ParameterSet & pset,
    std::string label = "") = 0;

  virtual void CalcXSecs(
    ThinSliceDistHolder & holder,
    double scale = 1.) const = 0; 
private:

};

}
#endif
