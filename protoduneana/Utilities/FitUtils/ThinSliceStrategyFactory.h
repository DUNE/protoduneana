#ifndef THINSLICESTRATEGYFACTORY_H
#define THINSLICESTRATEGYFACTORY_H

#include "ThinSliceStrategyRegistry.h"
#include "ThinSliceStrategy.h"

#include <string>

namespace protoana {
class BaseThinSliceStrategyFactory {
 public:
  virtual ThinSliceStrategy * Instantiate(
      const fhicl::ParameterSet & extra_options) = 0;
};

template <typename T> class ThinSliceStrategyFactory
    : public BaseThinSliceStrategyFactory {
 public:

  ThinSliceStrategyFactory(const std::string name) {
    ThinSliceStrategyRegistry::Instance()->AddFactory(name, this);
  }

  virtual ThinSliceStrategy * Instantiate(
      const fhicl::ParameterSet & extra_options) {
    return new T(/*extra_options*/);
  }
};
}

#define DECLARE_THINSLICESTRATEGY_FACTORY(strategy) \
  const ThinSliceStrategyFactory<strategy>& strategy##Factory = ThinSliceStrategyFactory<strategy>(#strategy)

// support for strategys  defined within a namespace
// a bit tricky due to cpp macro expansion and the use of "::"
// use  DECLARE_THINSLICESTRATEGY_FACTORY_NS( myns::MyStrategy, myns, strategybase )  // without trailing ";"
#define DECLARE_THINSLICESTRATEGY_FACTORY_NS( strategy, nsname, strategybase )  \
  namespace nsname { \
    const ThinSliceStrategyFactory<strategy>& strategybase##Factory = ThinSliceStrategyFactory<strategy>(#strategy); \
  }

#endif
