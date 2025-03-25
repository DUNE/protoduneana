#ifndef THINSLICEStrategyREGISTRY_h
#define THINSLICEStrategyREGISTRY_h

#include <map>
#include "fhiclcpp/ParameterSet.h"

namespace protoana {

class ThinSliceStrategy;
class BaseThinSliceStrategyFactory;

class ThinSliceStrategyRegistry {
 public:
  static ThinSliceStrategyRegistry * Instance();
  ~ThinSliceStrategyRegistry();
  void AddFactory(std::string name, BaseThinSliceStrategyFactory * factory);
  void PrintAvailableStrategies() const;
  ThinSliceStrategy * GetStrategy(
      const std::string & name,
      const fhicl::ParameterSet & extra_options);

 private:
  ThinSliceStrategyRegistry();
  static ThinSliceStrategyRegistry * fInstance;
  std::map<std::string, BaseThinSliceStrategyFactory *> fFactories;
  
};
}
#endif
