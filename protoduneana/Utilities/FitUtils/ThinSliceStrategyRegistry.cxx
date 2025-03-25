#include "ThinSliceStrategyRegistry.h"
#include "ThinSliceStrategy.h"
#include "ThinSliceStrategyFactory.h"

#include <iostream>
#include <stdexcept>

protoana::ThinSliceStrategyRegistry * protoana::ThinSliceStrategyRegistry::fInstance = 0;

protoana::ThinSliceStrategyRegistry * protoana::ThinSliceStrategyRegistry::Instance() {
  if(!fInstance) {
    static protoana::ThinSliceStrategyRegistry * manager_ptr = 0 ;
    if (!manager_ptr) {
      manager_ptr = new protoana::ThinSliceStrategyRegistry;
    }
    protoana::ThinSliceStrategyRegistry & manager = * manager_ptr;
    fInstance = & manager;
  }
  return fInstance;
}

protoana::ThinSliceStrategyRegistry::ThinSliceStrategyRegistry() {}

protoana::ThinSliceStrategyRegistry::~ThinSliceStrategyRegistry() {
  //Clean();
}

void protoana::ThinSliceStrategyRegistry::PrintAvailableStrategies() const {
  std::cout << "####ThinSliceStrategyRegistry####" << std::endl;
  std::cout << "Available Strategys:" << std::endl;
  for (auto it = fFactories.begin(); it != fFactories.end(); ++it) {
    std::cout << "Strategy: " << it->first << std::endl;
  }
  std::cout << "###############################" << std::endl << std::endl;
}

void protoana::ThinSliceStrategyRegistry::AddFactory(
    std::string name, protoana::BaseThinSliceStrategyFactory * factory) {
  fFactories[name] = factory;
}

protoana::ThinSliceStrategy * protoana::ThinSliceStrategyRegistry::GetStrategy(
    const std::string & name,
    const fhicl::ParameterSet & extra_options) {
  if (fFactories.find(name) == fFactories.end()) {
    std::string message = "Could not find ThinSliceStrategy of type: " +
                          name;
    throw std::runtime_error(message);
  }
  else {
    return fFactories[name]->Instantiate(extra_options);
  }
}
