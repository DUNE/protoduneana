#ifndef THINSLICEDISTHOLDER_hh
#define THINSLICEDISTHOLDER_hh

#include "TH1.h"

#include <map>

#include "fhiclcpp/ParameterSet.h"

namespace protoana {

class ThinSliceDistHolder {
public:
  ThinSliceDistHolder() = default;

  ~ThinSliceDistHolder() = default;
private:

    //Each of the following is a vector corresponding the beam energy bins
        //Map of True Cat, Sel ID --> Vector of Bins (for signal) --> Hists
        std::vector<std::map<std::pair<int, int>, std::vector<TH1*>>> fSelectionHists;
        //Map of True Cat --> Hist representing interactions, incidents, xsecs
        std::vector<std::map<int, TH1*>> fInteractionHists;
        std::vector<std::map<int, TH1*>> fIncidentHists;
        std::vector<std::map<int, TH1*>> fXSecHists;
};

}
#endif
