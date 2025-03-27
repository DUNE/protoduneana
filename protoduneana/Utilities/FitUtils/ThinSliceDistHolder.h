#ifndef THINSLICEDISTHOLDER_hh
#define THINSLICEDISTHOLDER_hh

#include "TH1.h"

#include <map>

namespace protoana {
  using TrueCat_t = int;
  using SelID_t = int;
  using TrueCatHist_map = std::map<TrueCat_t, TH1*>;
  using TrueCatSelID_t = std::pair<TrueCat_t, SelID_t>;
  using TrueCatSelID_map = std::map<TrueCatSelID_t, std::vector<TH1*>>;
class ThinSliceDistHolder {


public:
  ThinSliceDistHolder() = default;
  ~ThinSliceDistHolder() = default;
  friend class ThinSliceDistBuilder; //This will provide access to the builder
  friend class ThinSliceDistBuilder1D; //This will provide access to the builder

  void Reset();
private:

    //Each of the following is a vector corresponding the beam energy bins
        //Map of True Cat, Sel ID --> Vector of Bins (for signal) --> Hists
        std::vector<TrueCatSelID_map> fSelectionHists;
        //Map of True Cat --> Hist representing interactions, incidents, xsecs
        std::vector<TrueCatHist_map> fInteractionHists;
        std::vector<TrueCatHist_map> fIncidentHists;
        std::vector<TrueCatHist_map> fXSecHists;
};

}
#endif
