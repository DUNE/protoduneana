#ifndef PROTOANA_THINSLICESTRATEGY_H
#define PROTOANA_THINSLICESTRATEGY_H


#include "ThinSliceEvent.h"
#include "TH1.h"

namespace protoana {

class ThinSliceStrategy {
public:

    ThinSliceStrategy() {};
    virtual ~ThinSliceStrategy() = default;

    // Pure virtual method(s) to be implemented by derived classes
    virtual void FillHistsFromEvent(const ThinSliceEvent & event) = 0;
    

protected:

    //Each of the following is a vector corresponding the beam energy bins
        //Map of True Cat, Sel ID --> Vector of Bins (for signal) --> Hists
        std::vector<std::map<std::pair<int, int>, std::vector<TH1*>>> fSelectionHists;
        //Map of True Cat --> Hist representing interactions, incidents, xsecs
        std::vector<std::map<int, TH1*>> fInteractionHists;
        std::vector<std::map<int, TH1*>> fIncidentHists;
        std::vector<std::map<int, TH1*>> fXSecHists;

    void ScaleHistsByBeamFlux(const std::vector<double> & factors);
};

} // namespace protoana

#endif // PROTOANA_THINSLICESTRATEGY_H