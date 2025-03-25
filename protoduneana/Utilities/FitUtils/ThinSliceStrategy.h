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


    void ScaleHistsByBeamFlux(const std::vector<double> & factors);
};

} // namespace protoana

#endif // PROTOANA_THINSLICESTRATEGY_H