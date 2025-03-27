#ifndef PROTOANA_THINSLICESTRATEGY_H
#define PROTOANA_THINSLICESTRATEGY_H


#include "ThinSliceEvent.h"
#include "ThinSliceDistHolder.h"
#include "TH1.h"

namespace protoana {

class ThinSliceStrategy {
public:

    ThinSliceStrategy() {};
    virtual ~ThinSliceStrategy() = default;

    // Pure virtual method(s) to be implemented by derived classes
    virtual void FillHistsFromEvent(const ThinSliceEvent & event, ThinSliceDistHolder & dists) = 0;
    

protected:


    
};

} // namespace protoana

#endif // PROTOANA_THINSLICESTRATEGY_H