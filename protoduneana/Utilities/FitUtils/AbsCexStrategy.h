#ifndef PROTOANA_ABSCEXSTRATEGY_H
#define PROTOANA_ABSCEXSTRATEGY_H

#include "ThinSliceStrategy.h"


namespace protoana {

class AbsCexStrategy : public ThinSliceStrategy {
public:
    AbsCexStrategy() = default;
    ~AbsCexStrategy() override = default;

    // Implementation of the pure virtual method from ThinSliceStrategy
    void FillHistsFromEvent(const ThinSliceEvent & event) const override {
        // Add logic to fill histograms from an event
        // Example:
        // for (auto & histMap : fSelectionHists) {
        //     // Fill histograms based on event data
        // }
    }
};

} // namespace protoana

#endif // PROTOANA_ABSCEXSTRATEGY_H
