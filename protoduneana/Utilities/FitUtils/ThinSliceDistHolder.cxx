#include "ThinSliceDistHolder.h"

void protoana::ThinSliceDistHolder::Reset() {
    for (auto & hist_map : fSelectionHists) {
        for (auto & [key, vec] : hist_map) {
            for (auto & hist : vec) {
                hist->Reset();
            }
        }
    }

    for (auto & hist_map : fInteractionHists) {
        for (auto & [key, val] : hist_map) {
            val->Reset();
        }
    }
    for (auto & hist_map : fIncidentHists) {
        for (auto & [key, val] : hist_map) {
            val->Reset();
        }
    }
    for (auto & hist_map : fXSecHists) {
        for (auto & [key, val] : hist_map) {
            val->Reset();
        }
    }

}