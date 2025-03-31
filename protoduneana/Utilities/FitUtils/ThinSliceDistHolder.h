#ifndef THINSLICEDISTHOLDER_hh
#define THINSLICEDISTHOLDER_hh

#include "TH1.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include <map>

namespace protoana {
  using TrueCat_t = int;
  using SelID_t = int;
  using TrueCatHist_map = std::map<TrueCat_t, std::shared_ptr<TH1>>;
  using SelIDHist_map = std::map<SelID_t, std::shared_ptr<TH1>>;
  using TrueCatSelID_t = std::pair<TrueCat_t, SelID_t>;
  using TrueCatSelID_map = std::map<TrueCatSelID_t, std::vector<std::shared_ptr<TH1>>>;
  using Range_t = std::pair<double, double>;
  using Range_vec = std::vector<Range_t>;
  using Range_map = std::map<TrueCat_t, Range_vec>;
class ThinSliceDistHolder {


public:
  ThinSliceDistHolder() = default;
  ~ThinSliceDistHolder() = default;
  friend class ThinSliceDistBuilder; //This will provide access to the builder
  friend class ThinSliceDistBuilder1D; //This will provide access to the builder
  void SetIDs(const fhicl::ParameterSet & pset);
  void Reset();

  void ScaleInBeamBin(int beam_bin, double value) {
    for (auto & [true_id, hist] : fInteractionHists[beam_bin]) {
      hist->Scale(value);
    }
    for (auto & [true_id, hist] : fIncidentHists[beam_bin]) {
      hist->Scale(value);
    }


    //Extra dimension here
    for (auto & [ids, hists] : fSelectionHists[beam_bin]) {
      for (auto & hist : hists) {
        hist->Scale(value);
      }
    }

  }

  void CalcTotalDists();


  const std::vector<TrueCatSelID_map> & GetSelectionHists() const { return fSelectionHists; }
  const SelIDHist_map & GetTotalSelectionHists() const { return fTotalSelectionHists; }

  const std::vector<TrueCatHist_map> &  GetInteractionHists() const { return fInteractionHists; }
  const std::vector<TrueCatHist_map> &  GetIncidentHists() const { return fIncidentHists; }
  const TrueCatHist_map &  GetXSecHists() const { return fXSecHists; }
  const TrueCatHist_map &  GetTotalIncidentHists() const { return fTotalIncidentHists; }


  std::shared_ptr<TH1> GetInteractionHist(TrueCat_t true_id, size_t beam_bin) {
    if (beam_bin >= fInteractionHists.size()) {
      throw cet::exception("ThinSliceDistHolder")
        << "Requested out of range beam bin "
        << beam_bin << " only have: " << fInteractionHists.size() << "\n"; 
    }


    if (fInteractionHists.at(beam_bin).find(true_id) == fInteractionHists.at(beam_bin).end()) {
      throw cet::exception("ThinSliceDistHolder")
        << "Requested bad true id " << true_id << "\n"; 
    }
    return fInteractionHists.at(beam_bin).at(true_id);
  }

  std::shared_ptr<TH1> GetIncidentHist(TrueCat_t true_id, size_t beam_bin) {
    if (beam_bin >= fIncidentHists.size()) {
      throw cet::exception("ThinSliceDistHolder")
        << "Requested out of range beam bin "
        << beam_bin << " only have: " << fIncidentHists.size() << "\n"; 
    }


    if (fIncidentHists.at(beam_bin).find(true_id) == fIncidentHists.at(beam_bin).end()) {
      throw cet::exception("ThinSliceDistHolder")
        << "Requested bad true id " << true_id << "\n"; 
    }
    return fIncidentHists.at(beam_bin).at(true_id);
  }

  std::shared_ptr<TH1> GetSelectionHist(TrueCat_t true_id, SelID_t sel_id, size_t beam_bin, size_t signal_bin) {
    if (beam_bin >= fSelectionHists.size()) {
      throw cet::exception("ThinSliceDistHolder")
        << "Requested out of range beam bin "
        << beam_bin << " only have: " << fSelectionHists.size() << "\n"; 
    }

    if (fSelectionHists.at(beam_bin).find({true_id, sel_id}) == fSelectionHists.at(beam_bin).end()) {
      throw cet::exception("ThinSliceDistHolder")
        << "Requested bad selection id "
        << sel_id << " or true id " << true_id << "\n"; 
    }

    // if (std::find(fSelectionHists.at(beam_bin).begin(), fSelectionHists.at(beam_bin).begin(), sel_id) == fSelectionHists.at(beam_bin).begin().end()) {
    //   throw cet::exception("ThinSliceDistHolder")
    //     << "Requested bad selection id"
    //     << sel_id << "\n"; 
    // }

    return fSelectionHists.at(beam_bin).at({true_id, sel_id}).at(signal_bin);
  }

  bool IsInteraction(TrueCat_t true_id) const {
    return (std::find(fMeasurementIDs.begin(), fMeasurementIDs.end(), true_id) != fMeasurementIDs.end());
  }
  bool IsIncident(TrueCat_t true_id) const {
    return (std::find(fIncidentIDs.begin(), fIncidentIDs.end(), true_id) != fIncidentIDs.end());
  }

  const Range_map & GetSignalRanges() const { return fSignalRanges; }
  const std::vector<TrueCat_t> & GetTrueCatIDs() const { return fTrueCatIDs; }
  const std::vector<TrueCat_t> & GetMeasurementIDs() const { return fMeasurementIDs; }
  const std::vector<TrueCat_t> & GetIncidentIDs() const { return fIncidentIDs; }
  const std::vector<TrueCat_t> & GetSignalIDs() const { return fSignalIDs; }
  const std::vector<SelID_t> & GetSelIDs() const { return fSelIDs; }



private:

    //Each of the following is a vector corresponding the beam energy bins
      //Map of True Cat, Sel ID --> Vector of Bins (for signal) --> Hists
      std::vector<TrueCatSelID_map> fSelectionHists;

      //Map of True Cat --> Hist representing interactions, incidents, xsecs
      std::vector<TrueCatHist_map> fInteractionHists;
      std::vector<TrueCatHist_map> fIncidentHists;

    //XSec hists are collapsed over the beam bins, so just a single map
    TrueCatHist_map fXSecHists;
    TrueCatHist_map fTotalIncidentHists;

    //Total selection hists + stacks are collapsed over beam bins and true cats
    SelIDHist_map fTotalSelectionHists;
    // SelIDStack_map fSelectionStacks;

    std::vector<TrueCat_t> fTrueCatIDs;
    std::vector<SelID_t> fSelIDs;
    std::vector<TrueCat_t> fSignalIDs;
    std::vector<TrueCat_t> fMeasurementIDs;
    std::vector<TrueCat_t> fIncidentIDs;
    Range_map fSignalRanges;
};

}
#endif
