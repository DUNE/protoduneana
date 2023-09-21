#ifndef KaonExcReweighter_h
#define KaonExcReweighter_h

#include "geant4reweight/ReweightBase/G4Reweighter.hh"

class G4KaonPlus;

class KaonExcReweighter : public G4Reweighter {
  public:

    KaonExcReweighter(TFile *, const std::map<std::string, TH1D*> &,
                       const fhicl::ParameterSet &,
                       G4ReweightManager * rw_manager,
                       TH1D * inputElasticBiasHist = nullptr, bool fix = false);
    virtual ~KaonExcReweighter();
    std::string GetInteractionSubtype(const G4ReweightTraj &) override;

  protected:
    G4KaonPlus * kplus;
};

#endif
