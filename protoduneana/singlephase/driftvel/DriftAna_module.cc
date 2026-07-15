////////////////////////////////////////////////////////////////////////
// Class:       DriftAna
// Plugin Type: analyzer (art v3_03_01)
// File:        DriftAna_module.cc
//
// Generated at Sun Dec  1 22:31:29 2019 by Tingjun Yang using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "TH1D.h"
#include "TGraph.h"

class DriftAna;

class DriftAna : public art::EDAnalyzer {
public:
  explicit DriftAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DriftAna(DriftAna const&) = delete;
  DriftAna(DriftAna&&) = delete;
  DriftAna& operator=(DriftAna const&) = delete;
  DriftAna& operator=(DriftAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  TH1D *drifttimenosce;
  TH1D *drifttimeave;
  TH1D *drifttimeint;
  TH1D *fEfield;

  TH1D *fdeltatnosce;
  TH1D *fdeltatave;
  TH1D *fdeltatint;

  TH1D *fdriftvnosce;
  TH1D *fdriftvave;
  TH1D *fdriftvint;

  TH1D *fefieldnosce;
  TH1D *fefieldave;
  TH1D *fefieldint;
};


DriftAna::DriftAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void DriftAna::analyze(art::Event const& e)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
  double driftv = detProp.DriftVelocity(detProp.PerPlaneEfield(),
                                        detProp.Temperature());
  double efield = detProp.PerPlaneEfield();
  std::cout<<"E field = "<<efield<<" Nominal drift velocity = "<<driftv<<std::endl;

  const int n = 100;
  const double e0 = 0.4;
  const double e1 = 0.7;
  std::vector<double> vecv;
  std::vector<double> vece;
  for (int i = 0; i<n; ++i){
    double e_field = e0 + i*(e1-e0)/n;
    double drift_v = detProp.DriftVelocity(e_field,
                                           detProp.Temperature());
    vecv.push_back(drift_v);
    vece.push_back(e_field);
  }
  TGraph *gr = new TGraph(n, &vecv[0], &vece[0]);
  for (int i = 1; i<=drifttimeave->GetNbinsX(); ++i){
    double xyz[3] = {0};
    xyz[0] = drifttimeave->GetBinCenter(i);
    xyz[1] = 300;
    xyz[2] = 350;
    std::cout<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;
    double driftdis = 358.6 - std::abs(xyz[0]);
    double drifttime = driftdis/driftv;
    drifttimenosce->SetBinContent(i, drifttime);
    geo::Vector_t posOffsets;
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    if (SCE->EnableSimSpatialSCE() == true){
      posOffsets = SCE->GetPosOffsets({ xyz[0], xyz[1], xyz[2] });
    }
    posOffsets.SetX(-posOffsets.X());
    driftdis += posOffsets.X();
    drifttime = driftdis/driftv;
    drifttimeave->SetBinContent(i, drifttime);
    
    double dir = 1;
    if (xyz[0]<0) dir = -1;
    double xyz1[3];
    xyz1[0] = xyz[0];
    xyz1[1] = 300;
    xyz1[2] = 350;
    double deltax = 1;
    double dt = 0;
    while (std::abs(xyz1[0])<358.6){
      geo::Vector_t fEfieldOffsets = SCE->GetEfieldOffsets({ xyz1[0], xyz1[1], xyz1[2] });
      double EField = std::sqrt( (efield + efield*fEfieldOffsets.X())*(efield + efield*fEfieldOffsets.X()) +
                                 (efield*fEfieldOffsets.Y()*efield*fEfieldOffsets.Y()) +
                                 (efield*fEfieldOffsets.Z()*efield*fEfieldOffsets.Z()) );
      double v = detProp.DriftVelocity(EField,
                                       detProp.Temperature());
      dt += deltax/v;
      xyz1[0] += deltax*dir;
    }
    drifttimeint->SetBinContent(i, dt);

    geo::Vector_t fEfieldOffsets = SCE->GetEfieldOffsets({ xyz[0], xyz[1], xyz[2] });
    double EField = std::sqrt( (efield + efield*fEfieldOffsets.X())*(efield + efield*fEfieldOffsets.X()) +
                               (efield*fEfieldOffsets.Y()*efield*fEfieldOffsets.Y()) +
                               (efield*fEfieldOffsets.Z()*efield*fEfieldOffsets.Z()) );
    fEfield->SetBinContent(i, EField);
  }

  for (int i = 1; i<=drifttimeave->GetNbinsX() - 1; ++i){
    double x = drifttimeave->GetBinCenter(i);
    if (x<0){
      fdeltatnosce->SetBinContent(i, drifttimenosce->GetBinContent(i) - drifttimenosce->GetBinContent(i-1));
      fdeltatave->SetBinContent(i, drifttimeave->GetBinContent(i) - drifttimeave->GetBinContent(i-1));
      fdeltatint->SetBinContent(i, drifttimeint->GetBinContent(i) - drifttimeint->GetBinContent(i-1));
    }
    else{
      fdeltatnosce->SetBinContent(i, drifttimenosce->GetBinContent(i) - drifttimenosce->GetBinContent(i+1));
      fdeltatave->SetBinContent(i, drifttimeave->GetBinContent(i) - drifttimeave->GetBinContent(i+1));
      fdeltatint->SetBinContent(i, drifttimeint->GetBinContent(i) - drifttimeint->GetBinContent(i+1));
    }
    fdriftvnosce->SetBinContent(i, fdeltatnosce->GetBinWidth(i)/fdeltatnosce->GetBinContent(i));
    fdriftvave->SetBinContent(i, fdeltatave->GetBinWidth(i)/fdeltatave->GetBinContent(i));
    fdriftvint->SetBinContent(i, fdeltatint->GetBinWidth(i)/fdeltatint->GetBinContent(i));
    fefieldnosce->SetBinContent(i, gr->Eval(fdeltatnosce->GetBinWidth(i)/fdeltatnosce->GetBinContent(i)));
    fefieldave->SetBinContent(i, gr->Eval(fdeltatave->GetBinWidth(i)/fdeltatave->GetBinContent(i)));
    fefieldint->SetBinContent(i, gr->Eval(fdeltatint->GetBinWidth(i)/fdeltatint->GetBinContent(i)));
  }
}

void DriftAna::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;
    drifttimenosce = tfs->make<TH1D>("drifttimenosce",";x (cm);Drift time (#mus)",120,-360,360);
    drifttimeave = tfs->make<TH1D>("drifttimeave",";x (cm);Drift time (#mus)",120,-360,360);
    drifttimeint = tfs->make<TH1D>("drifttimeint",";x (cm);Drift time (#mus)",120,-360,360);
    fEfield = tfs->make<TH1D>("fEfield","E field; x (cm);E field (kV/cm)",120,-360,360);

    fdeltatnosce = tfs->make<TH1D>("fdeltatnosce",";x (cm);#Delta t (#mus)", 120, -360, 360);
    fdeltatave = tfs->make<TH1D>("fdeltatave",";x (cm);#Delta t (#mus)", 120, -360, 360);
    fdeltatint = tfs->make<TH1D>("fdeltatint",";x (cm);#Delta t (#mus)", 120, -360, 360);

    fdriftvnosce = tfs->make<TH1D>("fdriftvnosce",";x (cm);Drift velocity (cm/#mus)", 120, -360, 360);
    fdriftvave = tfs->make<TH1D>("fdriftvave",";x (cm);Drift velocity (cm/#mus)", 120, -360, 360);
    fdriftvint = tfs->make<TH1D>("fdriftvint",";x (cm);Drift velocity (cm/#mus)", 120, -360, 360);

    fefieldnosce = tfs->make<TH1D>("fefieldnosce",";x (cm);E field (kV/cm)", 120, -360, 360);
    fefieldave = tfs->make<TH1D>("fefieldave",";x (cm);E field (kV/cm)", 120, -360, 360);
    fefieldint = tfs->make<TH1D>("fefieldint",";x (cm);E field (kV/cm)", 120, -360, 360);

}

DEFINE_ART_MODULE(DriftAna)
