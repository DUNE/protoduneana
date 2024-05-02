#ifndef THINSLICEEVENT_hh
#define THINSLICEEVENT_hh
#include "TSpline.h"
namespace protoana {
class ThinSliceEvent {
 public:
  ThinSliceEvent(int event, int subrun, int run)
    : event_ID(event), subrun_ID(subrun), run_ID(run) {
    sample_ID = -999;
    selection_ID = -999;
    true_beam_interactingEnergy = -999;
    true_beam_initEnergy = -999.;
    reco_beam_interactingEnergy = -999;
    true_beam_endP = -999;
    true_beam_startP = -999;
    true_beam_endZ = -999.;
    true_beam_mass = -999;
    reco_beam_endZ = -999;
    reco_beam_startY = -999.;
    reco_beam_startX_SCE = -999.;
    reco_beam_startY_SCE = -999.;
    reco_beam_startZ_SCE = -999.;
    beam_inst_P = -999;
    pdg = -999;
    is_beam_scraper = false;
    reco_beam_incidentEnergies = std::vector<double>();
    true_beam_incidentEnergies = std::vector<double>();
    true_beam_traj_Z = std::vector<double>();
    true_beam_traj_KE = std::vector<double>();
    true_beam_slices = std::vector<int>();
    calibrated_dQdX = std::vector<double>();
    beam_EField = std::vector<double>();
    track_pitch = std::vector<double>();
    g4rw_weights = std::map<std::string, std::vector<double>>();
    g4rw_splines = std::map<std::string, TSpline3*>();
    reco_daughter_track_thetas = std::vector<double>();
    reco_daughter_track_scores = std::vector<double>();
    reco_daughter_track_dQdX = std::vector<std::vector<double>>();
    reco_daughter_track_res_range = std::vector<std::vector<double>>();
    reco_daughter_efield = std::vector<std::vector<double>>();
    has_pi0_shower = false;
    true_daughter_PDGs = std::vector<int>();
    true_n_neutrons = -999;
    true_n_protons = -999;
    reco_beam_origin = -999;
    reco_daughter_truncated_dEdX = std::vector<double>();
    reco_daughter_chi2s_perhit = std::vector<double>();
    stored_reco_energy = -999.;

  };
  /*
  ~ThinSliceEvent() {
    for (auto s : g4rw_splines) {
      delete s.second;
    }
  };*/

  /*
  int GetEventID() const {
    return event_ID;
  };

  int GetSubrunID() const {
    return subrun_ID;
  };

  int GetRunID() const {
    return run_ID;
  };*/

  int GetSampleID() const {
    return sample_ID;
  };
  void SetSampleID(int s) {
    sample_ID = s;
  };

  int GetSelectionID() const {
    return selection_ID;
  };
  void SetSelectionID(int s) {
    selection_ID = s;
  };

  int GetCalUpSelectionID() const {
    return cal_up_selection_ID;
  };
  void SetCalUpSelectionID(int s) {
    cal_up_selection_ID = s;
  };

  int GetCalDownSelectionID() const {
    return cal_down_selection_ID;
  };
  void SetCalDownSelectionID(int s) {
    cal_down_selection_ID = s;
  };


  bool GetHasPi0Shower() const {
    return has_pi0_shower;
  };
  void SetHasPi0Shower(bool s) {
    has_pi0_shower = s;
  };

  double GetTrueInteractingEnergy() const {
    return true_beam_interactingEnergy;
  };
  void SetTrueInteractingEnergy(double e) {
    true_beam_interactingEnergy = e;
  };

  double GetTrueInitEnergy() const {
    return true_beam_initEnergy;
  };
  void SetTrueInitEnergy(double e) {
    true_beam_initEnergy = e;
  };

  double GetRecoInteractingEnergy() const {
    return reco_beam_interactingEnergy;
  };
  void SetRecoInteractingEnergy(double e) {
    reco_beam_interactingEnergy = e;
  };

  double GetRecoStartY() const {
    return reco_beam_startY;
  };
  void SetRecoStartY(double y) {
    reco_beam_startY = y;
  };

  double GetRecoStartX_SCE() const {
    return reco_beam_startX_SCE;
  };
  void SetRecoStartX_SCE(double x) {
    reco_beam_startX_SCE = x;
  };
  double GetRecoStartY_SCE() const {
    return reco_beam_startY_SCE;
  };
  void SetRecoStartY_SCE(double y) {
    reco_beam_startY_SCE = y;
  };
  double GetRecoStartZ_SCE() const {
    return reco_beam_startZ_SCE;
  };
  void SetRecoStartZ_SCE(double z) {
    reco_beam_startZ_SCE = z;
  };

  double GetTrueEndP() const {
    return true_beam_endP;
  };
  void SetTrueEndP(double p) {
    true_beam_endP = p;
  };

  double GetTrueEndZ() const {
    return true_beam_endZ;
  };
  void SetTrueEndZ(double z) {
    true_beam_endZ = z;
  };

  double GetRecoEndZ() const {
    return reco_beam_endZ;
  };
  void SetRecoEndZ(double p) {
    reco_beam_endZ = p;
  };

  double GetTrueStartP() const {
    return true_beam_startP;
  };
  void SetTrueStartP(double p) {
    true_beam_startP = p;
  };

  double GetTrueMass() const {
    return true_beam_mass;
  };
  void SetTrueMass(double m) {
    true_beam_mass = m;
  };

  double GetVertexMichelScore() const {
    return vertex_michel_score;
  };
  void SetVertexMichelScore(double v) {
    vertex_michel_score = v;
  };
  double GetVertexNHits() const {
    return vertex_nhits;
  };
  void SetVertexNHits(int v) {
    vertex_nhits = v;
  };

  const std::vector<double> & GetRecoIncidentEnergies() const {
    return reco_beam_incidentEnergies;
  };
  void SetRecoIncidentEnergies(std::vector<double> v) {
    reco_beam_incidentEnergies = v;
  };

  const std::vector<double> & GetTrueIncidentEnergies() const {
    return true_beam_incidentEnergies;
  };
  void SetTrueIncidentEnergies(std::vector<double> v) {
    true_beam_incidentEnergies = v;
  };

  const std::vector<double> & GetTrueTrajZ() const {
    return true_beam_traj_Z;
  };
  void SetTrueTrajZ(std::vector<double> v) {
    true_beam_traj_Z = v;
  };

  const std::vector<double> & GetTrueTrajKE() const {
    return true_beam_traj_KE;
  };
  void SetTrueTrajKE(std::vector<double> v) {
    true_beam_traj_KE = v;
  };

  const std::vector<double> & GetRecoDaughterTrackThetas() const {
    return reco_daughter_track_thetas;
  };
  void SetRecoDaughterTrackThetas(std::vector<double> v) {
    reco_daughter_track_thetas = v;
  };

  const std::vector<double> & GetRecoDaughterTrackScores() const {
    return reco_daughter_track_scores;
  };
  void SetRecoDaughterTrackScores(std::vector<double> v) {
    reco_daughter_track_scores = v;
  };

  const std::vector<std::vector<double>>
      & GetRecoDaughterTrackResRanges() const {
    return reco_daughter_track_res_range;
  };
  void AddRecoDaughterTrackResRange(std::vector<double> v) {
    reco_daughter_track_res_range.push_back(v);
  };

  const std::vector<std::vector<double>>
      & GetRecoDaughterTrackdQdXs() const {
    return reco_daughter_track_dQdX;
  };
  void AddRecoDaughterTrackdQdX(std::vector<double> v) {
    reco_daughter_track_dQdX.push_back(v);
  };

  const std::vector<std::vector<double>>
      & GetRecoDaughterEFields() const {
    return reco_daughter_efield;
  };
  void AddRecoDaughterEField(std::vector<double> v) {
    reco_daughter_efield.push_back(v);
  };

  const std::vector<int> & GetTrueSlices() const {
    return true_beam_slices;
  };
  void SetTrueSlices(std::vector<int> v) {
    true_beam_slices = v;
  };

  const std::vector<int> & GetTrueDaughterPDGs() const {
    return true_daughter_PDGs;
  };
  void SetTrueDaughterPDGs(std::vector<int> v) {
    true_daughter_PDGs = v;
  };

  const std::vector<double> & GetTrueDaughterStartPs() const {
    return true_daughter_startPs;
  };
  void SetTrueDaughterStartPs(std::vector<double> v) {
    true_daughter_startPs = v;
  };

  const std::vector<double> & GetdQdXCalibrated() const {
    return calibrated_dQdX;
  };
  void SetdQdXCalibrated(std::vector<double> v) {
    calibrated_dQdX = v;
  };

  const std::vector<double> & GetEField() const {
    return beam_EField;
  };
  void SetEField(std::vector<double> v) {
    beam_EField = v;
  };

  const std::vector<double> & GetTrackPitch() const {
    return track_pitch;
  };
  void SetTrackPitch(std::vector<double> v) {
    track_pitch = v;
  };

  void SetBeamInstP(double p) {
    beam_inst_P = p;
  };
  double GetBeamInstP() const {
    return beam_inst_P;
  };

  void SetPDG(int p) {
    pdg = p;
  };
  int GetPDG() const {
    return pdg;
  };

  void SetStoredRecoEnergy(double e) {stored_reco_energy = e;};
  double GetStoredRecoEnergy() {return stored_reco_energy;};

  void MakeG4RWBranch(const std::string & br, const std::vector<double> & ws) {
    g4rw_weights[br] = ws;
  };
  void MakeG4RWCoeff(const std::string & br, const std::vector<double> & cs) {
    g4rw_coeffs[br] = cs;

    //MakeG4RWExtendCoeffs(br, cs);
  };

  /*
  void MakeG4RWExtendCoeffs(const std::string & br, const std::vector<double> & cs) {
    g4rw_lower_coeffs[br] = GetExpCoeffs(.1, cs);

    double wp2 = GetPolPrime2(2, coeffs);
    double wp = GetPolPrime(2, coeffs);

    if (wp2 < 0. && wp > 0.) {
      double b = -1.*wp2/wp;
      double a = GetPol(2., coeffs)/(1. - exp(-1.*b*2.));
      g4rw_upper_coeffs[br] = {a, b};
      //return a*(1. - exp(-1.*b*x));
    }
    else {
      g4rw_upper_coeffs[br] = GetExpCoeffs(2., cs);
    }
  };*/

  //For extending
  double GetPol(double x, const std::vector<double> & coeffs) const {
    double result = 0.; 
    for (size_t i = 0; i < coeffs.size(); ++i) {
      result += coeffs[i]*std::pow(x, i);
    }
    return result;
  };
  double GetPolPrime(double x, const std::vector<double> & coeffs) const {
    double result = 0.; 
    for (size_t i = 1; i < coeffs.size(); ++i) {
      result += i*coeffs[i]*std::pow(x, i-1);
    }
    return result;
  };
  double GetPolPrime2(double x, const std::vector<double> & coeffs) const {
    double result = 0.; 
    for (size_t i = 2; i < coeffs.size(); ++i) {
      result += i*(i-1)*coeffs[i]*std::pow(x, i-2);
    }
    return result;
  };

  std::pair<double, double> GetExpCoeffs(double x, const std::vector<double> & coeffs) const {
    //std::cout << "Exp coeffs " << x << "\n";
    //std::cout << GetPolPrime(x, coeffs) << " " << GetPol(x, coeffs) << "\n";
    double b = GetPolPrime(x, coeffs)/GetPol(x, coeffs); 
    //std::cout << exp(x*b) << "\n";
    double a = GetPol(x, coeffs)/exp(x*b);
    //std::cout << a << " " << b << std::endl;

    return {a, b};
  };

  double GetLowerWeight(double x, const std::vector<double> & coeffs, double test_pt=.1) const {
    auto exp_coeffs = GetExpCoeffs(test_pt, coeffs);
    return exp_coeffs.first*exp(exp_coeffs.second*x);
  };

  double GetUpperWeight(double x, const std::vector<double> & coeffs) const {
    double wp2 = GetPolPrime2(2, coeffs);
    double wp = GetPolPrime(2, coeffs);

    if (wp2 < 0. && wp > 0.) {
      double b = -1.*wp2/wp;
      double a = GetPol(2., coeffs)/(1. - exp(-1.*b*2.));
      return a*(1. - exp(-1.*b*x));
    }
    else {
      return GetLowerWeight(x, coeffs, 2.);
    }
  }

  //Get Weight from the polynomial defined by the coeffs
  double GetG4RWCoeffWeight(const std::string & br, double input, bool extend=false) const {
    auto & coeffs = g4rw_coeffs.at(br);
    if (coeffs.size() == 0) {
      return 1.;
    }

    //std::cout << GetLowerWeight(.1, coeffs) << " " << GetPol(.1, coeffs) << std::endl;
    //std::cout << GetUpperWeight(2, coeffs) << " " << GetPol(2, coeffs) << std::endl;

    if (extend && input < .1) {
      double weight = GetLowerWeight(input, coeffs);
      //std::cout << "Extended low " << std::endl;
      //std::cout << weight << std::endl;
      return weight;
    }
    else if (extend && input > 2.) {
      //std::cout << "Extended high " << GetUpperWeight(input, coeffs) << std::endl;
      return GetUpperWeight(input, coeffs);
    }

    /*double results = 0.;
    for (size_t i = 0; i < coeffs.size(); ++i) {
      results += coeffs[i]*std::pow(input, i); 
    }
    std::cout << results << " " << GetPol(input, coeffs) << std::endl;
    return results;*/
    return GetPol(input, coeffs);
  };

  double GetG4RWWeight(const std::string & br, size_t i) const {
    if (g4rw_weights.at(br).size() == 0) return 1.;
    return g4rw_weights.at(br).at(i); 
  };
  const std::map<std::string, std::vector<double>> & GetG4RWWeightMap() const {
    return g4rw_weights; 
  };
  const std::vector<double> & GetG4RWBranch(const std::string & br) const {
    return g4rw_weights.at(br);
  };
  bool HasG4RWBranch(const std::string & br) const {
    return (g4rw_weights.find(br) != g4rw_weights.end());
  };

  void MakeG4RWSpline(const std::string & br) {
    std::vector<double> vars;
    if (!g4rw_weights[br].size()) return;
    for (size_t i = 0; i < g4rw_weights[br].size(); ++i) {
      vars.push_back(.1*(1+i));
    }
    std::string name = "g4rw_spline_event_" + std::to_string(event_ID) + "_" +
                       std::to_string(subrun_ID);
    g4rw_splines[br] = new TSpline3(name.c_str(), &vars[0],
                                    &g4rw_weights[br][0], vars.size());
  };

  int GetEventID() const {return event_ID;};
  int GetSubrunID() const {return subrun_ID;};
  int GetRunID() const {return run_ID;};

  int GetTrueID() const {return true_beam_ID;};
  int GetRecoToTrueID() const {return reco_beam_true_byHits_ID;};

  void SetTrueID(int id) {true_beam_ID = id;};
  void SetRecoToTrueID(int id) {reco_beam_true_byHits_ID = id;};

  double GetDeltaEToTPC() const {return delta_e_to_tpc;};
  void SetDeltaEToTPC(double delta_e) {delta_e_to_tpc = delta_e;};

  double GetLeadingPCostheta() const {return leading_p_costheta;};
  double GetLeadingPiPlusCostheta() const {return leading_piplus_costheta;};
  double GetLeadingPi0Costheta() const {return leading_pi0_costheta;};

  double GetLeadingPMomentum() const {return leading_p_momentum;};
  double GetLeadingPiPlusMomentum() const {return leading_piplus_momentum;};
  double GetLeadingPi0Momentum() const {return leading_pi0_momentum;};

  void SetLeadingPCostheta(double p) {leading_p_costheta = p;};
  void SetLeadingPiPlusCostheta(double p) {leading_piplus_costheta = p;};
  void SetLeadingPi0Costheta(double p) {leading_pi0_costheta = p;};
  void SetLeadingPMomentum(double p) {leading_p_momentum = p;};
  void SetLeadingPiPlusMomentum(double p) {leading_piplus_momentum = p;};
  void SetLeadingPi0Momentum(double p) {leading_pi0_momentum = p;};

  void SetRecoOrigin(int origin) {reco_beam_origin = origin;};
  int GetRecoOrigin() const {return reco_beam_origin;};

  void SetIsBeamScraper(bool val) {is_beam_scraper = val;};
  bool GetIsBeamScraper() const {return is_beam_scraper;};

  void AddOneChi2PerHit(double val) {
    reco_daughter_chi2s_perhit.push_back(val);
  };
  void SetChi2PerHit(std::vector<double> & vals) {reco_daughter_chi2s_perhit = vals;};
  const std::vector<double> & GetChi2sPerHit() const {return reco_daughter_chi2s_perhit;};

  void AddOneTrunc_dEdX(double val) {
    reco_daughter_truncated_dEdX.push_back(val);
  };
  void SetTrunc_dEdXs(std::vector<double> & vals) {reco_daughter_truncated_dEdX = vals;};
  const std::vector<double> & GetTrunc_dEdXs() const {return reco_daughter_truncated_dEdX;};

  void SetMCStatVarWeight(double w) {mc_stat_var_weight = w;};
  double GetMCStatVarWeight() const {return mc_stat_var_weight;};

  void SetTrueNNeutrons(int i) {true_n_neutrons = i;};
  double GetTrueNNeutrons() const {return true_n_neutrons;};

  void SetTrueNProtons(int i) {true_n_protons = i;};
  double GetTrueNProtons() const {return true_n_protons;};

  void SetTrueNPiPlus(int i) {true_n_piplus = i;};
  double GetTrueNPiPlus() const {return true_n_piplus;};

  void SetTrueNPiMinus(int i) {true_n_piminus = i;};
  double GetTrueNPiMinus() const {return true_n_piminus;};

  void SetTrueNPi0(int i) {true_n_pi0 = i;};
  double GetTrueNPi0() const {return true_n_pi0;};



 private:
  int event_ID, subrun_ID, run_ID;
  int sample_ID;
  int selection_ID, cal_up_selection_ID, cal_down_selection_ID;
  int pdg;
  double true_beam_interactingEnergy, reco_beam_interactingEnergy;
  double true_beam_initEnergy;
  double true_beam_endP, true_beam_mass;
  double reco_beam_endZ, true_beam_startP, true_beam_endZ;
  double reco_beam_startY;
  double reco_beam_startX_SCE, reco_beam_startY_SCE, reco_beam_startZ_SCE;
  double beam_inst_P;
  bool has_pi0_shower;
  std::vector<double> reco_beam_incidentEnergies,
                      true_beam_incidentEnergies,
                      true_beam_traj_Z,
                      true_beam_traj_KE,
                      reco_daughter_track_thetas,
                      reco_daughter_track_scores;
  std::vector<std::vector<double>> reco_daughter_track_dQdX,
                                   reco_daughter_track_res_range,
                                   reco_daughter_efield;

  std::vector<int> true_beam_slices, true_daughter_PDGs;
  std::vector<double> true_daughter_startPs;
  std::vector<double> calibrated_dQdX, beam_EField,
                      track_pitch;
  std::map<std::string, std::vector<double>> g4rw_weights;
  std::map<std::string, std::vector<double>> g4rw_coeffs;
  std::map<std::string, std::pair<double, double>> g4rw_upper_coeffs,
                                                   g4rw_lower_coeffs;
  std::map<std::string, TSpline3*> g4rw_splines;
  int true_beam_ID;
  int reco_beam_true_byHits_ID;
  double delta_e_to_tpc;
  double leading_p_costheta, leading_piplus_costheta, leading_pi0_costheta;
  double leading_p_momentum, leading_piplus_momentum, leading_pi0_momentum;
  int reco_beam_origin = -999;
  bool is_beam_scraper;

  std::vector<double> reco_daughter_truncated_dEdX,
                      reco_daughter_chi2s_perhit;
  double vertex_michel_score;
  int vertex_nhits;
  double stored_reco_energy = -999.;
  double mc_stat_var_weight = 1.;
  int true_n_neutrons = -999, true_n_protons = -999, true_n_piplus = -999,
      true_n_piminus = -999, true_n_pi0 = -999;
};
}
#endif
