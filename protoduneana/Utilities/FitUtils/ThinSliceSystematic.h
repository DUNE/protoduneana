#ifndef THINSLICESYSTEMATIC_hh
#define THINSLICESYSTEMATIC_hh
#include "fhiclcpp/ParameterSet.h"
namespace protoana {
class ThinSliceSystematic {
 public:
  ThinSliceSystematic(const fhicl::ParameterSet & pset,
                      int selection_id=-1, int selection_bin=-1)
    : fType(pset.get<std::string>("Name")),
      fVal(pset.get<double>("Central")),
      fCentral(pset.get<double>("Central")),
      fUpperLimit(pset.get<double>("UpperLimit")),
      fLowerLimit(pset.get<double>("LowerLimit")),
      fSetLimits(pset.get<bool>("SetLimits", true)),
      fThrowLimit(pset.get<double>("ThrowLimit")),
      fThrowLimitUp(pset.get<double>("ThrowLimitUp")),
      fGenThrowLimit(pset.get<double>("GenThrowLimit")),
      fGenThrowLimitUp(pset.get<double>("GenThrowLimitUp")),
      fSigma(pset.get<double>("Sigma", 1.)),
      fIsG4RWCoeff(pset.get<bool>("IsG4RWCoeff", false)),
      fIsTiedG4RWCoeff(pset.get<bool>("IsTiedG4RWCoeff", false)),
      fIsSelVar(pset.get<bool>("IsSelVar", false)),
      fIsCosThetaEff(pset.get<bool>("IsCosThetaEff", false)),
      fSelectionID(selection_id),
      fSelectionBin(selection_bin),
      fOptions(pset.get<fhicl::ParameterSet>("Options")) {
    fName = "par_" + fType + "_syst";
    if (fIsG4RWCoeff) {
      fG4RWCoeffBranch = fOptions.get<std::string>("Branch");
      fG4RWExtend = fOptions.get<bool>("Extend", false);
    }
    if (fIsTiedG4RWCoeff) {
      fTiedG4RWCoeffBranches
          = fOptions.get<std::vector<std::string>>("Branches");
    }

    if (fIsSelVar) {
      fName = "par_" + fType + "_syst_" + std::to_string(fSelectionID) + "_" +
              std::to_string(fSelectionBin);
      std::cout << "SelectionVar: " << fName << std::endl;
    }
  };

  ThinSliceSystematic(){};

  ~ThinSliceSystematic(){};

  const double GetValue() const {
    return fVal;
  };

  const std::string & GetG4RWCoeffBranch() const {
    return fG4RWCoeffBranch;
  };

  const bool & GetG4RWExtend() const {
    return fG4RWExtend;
  };

  const std::vector<std::string> & GetTiedG4RWBranches() const {
    return fTiedG4RWCoeffBranches;
  };

  const double GetCentral() const {
    return fCentral;
  };
  void SetCentral(double val) {
    fCentral = val;
  };

  const double GetUpperLimit() const {
    return fUpperLimit;
  };

  const bool GetSetLimits() const {
    return fSetLimits;
  };

  const double GetThrowLimit() const {
    return fThrowLimit;
  };
  const double GetThrowLimitUp() const {
    return fThrowLimitUp;
  };

  const double GetGenThrowLimit() const {
    return fGenThrowLimit;
  };
  const double GetGenThrowLimitUp() const {
    return fGenThrowLimitUp;
  };

  const double GetLowerLimit() const {
    return fLowerLimit;
  };

  const std::string GetName() const {
    return fName;
  };

  void SetValue(double v) {
    fVal = v;
  };

  void SetToNSigma(int n = 1) {
    SetValue(n*fSigma);
  }

  void SetToCentral() {
    SetValue(fCentral); 
  }

  template <typename T> const T GetOption(std::string name) const {
    return fOptions.get<T>(name);
  }

  const bool GetIsG4RWCoeff() const {
    return fIsG4RWCoeff;
  }

  const bool GetIsCosThetaEff() const {
    return fIsCosThetaEff;
  }

  const bool GetIsTiedG4RWCoeff() const {
    return fIsTiedG4RWCoeff;
  }

  const bool GetIsSelVar() const {
    return fIsSelVar;
  }

  const int GetSelectionID() const {
    return fSelectionID;
  }

  const int GetSelectionBin() const {
    return fSelectionBin;
  }

 private:

  std::string fType;
  double fVal;
  double fCentral;
  double fUpperLimit; 
  double fLowerLimit;
  bool fSetLimits;
  double fThrowLimit; 
  double fThrowLimitUp; 
  double fGenThrowLimit; 
  double fGenThrowLimitUp; 
  double fSigma;

  bool fIsG4RWCoeff, fIsTiedG4RWCoeff, fIsSelVar, fIsCosThetaEff;
  int fSelectionID, fSelectionBin;
  std::string fG4RWCoeffBranch = "";
  bool fG4RWExtend = false;
  std::vector<std::string> fTiedG4RWCoeffBranches;

  fhicl::ParameterSet fOptions;

  std::string fName;
};
}
#endif
