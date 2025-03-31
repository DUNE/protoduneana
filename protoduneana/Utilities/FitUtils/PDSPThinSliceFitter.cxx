// C++ language includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include "math.h"
#include "stdio.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"

#include "TStyle.h"
#include "TROOT.h"
#include <ROOT/TSeq.hxx>
#include "TCanvas.h"
#include "TLine.h"
#include "TList.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TDecompChol.h"
#include "TF1.h"
#include "TPad.h"
#include "TVector.h"

#include "PDSPThinSliceFitter.h"
#include "AbsCexDriver.h"

#include "ThinSliceDriverRegistry.h"
#include "ThinSliceStrategyRegistry.h"


protoana::PDSPThinSliceFitter::PDSPThinSliceFitter(
    std::string fcl_file, std::string output_file,
    std::string mc_file, std::string data_file, std::string refit_file,
    std::string tune_file, bool retune, int njobs)
    : fOutputFile(output_file.c_str(), "RECREATE"),
      fMCFileName(mc_file), fDataFileName(data_file),
      fRefitFile(refit_file), fRetune(retune), fNJobs(njobs) {

  Configure(fcl_file);

  if (tune_file.size()) {
    Tune(tune_file);
  }

  gStyle->SetOptStat(0);
  gROOT->SetBatch(1);
  
  try {
    fThinSliceDriver = protoana::ThinSliceDriverRegistry::Instance()->GetDriver(
        fDriverName, fAnalysisOptions);
  }
  catch (const std::runtime_error & e) {
    protoana::ThinSliceDriverRegistry::Instance()->PrintAvailableDrivers();
    throw e;
  }

  InitializeMCSamples();
  fDataSet = ThinSliceDataSet(fSelectionSets, fBeamEnergyBins);
  // SetupExtraHists();
}

protoana::PDSPThinSliceFitter::~PDSPThinSliceFitter() {
  fOutputFile.Close();
}

void protoana::PDSPThinSliceFitter::ResetParameters() {
  size_t n_par = 0;
  for (auto it = fSignalParameters.begin();
       it != fSignalParameters.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (fRandomStart) {
        it->second[i] = fRNG.Gaus(1., .1)*it->second[i];
      }
      fMinimizer->SetVariable(n_par, fSignalParameterNames[it->first][i].c_str(),
                              it->second[i], 0.01);
      if (fSetSigLimits)
        fMinimizer->SetVariableLimits(n_par, fLowerLimit, fUpperLimit);
      ++n_par;
    }
  }

  for (auto it = fFluxParameters.begin();
       it != fFluxParameters.end(); ++it) {
    if (fRandomStart) {
      it->second = fRNG.Gaus(1., .1)*it->second;
    }

    if (fFixFluxes.find(it->first) != fFixFluxes.end()) {
      std::cout << "Fixing flux " << it->first << " " << it->second <<
                   std::endl;
      fMinimizer->FixVariable(n_par);
    }

    fMinimizer->SetVariable(n_par, fFluxParameterNames[it->first].c_str(),
                            it->second, 0.01);
    if (fSetSigLimits) {
      fMinimizer->SetVariableLimits(n_par, fLowerLimit, fUpperLimit);
    }

    ++n_par;
  }

  for (auto it = fSystParameters.begin(); it != fSystParameters.end(); ++it) {
    std::cout << "Adding parameter " << it->second.GetName().c_str() <<
                 " " << it->second.GetCentral() << " " <<
                 it->second.GetLowerLimit() << " " <<
                 it->second.GetUpperLimit() << std::endl;
    double val = (fThrowAndFixSysts ? it->second.GetValue() :
                  it->second.GetCentral());
    if (fRandomStart) {
      it->second.SetValue(fRNG.Gaus(1., .1)*it->second.GetCentral());
      val = it->second.GetValue();
    }
    fMinimizer->SetVariable(n_par, it->second.GetName().c_str(),
                            val, 0.01);
    if (fSetSystLimits && it->second.GetSetLimits()) {
      std::cout << "Setting limit " << it->second.GetLowerLimit() << " " <<
                   it->second.GetUpperLimit() << std::endl;
      fMinimizer->SetVariableLimits(n_par, it->second.GetLowerLimit(),
                                    it->second.GetUpperLimit());
    }
    if (fFixVariables &&
        fSystsToFix.find(it->first) != fSystsToFix.end()) {
      std::cout << "Fixing variable " << it->second.GetName() <<
                   " to value " << fSystsToFix.at(it->first) <<
                   std::endl;
      fMinimizer->SetVariableValue(n_par, fSystsToFix.at(it->first));
      fMinimizer->FixVariable(n_par);
      it->second.SetValue(fSystsToFix.at(it->first));
    }
    else if (fThrowAndFixSysts) {
      //This is set above, just fix it here
      fMinimizer->FixVariable(n_par);
    }
    ++n_par;
  }


  for (auto it = fG4RWParameters.begin(); it != fG4RWParameters.end(); ++it) {
    std::cout << "Adding parameter " << it->second.GetName().c_str() <<
                 " " << it->second.GetCentral() << " " <<
                 it->second.GetLowerLimit() << " " <<
                 it->second.GetUpperLimit() << std::endl;
    double val = it->second.GetCentral();
    if (fRandomStart) {
      it->second.SetValue(fRNG.Gaus(1., .1)*it->second.GetCentral());
      val = it->second.GetValue();
    }
    fMinimizer->SetVariable(n_par, it->second.GetName().c_str(),
                            val, 0.01);
    fMinimizer->SetVariableLimits(n_par, it->second.GetLowerLimit(),
                                  it->second.GetUpperLimit());
 
    if (fTuneG4RWPars) {
      fMinimizer->FixVariable(n_par);
      std::cout << "fixed " << n_par << " to " << val << std::endl;
    }
    ++n_par;
  }
}


//Good
void protoana::PDSPThinSliceFitter::MakeMinimizer() {
  fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
    (ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));

  fMinimizer->SetMaxFunctionCalls(fMaxCalls);
  fMinimizer->SetMaxIterations(fMaxIterations);
  fMinimizer->SetTolerance(fTolerance);
  fMinimizer->SetPrintLevel(fPrintLevel);

  size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters +
                            fTotalSystParameters + fTotalG4RWParameters;
  TH1D parsHist("preFitPars", "", total_parameters, 0, total_parameters);
  fPreFitParsNormal = TH1D("preFitParsNormal", "", 
                           total_parameters, 0, total_parameters);

  bool sig_throw_limits = (fSetSigLimits && !fRemoveSigThrowLimits);
  size_t n_par = 0;
  for (auto it = fSignalParameters.begin();
       it != fSignalParameters.end(); ++it) {
    const auto & fixed_bins = fSignalFixedBins[it->first];
    for (size_t i = 0; i < it->second.size(); ++i) {
      if (fRandomStart) {
        it->second[i] = fRNG.Gaus(1., .1);
      }
      fMinimizer->SetVariable(n_par, fSignalParameterNames[it->first][i].c_str(),
                              it->second[i], 0.01);
      if (fSetSigLimits)
        fMinimizer->SetVariableLimits(n_par, fLowerLimit, fUpperLimit);

      if (std::find(fixed_bins.begin(), fixed_bins.end(), i) != fixed_bins.end())
        fMinimizer->FixVariable(n_par);
      fMinimizerInitVals.push_back(it->second[i]);
      parsHist.SetBinContent(n_par+1, it->second[i]);
      parsHist.SetBinError(n_par+1, 0.);
      //parsHist.GetXaxis()->SetBinLabel(
      //    n_par+1, fSignalParameterNames[it->first][i].c_str());

      fPreFitParsNormal.SetBinContent(n_par+1, it->second[i]);
      fPreFitParsNormal.SetBinError(n_par+1, 0.);
      fPreFitParsNormal.GetXaxis()->SetBinLabel(
          n_par+1, fSignalParameterNames[it->first][i].c_str());
      fParLimits.push_back((sig_throw_limits ? 0. : -1.e6)); //Set to something very low if not setting limits
      fParLimitsUp.push_back(100.);
      ++n_par;

    }
  }

  for (auto it = fFluxParameters.begin();
       it != fFluxParameters.end(); ++it) {
    if (fRandomStart) {
      it->second = fRNG.Gaus(1., .1);
    }

    fMinimizer->SetVariable(n_par, fFluxParameterNames[it->first].c_str(),
                            it->second, 0.01);
    if (fFixFluxes.find(it->first) != fFixFluxes.end()) {
      std::cout << "Fixing flux " << it->first << " " << it->second <<
                   std::endl;
      fMinimizer->FixVariable(n_par);
    }
    fMinimizerInitVals.push_back(it->second);
    if (fSetSigLimits)
      fMinimizer->SetVariableLimits(n_par, fLowerLimit, fUpperLimit);
    parsHist.SetBinContent(n_par+1, it->second);
    parsHist.SetBinError(n_par+1, 0.);
    //parsHist.GetXaxis()->SetBinLabel(
    //    n_par+1, fFluxParameterNames[it->first].c_str());
    fPreFitParsNormal.SetBinContent(n_par+1, it->second);
    fPreFitParsNormal.SetBinError(n_par+1, 0.);
    fPreFitParsNormal.GetXaxis()->SetBinLabel(
        n_par+1, fFluxParameterNames[it->first].c_str());
    //fParLimits.push_back(0.);
    fParLimits.push_back((sig_throw_limits ? 0. : -1.e6));  //Set to something very low if not setting limits
    fParLimitsUp.push_back(100.);
    ++n_par;
  }


  size_t cov_bin_counter = 0;
  for (auto it = fSystParameters.begin(); it != fSystParameters.end(); ++it) {
    std::cout << "Adding parameter " << it->second.GetName().c_str() <<
                 " " << it->second.GetCentral() << " " <<
                 it->second.GetLowerLimit() << " " <<
                 it->second.GetUpperLimit() << std::endl;
    double val = it->second.GetCentral();
    if (fRandomStart) {
      it->second.SetValue(fRNG.Gaus(1., .1)*it->second.GetCentral());
      val = it->second.GetValue();
    }
    fMinimizer->SetVariable(n_par, it->second.GetName().c_str(),
                            val, 0.01);
    fMinimizerInitVals.push_back(val);
    if (fSetSystLimits && it->second.GetSetLimits()) {
      std::cout << "Setting limit " << it->second.GetLowerLimit() << " " <<
                   it->second.GetUpperLimit() << std::endl;
      fMinimizer->SetVariableLimits(n_par, it->second.GetLowerLimit(),
                                    it->second.GetUpperLimit());
    }

        ;

    if (fFixVariables &&
        fSystsToFix.find(it->first) != fSystsToFix.end()) {
      std::cout << "Fixing variable " << it->second.GetName() <<
                   " to value " << fSystsToFix.at(it->first) <<
                   std::endl;
      fMinimizer->SetVariableValue(n_par, fSystsToFix.at(it->first));
      fMinimizer->FixVariable(n_par);
      it->second.SetValue(fSystsToFix.at(it->first));
    }
    //Configured to throw the systematics parameters and fix them 
    //TODO -- EXPAND TO CORRELATED
    else if (fThrowAndFixSysts) {
      //If we're throwing this, generate a random number --
      //elsewise, just use the central
      if (std::find(fSystsToThrowAndFix.begin(),
                    fSystsToThrowAndFix.end(),
                    it->first) != fSystsToThrowAndFix.end()) {
                    //cov_bin_counter) != fSystsToThrowAndFix.end()) {
        //Get the central value & error and throw 1D gaussian
        size_t new_cov_bin = fCovarianceBinsSimple[cov_bin_counter];
        double err = sqrt((*fCovMatrixDisplay)[new_cov_bin][new_cov_bin]);
        /*double */val = fRNG.Gaus(it->second.GetCentral(), err);

        //Check it's in valid range
        //TODO -- RETHROW IF OUT OF RANGE OR JUST SET AT LIMIT?
        if (val < it->second.GetLowerLimit()) val = it->second.GetLowerLimit();
        else if (val > it->second.GetUpperLimit()) val = it->second.GetUpperLimit();
      }
      else {
        val = it->second.GetCentral();
      }

      //Set it and fix
      fMinimizer->SetVariableValue(n_par, val);
      fMinimizer->FixVariable(n_par);
      it->second.SetValue(val);
      std::cout << "Fixing variable (All, Systs): (" << n_par << ", " <<
                   cov_bin_counter << ") to " << val << std::endl;
    }

    fSystParameterIndices[it->first] = n_par;

    fParLimits.push_back(it->second.GetThrowLimit());
    fParLimitsUp.push_back(it->second.GetThrowLimitUp());
    if (abs(it->second.GetCentral()) < 1.e-5) {
      parsHist.SetBinContent(n_par+1, val + 1.);
      if (fAddSystTerm) {
        //size_t cov_bin = fCovarianceBins[it->first];
        size_t new_cov_bin = fCovarianceBinsSimple[cov_bin_counter];
        //parsHist.SetBinError(n_par+1, sqrt((*fCovMatrixDisplay)[cov_bin][cov_bin]));
        parsHist.SetBinError(n_par+1, sqrt((*fCovMatrixDisplay)[new_cov_bin][new_cov_bin]));
        std::cout << "1 Set Error: " << it->first << " " <<
                     sqrt((*fCovMatrixDisplay)[new_cov_bin][new_cov_bin]) << std::endl;
      }
    }
    else {
      parsHist.SetBinContent(n_par+1,
          val/it->second.GetCentral());
      if (fAddSystTerm) {
        //size_t cov_bin = fCovarianceBins[it->first];
        size_t new_cov_bin = fCovarianceBinsSimple[cov_bin_counter];
        parsHist.SetBinError(n_par+1,
            //sqrt((*fCovMatrixDisplay)[cov_bin][cov_bin])/it->second.GetCentral());
            sqrt((*fCovMatrixDisplay)[new_cov_bin][new_cov_bin])/it->second.GetCentral());
        std::cout << "2 Set Error: " << it->first << " " <<
                     sqrt((*fCovMatrixDisplay)[new_cov_bin][new_cov_bin]) << " " << it->second.GetCentral() << " " <<
                     sqrt((*fCovMatrixDisplay)[new_cov_bin][new_cov_bin])/it->second.GetCentral() << std::endl;
      }
    }

    fPreFitParsNormal.SetBinContent(n_par+1, it->second.GetValue());
    fPreFitParsNormal.GetXaxis()->SetBinLabel(
        n_par+1, it->second.GetName().c_str());
    if (fAddSystTerm) {
      //size_t cov_bin = fCovarianceBins[it->first];
      size_t new_cov_bin = fCovarianceBinsSimple[cov_bin_counter];
      //std::cout << cov_bin << " " << new_cov_bin << std::endl;
      fPreFitParsNormal.SetBinError(
          n_par+1, sqrt((*fCovMatrixDisplay)[new_cov_bin][new_cov_bin]));
    }
    ++cov_bin_counter;
    ++n_par;
  }


  for (auto it = fG4RWParameters.begin(); it != fG4RWParameters.end(); ++it) {
    std::cout << "Adding g4rw parameter " << it->second.GetName().c_str() <<
                 " " << it->second.GetCentral() << " " <<
                 it->second.GetLowerLimit() << " " <<
                 it->second.GetUpperLimit() << std::endl;
    double val = it->second.GetCentral();
    if (fRandomStart) {
      it->second.SetValue(fRNG.Gaus(1., .1)*it->second.GetCentral());
      val = it->second.GetValue();
    }
    fMinimizer->SetVariable(n_par, it->second.GetName().c_str(),
                            val, 0.01);
    fMinimizerInitVals.push_back(val);
    if (fSetSystLimits && it->second.GetSetLimits()) {
      std::cout << "Setting limit " << it->second.GetLowerLimit() << " " <<
                   it->second.GetUpperLimit() << std::endl;
      fMinimizer->SetVariableLimits(n_par, it->second.GetLowerLimit(),
                                    it->second.GetUpperLimit());
    }

    fParLimits.push_back(it->second.GetThrowLimit());
    fParLimitsUp.push_back(it->second.GetThrowLimitUp());
    if (abs(it->second.GetCentral()) < 1.e-5) {
      parsHist.SetBinContent(n_par+1, val + 1.);
    }
    else {
      parsHist.SetBinContent(n_par+1,
          val/it->second.GetCentral());
    }

    fPreFitParsNormal.SetBinContent(n_par+1, it->second.GetValue());
    fPreFitParsNormal.GetXaxis()->SetBinLabel(
        n_par+1, it->second.GetName().c_str());
    if (fTuneG4RWPars) {
      fMinimizer->FixVariable(n_par);
      std::cout << "fixed " << n_par << " to " << val << std::endl;
    }
    ++n_par;
  }

  parsHist.SetMarkerColor(kBlue);
  parsHist.SetMarkerStyle(20);
  fOutputFile.cd();
  parsHist.Write();
  fPreFitParsNormal.Write();

}

void protoana::PDSPThinSliceFitter::Tune(std::string tune_file) {

  //TODO throw exception if not set fTuneG4RWPars or no G4RWParameters

  std::string line;
  std::ifstream myfile(tune_file);


  auto par = fG4RWParameters.begin();
  TVectorD tune_pars;
  size_t ipar = 0;
  if (myfile.is_open()) {
    while (std::getline(myfile, line)) {
      std::cout << line << '\n';
      par->second.SetValue(std::stof(line));
      par->second.SetCentral(std::stof(line));
      tune_pars.ResizeTo(ipar+1);
      tune_pars[ipar] = std::stof(line);
      ++par;
      ++ipar;
      std::cout << "set" << std::endl;
    }
    myfile.close();
  }
  else {
    std::exception e;
    throw e;
  }

  fOutputFile.cd();
  auto * dir = (TDirectory*)fOutputFile.mkdir("TunePars");
  dir->cd();
  tune_pars.Write("tune_pars");
  fOutputFile.cd();

}

void protoana::PDSPThinSliceFitter::InitializeMCSamples() {

  
  for (size_t i = 0; i < fSampleSets.size(); ++i) {
    fhicl::ParameterSet sample_set = fSampleSets[i];
    std::string sample_name = sample_set.get<std::string>("Name");
    int sample_ID = sample_set.get<int>("ID");

    if (fSamples.find(sample_ID) == fSamples.end()) {
      fSamples[sample_ID] = std::vector<std::vector<ThinSliceSample>>();
      fFakeSamples[sample_ID] = std::vector<std::vector<ThinSliceSample>>();

      if (fFitType == "ConstructCovariance") {
        fCovSamples[sample_ID] = std::vector<std::vector<ThinSliceSample>>();
      }
      fFluxesBySample[sample_ID] = std::vector<std::vector<double>>();
      fFakeFluxesBySample[sample_ID] = std::vector<std::vector<double>>();
    }
    fIsSignalSample[sample_ID] = sample_set.get<bool>("IsSignal");

    int flux_type = sample_set.get<int>("FluxType");
    std::vector<double> bins
        = sample_set.get<std::vector<double>>("SignalBins");
    fSignalBins[sample_ID] = bins;

    fSignalFixedBins[sample_ID]
        = sample_set.get<std::vector<size_t>>("FixedBins", {});

    bool single_signal_bin = sample_set.get<bool>("SingleSignalBin", false);
    auto signal_val = sample_set.get<double>("PrefitVal", 1.);
    std::cout << "GOT SIGNAL VAL: " << signal_val << std::endl;
    //

    for (size_t j = 1; j < fBeamEnergyBins.size(); ++j) {
      fSamples[sample_ID].push_back(std::vector<ThinSliceSample>());
      fFakeSamples[sample_ID].push_back(std::vector<ThinSliceSample>());
      if (fFitType == "ConstructCovariance") {
        fCovSamples[sample_ID].push_back(std::vector<ThinSliceSample>());
      }
      fFluxesBySample[sample_ID].push_back(std::vector<double>());
      fFakeFluxesBySample[sample_ID].push_back(std::vector<double>());
      if (sample_set.get<bool>("IsSignal")) {

        if (j == 1) {
          fSignalParameters[sample_ID] = std::vector<double>();
          fSignalParameterNames[sample_ID] = std::vector<std::string>();
        }

        if (!single_signal_bin) {
          ThinSliceSample underflow_sample(
              sample_name + "Underflow",
              flux_type, fSelectionSets,
              fIncidentRecoBins, fTrueIncidentBins, j);
          fSamples[sample_ID].back().push_back(underflow_sample);

          ThinSliceSample fake_underflow_sample(
              sample_name + "FakeUnderflow",
              flux_type, fSelectionSets,
              fIncidentRecoBins, fTrueIncidentBins, j);
          fFakeSamples[sample_ID].back().push_back(fake_underflow_sample);

          if (fFitType == "ConstructCovariance") {
            ThinSliceSample cov_underflow_sample(
                sample_name + "CovUnderflow",
                flux_type, fSelectionSets,
                fIncidentRecoBins, fTrueIncidentBins, j);
            fCovSamples[sample_ID].back().push_back(cov_underflow_sample);
          }

          if (j == 1 && fFitUnderOverflow) {
            fSignalParameters[sample_ID].push_back(signal_val);

            std::string par_name = "par_" + sample_name + "_underflow";
            fSignalParameterNames[sample_ID].push_back(par_name);
            ++fTotalSignalParameters;
          }

          fFluxesBySample[sample_ID].back().push_back(0.);
          fFakeFluxesBySample[sample_ID].back().push_back(0.);
        }
        else {
          ThinSliceSample single_sample(
              sample_name + "Single",
              flux_type, fSelectionSets,
              fIncidentRecoBins, fTrueIncidentBins, j);
          fSamples[sample_ID].back().push_back(single_sample);

          ThinSliceSample fake_single_sample(
              sample_name + "FakeSingle",
              flux_type, fSelectionSets,
              fIncidentRecoBins, fTrueIncidentBins, j);
          fFakeSamples[sample_ID].back().push_back(fake_single_sample);

          if (fFitType == "ConstructCovariance") {
            ThinSliceSample cov_single_sample(
                sample_name + "CovSingle",
                flux_type, fSelectionSets,
                fIncidentRecoBins, fTrueIncidentBins, j);
            fCovSamples[sample_ID].back().push_back(cov_single_sample);
          }

          if (j == 1) {
            fSignalParameters[sample_ID].push_back(signal_val);

            std::string par_name = "par_" + sample_name + "_single";
            fSignalParameterNames[sample_ID].push_back(par_name);
            ++fTotalSignalParameters;
          }

          fFluxesBySample[sample_ID].back().push_back(0.);
          fFakeFluxesBySample[sample_ID].back().push_back(0.);
        }

        for (size_t k = 1; k < bins.size(); ++k) {
          ThinSliceSample sample(sample_name, flux_type, fSelectionSets,
                                 fIncidentRecoBins, fTrueIncidentBins,
                                 j, true,
                                 {bins[k-1], bins[k]});

          std::string fake_name = sample_name + "Fake";
          ThinSliceSample fake_sample(fake_name, flux_type, fSelectionSets,
                                 fIncidentRecoBins, fTrueIncidentBins,
                                 j, true,
                                 {bins[k-1], bins[k]});
          fSamples[sample_ID].back().push_back(sample);
          fFakeSamples[sample_ID].back().push_back(fake_sample);

          if (fFitType == "ConstructCovariance") {
            std::string cov_name = sample_name + "Cov";
            ThinSliceSample cov_sample(cov_name, flux_type, fSelectionSets,
                                   fIncidentRecoBins, fTrueIncidentBins,
                                   j, true,
                                   {bins[k-1], bins[k]});
            fCovSamples[sample_ID].back().push_back(cov_sample);
          }

          if (j == 1) {
            fSignalParameters[sample_ID].push_back(signal_val);

            std::string par_name = "par_" + sample_name + "_" +
                                   PreciseToString(bins[k-1]) + "_" +
                                   PreciseToString(bins[k]);
            fSignalParameterNames[sample_ID].push_back(par_name);
            ++fTotalSignalParameters;
          }
          fFluxesBySample[sample_ID].back().push_back(0.);
          fFakeFluxesBySample[sample_ID].back().push_back(0.);
        }

        if (!single_signal_bin) {
        ThinSliceSample overflow_sample(sample_name + "Overflow",
                               flux_type, fSelectionSets,
                               fIncidentRecoBins, fTrueIncidentBins, j);
        fSamples[sample_ID].back().push_back(overflow_sample);

        ThinSliceSample fake_overflow_sample(sample_name + "FakeOverflow",
                               flux_type, fSelectionSets,
                               fIncidentRecoBins, fTrueIncidentBins, j);
        if (j == 1 && fFitUnderOverflow && !fTieUnderOver) {
          fSignalParameters[sample_ID].push_back(signal_val);

          std::string par_name = "par_" + sample_name + "_overflow";
          fSignalParameterNames[sample_ID].push_back(par_name);
          ++fTotalSignalParameters;
        }

        fFakeSamples[sample_ID].back().push_back(fake_overflow_sample);
        if (fFitType == "ConstructCovariance") {
          ThinSliceSample cov_overflow_sample(sample_name + "CovOverflow",
                                 flux_type, fSelectionSets,
                                 fIncidentRecoBins, fTrueIncidentBins, j);
          fCovSamples[sample_ID].back().push_back(cov_overflow_sample);
        }

        fFluxesBySample[sample_ID].back().push_back(0.);
        fFakeFluxesBySample[sample_ID].back().push_back(0.);
        }
      }
      else {
        ThinSliceSample sample(sample_name, flux_type, fSelectionSets,
                                       fIncidentRecoBins, fTrueIncidentBins, j);
        fSamples[sample_ID].back().push_back(sample);

        std::string fake_name = sample_name + "Fake";
        ThinSliceSample fake_sample(fake_name, flux_type, fSelectionSets,
                                       fIncidentRecoBins, fTrueIncidentBins, j);
        fFakeSamples[sample_ID].back().push_back(fake_sample);
        if (fFitType == "ConstructCovariance") {
          std::string cov_name = sample_name + "Cov";
          ThinSliceSample cov_sample(cov_name, flux_type, fSelectionSets,
                                         fIncidentRecoBins, fTrueIncidentBins, j);
          fCovSamples[sample_ID].back().push_back(cov_sample);
        }
        fFluxesBySample[sample_ID].back().push_back(0.);
        fFakeFluxesBySample[sample_ID].back().push_back(0.);
        if (fFluxParameters.find(flux_type) != fFluxParameters.end() &&
            j == 1) {
          fFluxParsToSamples[flux_type].push_back(sample_ID);
        }
      }
    }
  }

  MakeMinimizer();
}

void protoana::PDSPThinSliceFitter::NewBuildMC() {

  OpenFile(fMCFileName, fTreeName);

  std::cout << "Filling MC Events" << std::endl;

  int sample_ID, selection_ID, event, run, subrun;
  int cal_up_selection_ID, cal_down_selection_ID;
  int selection_ID_front_shift_up, selection_ID_front_shift_down,
      selection_ID_back_shift_up, selection_ID_back_shift_down;
  int true_beam_PDG;
  double true_beam_interactingEnergy, reco_beam_interactingEnergy, reco_beam_alt_len;
  double true_beam_endP, true_beam_mass, true_beam_endZ;
  double reco_beam_endZ, true_beam_startP, reco_beam_startY;
  double reco_beam_startX_SCE, reco_beam_startY_SCE, reco_beam_startZ_SCE;
  double reco_beam_endZ_SCE;
  double vertex_michel_score;
  int vertex_nhits;
  double beam_inst_P;
  double leading_p_costheta, leading_piplus_costheta, leading_pi0_costheta;
  double leading_p_momentum, leading_piplus_momentum, leading_pi0_momentum;
  std::vector<double> * reco_beam_incidentEnergies = 0x0,
                      * true_beam_incidentEnergies = 0x0,
                      * true_beam_traj_Z = 0x0,
                      * true_beam_traj_KE = 0x0,
                      * reco_daughter_track_thetas = 0x0,
                      * reco_daughter_track_scores = 0x0;
  fMCTree->SetBranchAddress("event", &event); //good
  fMCTree->SetBranchAddress("subrun", &subrun); //good
  fMCTree->SetBranchAddress("run", &run); //good

  fMCTree->SetBranchAddress("true_beam_PDG", &true_beam_PDG); //good

  // if (!fInclusive) { //good
    fMCTree->SetBranchAddress("new_interaction_topology", &sample_ID);

    // if (fNewFVSelection) {
      fMCTree->SetBranchAddress("selection_ID", &selection_ID);
    // }
    // else {
    //   fMCTree->SetBranchAddress("selection_ID", &selection_ID);
    // }
  // }
  // else {
  //   fMCTree->SetBranchAddress("inclusive_topology", &sample_ID);
  //   fMCTree->SetBranchAddress("selection_ID_inclusive", &selection_ID);
  // }
  fMCTree->SetBranchAddress("cal_up_selection_ID", &cal_up_selection_ID);
  fMCTree->SetBranchAddress("cal_down_selection_ID", &cal_down_selection_ID);

  fMCTree->SetBranchAddress("selection_ID_front_shift_up",
                        &selection_ID_front_shift_up);
  fMCTree->SetBranchAddress("selection_ID_front_shift_down",
                        &selection_ID_front_shift_down);

  fMCTree->SetBranchAddress("selection_ID_back_shift_up",
                        &selection_ID_back_shift_up);
  fMCTree->SetBranchAddress("selection_ID_back_shift_down",
                        &selection_ID_back_shift_down);

  fMCTree->SetBranchAddress("true_beam_interactingEnergy", //good
                        &true_beam_interactingEnergy);
  fMCTree->SetBranchAddress("true_beam_endP", &true_beam_endP); //good
  fMCTree->SetBranchAddress("true_beam_endZ", &true_beam_endZ); //good
  fMCTree->SetBranchAddress("true_beam_mass", &true_beam_mass); //good
  fMCTree->SetBranchAddress("reco_beam_interactingEnergy", //good
                        &reco_beam_interactingEnergy);
  fMCTree->SetBranchAddress("reco_beam_alt_len", &reco_beam_alt_len);
  fMCTree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ); //good
  fMCTree->SetBranchAddress("reco_beam_calo_startX", &reco_beam_startX_SCE);
  fMCTree->SetBranchAddress("reco_beam_calo_startY", &reco_beam_startY_SCE);
  fMCTree->SetBranchAddress("reco_beam_calo_startZ", &reco_beam_startZ_SCE);
  fMCTree->SetBranchAddress("reco_beam_calo_endZ", &reco_beam_endZ_SCE);
  fMCTree->SetBranchAddress("reco_beam_vertex_michel_score", &vertex_michel_score);
  fMCTree->SetBranchAddress("reco_beam_vertex_nHits", &vertex_nhits);
  fMCTree->SetBranchAddress("reco_beam_startY", &reco_beam_startY); //good
  fMCTree->SetBranchAddress("reco_beam_incidentEnergies", //good
                        &reco_beam_incidentEnergies);
  fMCTree->SetBranchAddress("true_beam_incidentEnergies", //good
                        &true_beam_incidentEnergies);

  fMCTree->SetBranchAddress("true_beam_startP", &true_beam_startP); //good
  fMCTree->SetBranchAddress("true_beam_traj_Z", &true_beam_traj_Z); //good
  fMCTree->SetBranchAddress("true_beam_traj_KE", &true_beam_traj_KE); //good

  fMCTree->SetBranchAddress("beam_inst_P", &beam_inst_P); //good
  fMCTree->SetBranchAddress("reco_daughter_allTrack_Theta", &reco_daughter_track_thetas); //good
  fMCTree->SetBranchAddress("reco_daughter_PFP_trackScore_collection", //good
                            &reco_daughter_track_scores);
  std::vector<int> * true_beam_daughter_PDG = 0x0;
  fMCTree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);//good

  int true_beam_ID, reco_beam_true_byHits_ID;
  fMCTree->SetBranchAddress("true_beam_ID", &true_beam_ID); //good
  fMCTree->SetBranchAddress("reco_beam_true_byHits_ID", &reco_beam_true_byHits_ID); //good

  std::vector<std::vector<double>> * g4rw_full_grid_proton_coeffs = 0x0,
                                  * g4rw_full_grid_piplus_coeffs = 0x0,
                                  * g4rw_full_fine_piplus_coeffs = 0x0,
                                  * g4rw_full_grid_abscex_coeffs = 0x0,
                                  * g4rw_primary_grid_abscex_coeffs = 0x0,
                                  * g4rw_downstream_grid_piplus_coeffs = 0x0;
  fMCTree->SetBranchAddress("g4rw_full_grid_proton_coeffs", //good
                        &g4rw_full_grid_proton_coeffs);
  fMCTree->SetBranchAddress("g4rw_full_grid_piplus_coeffs",
                        &g4rw_full_grid_piplus_coeffs);
  fMCTree->SetBranchAddress("g4rw_full_fine_piplus_coeffs",
                        &g4rw_full_fine_piplus_coeffs);
  fMCTree->SetBranchAddress("g4rw_full_grid_abscex_coeffs",
                        &g4rw_full_grid_abscex_coeffs);
  fMCTree->SetBranchAddress("g4rw_primary_grid_abscex_coeffs",
                        &g4rw_primary_grid_abscex_coeffs);
  fMCTree->SetBranchAddress("g4rw_downstream_grid_piplus_coeffs", //good
                        &g4rw_downstream_grid_piplus_coeffs);

  std::vector<std::vector<double>> * daughter_dQdXs = 0x0,
                                  * daughter_resRanges = 0x0,
                                  * daughter_EFields = 0x0;
  fMCTree->SetBranchAddress(
      "reco_daughter_allTrack_dQdX_SCE", &daughter_dQdXs);//good
  fMCTree->SetBranchAddress(
      "reco_daughter_allTrack_resRange_SCE", &daughter_resRanges);//good
  fMCTree->SetBranchAddress(
      "reco_daughter_allTrack_EField_SCE", &daughter_EFields);//good
  bool has_pi0_shower;
  fMCTree->SetBranchAddress("has_shower_dist_energy", &has_pi0_shower);//good 
  std::vector<double> * true_beam_daughter_startP = 0x0;
  fMCTree->SetBranchAddress("true_beam_daughter_startP",  &true_beam_daughter_startP);//good
  fMCTree->SetBranchAddress("leading_p_costheta", &leading_p_costheta);//good
  fMCTree->SetBranchAddress("leading_piplus_costheta", &leading_piplus_costheta);//good
  fMCTree->SetBranchAddress("leading_pi0_costheta", &leading_pi0_costheta);//good
  fMCTree->SetBranchAddress("leading_p_momentum", &leading_p_momentum);//good
  fMCTree->SetBranchAddress("leading_piplus_momentum", &leading_piplus_momentum);//good
  fMCTree->SetBranchAddress("leading_pi0_momentum", &leading_pi0_momentum);//good
  bool is_beam_scraper;
  fMCTree->SetBranchAddress("true_beam_is_scraper", &is_beam_scraper);//good

  std::vector<double> * reco_daughter_chi2 = 0x0,
                      * reco_daughter_truncated_dEdX = 0x0;
  std::vector<int> * reco_daughter_chi2_nhits = 0x0;
  fMCTree->SetBranchAddress("reco_daughter_allTrack_Chi2_proton",
                        &reco_daughter_chi2);
  fMCTree->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof",
                        &reco_daughter_chi2_nhits);
  fMCTree->SetBranchAddress("reco_daughter_allTrack_truncLibo_dEdX_pos",
                        &reco_daughter_truncated_dEdX);
  int true_n_neutrons, true_n_protons, true_n_piplus, true_n_piminus,
      true_n_pi0;
  fMCTree->SetBranchAddress("true_daughter_nNeutron", &true_n_neutrons);
  fMCTree->SetBranchAddress("true_daughter_nProton", &true_n_protons);
  fMCTree->SetBranchAddress("true_daughter_nPiPlus", &true_n_piplus);
  fMCTree->SetBranchAddress("true_daughter_nPiMinus", &true_n_piminus);
  fMCTree->SetBranchAddress("true_daughter_nPi0", &true_n_pi0);

  std::vector<double> * reco_daughter_shower_energy = 0x0;
  fMCTree->SetBranchAddress("reco_daughter_allShower_energy",
                        &reco_daughter_shower_energy);

  int nentries = (fMaxEntries < 0 ? fMCTree->GetEntries() : fMaxEntries);
  if(fMaxEntries > fMCTree->GetEntries()) {
    nentries = fMCTree->GetEntries();
    std::string message = "Requested more entries than in MC fMCTree. Setting nentries to fMCTree's numer of entries.";
  }

  fSplitVal = nentries;

  int events_end = nentries;
  int fake_start = 0;

  std::vector<int> event_list, fake_event_list;

  if (fSplitMC) {
    fSplitVal = nentries/2;
    events_end = nentries/2;
    fake_start = events_end;
    std::cout << "Note: Splitting MC in half. " <<
                fSplitVal << "/" << nentries <<std::endl;
    if (fShuffle) {
      for (int i = 0; i < nentries; ++i) {
        event_list.push_back(i);
      }
      for (int i = 0; i < events_end; ++i) {
        int r = fRNG.Integer(event_list.size());
        fake_event_list.push_back(event_list[r]);
        event_list.erase(event_list.begin()+r);
      }

      std::sort(event_list.begin(), event_list.end());
      std::sort(fake_event_list.begin(), fake_event_list.end());

      std::cout << "Shuffled " << event_list.size() << " " <<
                  fake_event_list.size() << std::endl;
      for (int & i : event_list) {
        if (std::find(fake_event_list.begin(), fake_event_list.end(), i) !=
            fake_event_list.end()) {
          std::cout << "Incorrectly shuffled" << std::endl;
        }
      }
    }
    else {
      for (int i = 0; i < events_end; ++i) {
        event_list.push_back(i);
      }
      for (int i = fake_start; i < nentries; ++i) {
        fake_event_list.push_back(i);
      }
      std::cout << "Split " << event_list.size() << " " <<
                  fake_event_list.size() << std::endl;
    }
  }
  else {
    int nfake_entries = (fMaxDataEntries < 0 ? fMCTree->GetEntries() : fMaxDataEntries);
    if (fMaxDataEntries > fMCTree->GetEntries()) {
        std::string message = "Requested more fake_entries than in MC fMCTree";
        throw std::runtime_error(message);
    }
    for (int i = 0; i < events_end; ++i) {
      event_list.push_back(i);
    }
    for (int i = 0; i < nfake_entries; ++i) {
      fake_event_list.push_back(i);
    }
  }

  for (int & i : event_list) {
    if (!(i % 20000)) std::cout << i << "/" << fSplitVal << std::endl;
    fMCTree->GetEntry(i);
  
    //ThinSliceEvent thin_slice_event(event, subrun, run);
    fEvents.push_back(ThinSliceEvent(event, subrun, run));
    fEvents.back().SetMCStatVarWeight(fRNG.Poisson(1));
    fEvents.back().SetSampleID(sample_ID);
    fEvents.back().SetSelectionID(selection_ID);
    fEvents.back().SetCalUpSelectionID(cal_up_selection_ID);
    fEvents.back().SetCalDownSelectionID(cal_down_selection_ID);
  
    fEvents.back().SetFrontUpSelectionID(selection_ID_front_shift_up);
    fEvents.back().SetFrontDownSelectionID(selection_ID_front_shift_down);
    fEvents.back().SetBackUpSelectionID(selection_ID_back_shift_up);
    fEvents.back().SetBackDownSelectionID(selection_ID_back_shift_down);
  
    fEvents.back().SetTrueInteractingEnergy(true_beam_interactingEnergy);
    // if (fDoEnergyByLen) {
    //   fEvents.back().SetRecoInteractingEnergy((*reco_beam_incidentEnergies)[0] - 2.1*reco_beam_alt_len);
    // }
    // else {
      fEvents.back().SetRecoInteractingEnergy(reco_beam_interactingEnergy);
    // }
    fEvents.back().SetTrueEndP(true_beam_endP);
    fEvents.back().SetTrueEndZ(true_beam_endZ);
    fEvents.back().SetTrueStartP(true_beam_startP);
    fEvents.back().SetTrueMass(true_beam_mass);
    fEvents.back().SetRecoEndZ(reco_beam_endZ);
    fEvents.back().SetRecoStartY(reco_beam_startY);
    fEvents.back().SetRecoStartX_SCE(reco_beam_startX_SCE);
    fEvents.back().SetRecoStartY_SCE(reco_beam_startY_SCE);
    fEvents.back().SetRecoStartZ_SCE(reco_beam_startZ_SCE);
    fEvents.back().SetRecoEndZ_SCE(reco_beam_endZ_SCE);
  
    fEvents.back().SetVertexMichelScore(vertex_michel_score);
    fEvents.back().SetVertexNHits(vertex_nhits);
  
    fEvents.back().SetTrueID(true_beam_ID);
    fEvents.back().SetRecoToTrueID(reco_beam_true_byHits_ID);
  
    fEvents.back().SetRecoIncidentEnergies(*reco_beam_incidentEnergies);
    fEvents.back().SetTrueInitEnergy(
      (true_beam_incidentEnergies->size() > 0 ?
       true_beam_incidentEnergies->at(0) : -999.)); // TODO -- CHANGE THIS WITH NEW VAL 
    fEvents.back().SetTrueTrajZ(*true_beam_traj_Z);
    fEvents.back().SetTrueTrajKE(*true_beam_traj_KE);
    fEvents.back().SetTrueNNeutrons(true_n_neutrons);
    fEvents.back().SetTrueNProtons(true_n_protons);
    fEvents.back().SetTrueNPiPlus(true_n_piplus);
    fEvents.back().SetTrueNPiMinus(true_n_piminus);
    fEvents.back().SetTrueNPi0(true_n_pi0);
  
    fEvents.back().SetBeamInstP(beam_inst_P);
    fEvents.back().SetPDG(true_beam_PDG);
    fEvents.back().SetRecoDaughterTrackScores(*reco_daughter_track_scores);
    fEvents.back().SetHasPi0Shower(has_pi0_shower);
    fEvents.back().SetLeadingPCostheta(leading_p_costheta);
    fEvents.back().SetLeadingPiPlusCostheta(leading_piplus_costheta);
    fEvents.back().SetLeadingPi0Costheta(leading_pi0_costheta);
    fEvents.back().SetLeadingPMomentum(leading_p_momentum);
    fEvents.back().SetLeadingPiPlusMomentum(leading_piplus_momentum);
    fEvents.back().SetLeadingPi0Momentum(leading_pi0_momentum);
    fEvents.back().SetIsBeamScraper(is_beam_scraper);
  
    fEvents.back().SetRecoShowerEnergy(*reco_daughter_shower_energy);
  
    for (size_t j = 0; j < daughter_dQdXs->size(); ++j) {
      fEvents.back().AddRecoDaughterEField((*daughter_EFields)[j]);
      auto nhits = (*reco_daughter_chi2_nhits)[j];
      fEvents.back().AddOneChi2PerHit(
          (nhits > 0) ?
          (*reco_daughter_chi2)[j]/nhits :
          -999.
      );
      fEvents.back().AddOneTrunc_dEdX(
        (*reco_daughter_truncated_dEdX)[j]
      );
    }
  
    for (size_t j = 0; j < g4rw_downstream_grid_piplus_coeffs->size(); ++j) {
      std::string name_downstream = "g4rw_downstream_grid_piplus_coeffs_" +
                                   std::to_string(j);
      fEvents.back().MakeG4RWCoeff(name_downstream,
                                  (*g4rw_downstream_grid_piplus_coeffs)[j]);
  
      std::string name_full = "g4rw_full_grid_piplus_coeffs_" + std::to_string(j);
      fEvents.back().MakeG4RWCoeff(name_full, (*g4rw_full_grid_piplus_coeffs)[j]);
    }
  
    for (size_t j = 0; j < g4rw_full_fine_piplus_coeffs->size(); ++j) {
      std::string fine_full = "g4rw_full_fine_piplus_coeffs_" + std::to_string(j);
      fEvents.back().MakeG4RWCoeff(fine_full, (*g4rw_full_fine_piplus_coeffs)[j]);
    }
  
    for (size_t j = 0; j < g4rw_full_grid_abscex_coeffs->size(); ++j) {
      std::string abscex = "g4rw_full_grid_abscex_coeffs_" + std::to_string(j);
      fEvents.back().MakeG4RWCoeff(abscex, (*g4rw_full_grid_abscex_coeffs)[j]);
    }
  
  
    fEvents.back().MakeG4RWCoeff("g4rw_full_grid_proton_coeffs",
                                (*g4rw_full_grid_proton_coeffs)[0]);
  
    bool found_start = false;
    for (size_t j = 1; j < true_beam_traj_Z->size(); ++j) {
       
      if ((*true_beam_traj_Z)[j] < fTrajZStart) continue;
  
      double delta = (*true_beam_traj_KE)[0] - (*true_beam_traj_KE)[j];
      fEvents.back().SetDeltaEToTPC(delta);
      found_start = true;
      break;
    }
    if (!found_start) {
      fEvents.back().SetDeltaEToTPC(-999.);
    }
  
  }
  
  if (fDoFakeData) {
    std::cout << "Filling fake data " << fake_event_list.size() << std::endl;
    //for (int i = fake_start; i < fake_end; ++i) {
    auto start_time = std::chrono::high_resolution_clock::now();
    for (int & i : fake_event_list) {
      //std::cout << "Fake event " << i << std::endl;
      fMCTree->GetEntry(i);
  
      fFakeDataEvents.push_back(ThinSliceEvent(event, subrun, run));
      if (!(fFakeDataEvents.size() % 20000))
        std::cout << fFakeDataEvents.size() << std::endl;
      fFakeDataEvents.back().SetMCStatVarWeight(fRNG.Poisson(1));
      fFakeDataEvents.back().SetSampleID(sample_ID);
      fFakeDataEvents.back().SetSelectionID(selection_ID);
      fFakeDataEvents.back().SetCalUpSelectionID(cal_up_selection_ID);
      fFakeDataEvents.back().SetCalDownSelectionID(cal_down_selection_ID);
  
      fFakeDataEvents.back().SetFrontUpSelectionID(selection_ID_front_shift_up);
      fFakeDataEvents.back().SetFrontDownSelectionID(selection_ID_front_shift_down);
      fFakeDataEvents.back().SetBackUpSelectionID(selection_ID_back_shift_up);
      fFakeDataEvents.back().SetBackDownSelectionID(selection_ID_back_shift_down);
  
      fFakeDataEvents.back().SetTrueInteractingEnergy(true_beam_interactingEnergy);
      // if (fDoEnergyByLen) {
      //   fFakeDataEvents.back().SetRecoInteractingEnergy((*reco_beam_incidentEnergies)[0] - 2.1*reco_beam_alt_len);
      // }
      // else {
        fFakeDataEvents.back().SetRecoInteractingEnergy(reco_beam_interactingEnergy);
      // }
      fFakeDataEvents.back().SetTrueEndP(true_beam_endP);
      fFakeDataEvents.back().SetTrueEndZ(true_beam_endZ);
      fFakeDataEvents.back().SetTrueStartP(true_beam_startP);
      fFakeDataEvents.back().SetTrueMass(true_beam_mass);
      fFakeDataEvents.back().SetRecoEndZ(reco_beam_endZ);
      fFakeDataEvents.back().SetRecoStartX_SCE(reco_beam_startX_SCE);
      fFakeDataEvents.back().SetRecoStartY_SCE(reco_beam_startY_SCE);
      fFakeDataEvents.back().SetRecoStartZ_SCE(reco_beam_startZ_SCE);
      fFakeDataEvents.back().SetRecoEndZ_SCE(reco_beam_endZ_SCE);
      fFakeDataEvents.back().SetVertexMichelScore(vertex_michel_score);
      fFakeDataEvents.back().SetVertexNHits(vertex_nhits);
  
      fFakeDataEvents.back().SetRecoIncidentEnergies(*reco_beam_incidentEnergies);
      fFakeDataEvents.back().SetTrueInitEnergy(
          (true_beam_incidentEnergies->size() > 0 ?
           true_beam_incidentEnergies->at(0) : -999.));
      fFakeDataEvents.back().SetTrueTrajZ(*true_beam_traj_Z);
      fFakeDataEvents.back().SetTrueTrajKE(*true_beam_traj_KE);
  
      fFakeDataEvents.back().SetBeamInstP(beam_inst_P);
      fFakeDataEvents.back().SetPDG(true_beam_PDG);
      if (fFakeDataRoutine == "EffVar") {
        fFakeDataEvents.back().SetRecoDaughterTrackThetas(*reco_daughter_track_thetas); }
      fFakeDataEvents.back().SetTrueDaughterPDGs(*true_beam_daughter_PDG);
      fFakeDataEvents.back().SetTrueDaughterStartPs(*true_beam_daughter_startP);
      fFakeDataEvents.back().SetRecoDaughterTrackScores(*reco_daughter_track_scores);
      fFakeDataEvents.back().SetHasPi0Shower(has_pi0_shower);
      fFakeDataEvents.back().SetLeadingPCostheta(leading_p_costheta);
      fFakeDataEvents.back().SetLeadingPiPlusCostheta(leading_piplus_costheta);
      fFakeDataEvents.back().SetLeadingPi0Costheta(leading_pi0_costheta);
      fFakeDataEvents.back().SetLeadingPMomentum(leading_p_momentum);
      fFakeDataEvents.back().SetLeadingPiPlusMomentum(leading_piplus_momentum);
      fFakeDataEvents.back().SetLeadingPi0Momentum(leading_pi0_momentum);
      fFakeDataEvents.back().SetIsBeamScraper(is_beam_scraper);
      fFakeDataEvents.back().SetRecoShowerEnergy(*reco_daughter_shower_energy);
  
      fFakeDataEvents.back().SetTrueID(true_beam_ID);
      fFakeDataEvents.back().SetTrueNNeutrons(true_n_neutrons);
      fFakeDataEvents.back().SetTrueNProtons(true_n_protons);
      fFakeDataEvents.back().SetTrueNPiPlus(true_n_piplus);
      fFakeDataEvents.back().SetTrueNPiMinus(true_n_piminus);
      fFakeDataEvents.back().SetTrueNPi0(true_n_pi0);
      fFakeDataEvents.back().SetRecoToTrueID(reco_beam_true_byHits_ID);
  
      for (size_t j = 0; j < daughter_dQdXs->size(); ++j) {
  
        fFakeDataEvents.back().AddRecoDaughterEField((*daughter_EFields)[j]);
        auto nhits = (*reco_daughter_chi2_nhits)[j];
        fFakeDataEvents.back().AddOneChi2PerHit(
            (nhits > 0) ?
            (*reco_daughter_chi2)[j]/nhits :
            -999.
        );
        fFakeDataEvents.back().AddOneTrunc_dEdX(
          (*reco_daughter_truncated_dEdX)[j]
        );
      }
  
      for (size_t j = 0; j < g4rw_full_fine_piplus_coeffs->size(); ++j) {
        std::string fine_full = "g4rw_full_fine_piplus_coeffs_" + std::to_string(j);
        fFakeDataEvents.back().MakeG4RWCoeff(fine_full, (*g4rw_full_fine_piplus_coeffs)[j]);
      }
      for (size_t j = 0; j < g4rw_full_grid_abscex_coeffs->size(); ++j) {
        std::string abscex = "g4rw_full_grid_abscex_coeffs_" + std::to_string(j);
        fFakeDataEvents.back().MakeG4RWCoeff(abscex, (*g4rw_full_grid_abscex_coeffs)[j]);
  
      }
  
      for (size_t j = 0; j < g4rw_downstream_grid_piplus_coeffs->size(); ++j) {
        std::string name_full = "g4rw_full_grid_piplus_coeffs_" + std::to_string(j);
        fFakeDataEvents.back().MakeG4RWCoeff(name_full, (*g4rw_full_grid_piplus_coeffs)[j]);
  
  
        std::string name_downstream = "g4rw_downstream_grid_piplus_coeffs_" +
                                     std::to_string(j);
        fFakeDataEvents.back().MakeG4RWCoeff(name_downstream,
                                    (*g4rw_downstream_grid_piplus_coeffs)[j]);
      }
  
      fFakeDataEvents.back().MakeG4RWCoeff("g4rw_full_grid_proton_coeffs",
                                            (*g4rw_full_grid_proton_coeffs)[0]);
      bool found_start = false;
      for (size_t j = 1; j < true_beam_traj_Z->size(); ++j) {
        if ((*true_beam_traj_Z)[j-1] < fTrajZStart &&
            fTrajZStart < (*true_beam_traj_Z)[j+1]) {
          double delta = (*true_beam_traj_KE)[0] - (*true_beam_traj_KE)[j];
          fFakeDataEvents.back().SetDeltaEToTPC(delta);
          found_start = true;
          break;
        }
      }
      if (!found_start) fFakeDataEvents.back().SetDeltaEToTPC(-999.);
    }
    auto new_time = std::chrono::high_resolution_clock::now();
    auto delta =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            new_time - start_time).count();
    std::cout << "Filling " << fFakeDataEvents.size() <<
                 " fake data took " << delta << std::endl;
  }
  
  std::cout << "Filled MC Events" << std::endl;
 
  
  CloseFile();
}

//Good
void protoana::PDSPThinSliceFitter::BuildMCSamples() {
  NewBuildMC();
  //Open the MC file and set branches
  // fThinSliceDriver->OpenFile(fMCFileName, fTreeName);

  // std::cout << "FakeData: " << fDoFakeData << " Routine: " <<
  //              fFakeDataRoutine << std::endl;
  // fThinSliceDriver->FillMCEvents(fEvents, fFakeDataEvents, fSplitVal,
  //                                fSplitMC, fShuffle, fMaxEntries,
  //                                fMaxDataEntries, fDoFakeData);
  // fThinSliceDriver->BuildMCSamples(fEvents, fSamples, fIsSignalSample,
  //                                  fNominalFluxes, fFluxesBySample,
  //                                  fBeamEnergyBins, fScaleToDataBeamProfile);
  // fThinSliceDriver->BuildMCSamples(fFakeDataEvents, fFakeSamples, fIsSignalSample,
  //                                  fFakeFluxes, fFakeFluxesBySample,
  //                                  fBeamEnergyBins, fScaleToDataBeamProfile);
  // fThinSliceDriver->CloseFile();

  // if (fDoSysts) {
  //   std::cout << "Setting up systs" << std::endl;
  //   fThinSliceDriver->SetupSysts(fEvents, fSamples, fIsSignalSample,
  //                                fBeamEnergyBins, fSystParameters,
  //                                fG4RWParameters,
  //                                fOutputFile);
  // }

}

void protoana::PDSPThinSliceFitter::BuildDataHists() {
  TFile * inputFile = TFile::Open((!fDoFakeData ? fDataFileName.c_str() : fMCFileName.c_str()),
                                  "OPEN");
  TTree * tree = (TTree*)inputFile->Get(fTreeName.c_str());
  if (!fDoFakeData && fRefitFile == "") {
    fThinSliceDriver->BuildDataHists(
        tree, fDataSet, fDataFlux,
        fBeamEnergyBins,
        fDataBeamFluxes,
        fMaxDataEntries);
    std::cout << "Data Beam fluxes: ";
    for (const auto & f : fDataBeamFluxes) {
      std::cout << f << " ";
    }
    std::cout << std::endl;
  }
  else if (fRefitFile != "") {

    std::string line;
    std::ifstream myfile(fRefitFile);

    if (myfile.is_open()) {
      while (std::getline(myfile, line)) {
        std::cout << line << '\n';
        break;
      }
      myfile.close();
    }
    else {
      std::exception e;
      std::cerr << "Could not get line from " << fRefitFile << std::endl;
      throw e;
    }

    if (line.find("/pnfs") == 0) {
      std::cout << "Replacing" << std::endl;
      line.replace(0, 5, "root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr");
    }
    TFile * refit_file = TFile::Open(line.c_str());    
    auto & sel_hists = fDataSet.GetSelectionHists();
    for (const auto & sel_set : fSelectionSets) {
      std::string name = sel_set.get<std::string>("Name");
      TString hname = TString::Format("Data/Data_selected_%s_hist", name.c_str());
      int id = sel_set.get<int>("ID");
      std::cout << id << " " << name << " " <<
                   refit_file->Get(hname) << std::endl;

      sel_hists.at(id) = (TH1*)refit_file->Get(hname);
      sel_hists.at(id)->SetDirectory(0);

    }
    TVectorD * fluxes = (TVectorD*)refit_file->Get("Data/data_fluxes");
    for (int i = 0; i < fluxes->GetNrows(); ++i) {
      fDataBeamFluxes.push_back((*fluxes)[i]);
    }

    std::cout << "Data Beam fluxes: ";
    for (const auto & f : fDataBeamFluxes) {
      std::cout << f << " ";
    }
    std::cout << std::endl;

    if (fDoFakeData) {
      auto * dir = (TDirectory*)refit_file->Get("FakeDataXSecs");
      if (dir) {
        //std::vector<TH1*> fake_hists;
        for (auto s : fMeasurementSamples) {
          auto & samples_vec_2D = fSamples[s];
          std::string xsec_name = "FakeData" +
                                   samples_vec_2D[0][1].GetName() + "FakeXSec";
          std::cout << "Getting " << xsec_name << std::endl;
          fFakeDataXSecs[s] = (TH1*)dir->Get(xsec_name.c_str());
          fFakeDataXSecs[s]->SetDirectory(0);
        }
      }
    }

    std::vector<double> vals;
    for (auto it = fSignalParameters.begin();
         it != fSignalParameters.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        vals.push_back(it->second[i]);
      }
    }
    for (auto it = fFluxParameters.begin();
         it != fFluxParameters.end(); ++it) {
      vals.push_back(it->second);
    }

    for (auto it = fSystParameters.begin();
         it != fSystParameters.end(); ++it) {
      vals.push_back(it->second.GetValue());
    }

    for (auto it = fG4RWParameters.begin();
         it != fG4RWParameters.end(); ++it) {
      vals.push_back(it->second.GetValue());
    }
    fThinSliceDriver->TurnOffFakeData();
    fFillIncidentInFunction = true;
    fUseFakeSamples = false;
    fThinSliceDriver->SetStatVar(false);
    fThinSliceDriver->SetFillFakeInMain(false);
    //Refill the hists for comparisons
    fFitFunction(&vals[0]); 
    fFillIncidentInFunction = false;

    fOutputFile.cd();
    if (fDoFakeData) {
      TDirectory * out = (TDirectory *)fOutputFile.mkdir("FakeDataXSecs");
      out->cd();
      for (auto it = fFakeDataXSecs.begin(); it != fFakeDataXSecs.end(); ++it) {
        it->second->Write();
      }
    }

    //Also add previous best fit xsecs
    std::cout << "Saving past xsecs " << std::endl;
    fOutputFile.cd();
    auto * out = (TDirectory *)fOutputFile.mkdir("PastFitXSecs");
    std::cout << out << std::endl;
    auto * past_xsecs = (TDirectory *)refit_file->Get("PostFitXSec");
    auto * keys = past_xsecs->GetListOfKeys();
    
    for (int i = 0; i < keys->GetSize(); ++i) {
      auto * xsec = (TH1D*)past_xsecs->Get(keys->At(i)->GetName());
      std::cout << keys->At(i)->GetName() << " " << xsec << std::endl;
      xsec->SetDirectory(0);
      out->cd();
      xsec->Write();
      std::cout << "Wrote" << std::endl;
    }

    std::cout << "Done" << std::endl;
    refit_file->Close();
    fOutputFile.cd();

  }
  else if (fFakeDataRoutine == "Toy") {
    BuildDataFromToy();
  }
  else if (fFakeDataRoutine == "Asimov") {
    std::cout << "Doing Asimov" << std::endl;
    for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
      std::vector<ThinSliceSample> & samples_vec = it->second[0];
      fFakeDataScales[it->first] = std::vector<double>(samples_vec.size(), 1.);
    }

    std::vector<double> vals;
    for (auto it = fSignalParameters.begin();
         it != fSignalParameters.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        vals.push_back(it->second[i]);
      }
    }
    for (auto it = fFluxParameters.begin();
         it != fFluxParameters.end(); ++it) {
      vals.push_back(it->second);
    }

    for (auto it = fSystParameters.begin();
         it != fSystParameters.end(); ++it) {
      vals.push_back(it->second.GetValue());
    }

    std::vector<double> g4rw_pars;
    for (auto it = fG4RWParameters.begin();
         it != fG4RWParameters.end(); ++it) {
      std::cout << "G4RW Fake Data: " << (fNominalG4RWFakeData ? 1. : it->second.GetValue()) << std::endl;
      vals.push_back((fNominalG4RWFakeData ? 1. : it->second.GetValue()));
      g4rw_pars.push_back(it->second.GetValue());
    }

    fFillIncidentInFunction = true;
    fUseFakeSamples = true;
    fThinSliceDriver->TurnOnFakeData();
    std::cout << "Vary stats: " << fVaryMCStatsForFakeData << std::endl;
    std::cout << fDataBeamFluxes.size() << std::endl;
    fThinSliceDriver->SetStatVar(fVaryMCStatsForFakeData); //Force on stat vars in driver
    fThinSliceDriver->SetFillFakeInMain(true);
    fFitFunction(&vals[0]);
    fDataSet.FillHistsFromSamples(fFakeSamples, fDataFlux, fDataBeamFluxes,
                                  (fDoFluctuateStats && fFluctuateInSamples));
    std::cout << "Asimov fluxes: " << std::endl;
    for (const auto & f: fDataBeamFluxes) std::cout << f << " ";
    std::cout << std::endl;
    BuildFakeDataXSecs(false);

    //Reset the parameters
    vals.clear();
    for (auto it = fSignalParameters.begin();
         it != fSignalParameters.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        vals.push_back(it->second[i]);
      }
    }
    for (auto it = fFluxParameters.begin();
         it != fFluxParameters.end(); ++it) {
      vals.push_back(it->second);
    }

    for (auto it = fSystParameters.begin();
         it != fSystParameters.end(); ++it) {
      vals.push_back(it->second.GetValue());
    }

    size_t ipar = 0;
    for (auto it = fG4RWParameters.begin();
         it != fG4RWParameters.end(); ++it) {
      it->second.SetValue(g4rw_pars[ipar]);
      vals.push_back(it->second.GetValue());
      std::cout << "G4RW Fake Data: " << it->second.GetValue() << std::endl;
      ++ipar;
    }
    fThinSliceDriver->TurnOffFakeData();
    fUseFakeSamples = false;
    fThinSliceDriver->SetStatVar(false);
    fThinSliceDriver->SetFillFakeInMain(false);
    //Refill the hists for comparisons
    fFitFunction(&vals[0]); 
    fFillIncidentInFunction = false;
  }
  else if (fFakeDataRoutine == "G4RWGrid") { // For now just do this one. Replace with general
    std::cout << "Building G4RW Grid" << std::endl;
    for (auto it = fFakeSamples.begin(); it != fFakeSamples.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        for (size_t j = 0; j < it->second[i].size(); ++j) {
          it->second[i][j].Reset();
        }
      }
    }

    std::vector<double> nominal_vals; //Will eventually Set if needed
    for (auto it = fSignalParameters.begin();
         it != fSignalParameters.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        nominal_vals.push_back(it->second[i]);
      }
    }
    for (auto it = fFluxParameters.begin();
         it != fFluxParameters.end(); ++it) {
      nominal_vals.push_back(it->second);
    }

    for (auto it = fSystParameters.begin();
         it != fSystParameters.end(); ++it) {
      nominal_vals.push_back(it->second.GetValue());
    }

    std::pair<size_t, size_t> g4rw_range
        = {nominal_vals.size(), nominal_vals.size()};
    for (auto it = fG4RWParameters.begin();
         it != fG4RWParameters.end(); ++it) {
      nominal_vals.push_back(it->second.GetValue());
      g4rw_range.second++;
    }

    auto vals = nominal_vals;
    auto set_vals = fAnalysisOptions.get<std::vector<std::pair<int, double>>>(
      "SetValsFakeData", {}
    );
    //ADD CHECK HERE FOR SIZE
    for (const auto & par : set_vals) {
      std::cout << "Setting: " << par.first << " " << par.second << std::endl;
      vals[par.first] = par.second;
    }
    for (size_t i = 0; i < vals.size(); ++i) {
      std::cout << i << " " << vals[i] << std::endl;
    }

    if (fTuneG4RWPars) {//Reset these now
      std::cout << "Resetting" << std::endl;
      for (size_t i = g4rw_range.first; i < g4rw_range.second; ++i) {
        vals[i] = 1.;
        std::cout << i << " " << vals[i] << std::endl;
      }
    }

    //can replace w/ below
    auto throw_syst = fAnalysisOptions.get<int>("ThrowSystFakeData", -1);
    if (throw_syst >= 0) {

      size_t nonsyst_pars = fTotalSignalParameters + fTotalFluxParameters;
      std::cout << "nonsyst: " << nonsyst_pars << std::endl;
      size_t total_pars = nonsyst_pars + fTotalSystParameters;
      std::vector<double> all_systs(nominal_vals.begin()+nonsyst_pars,
                                    nominal_vals.begin()+total_pars);

      TMatrixD * cov_lower = (TMatrixD*)fInputChol->GetU().Clone();
      cov_lower->Transpose(fInputChol->GetU());

      size_t n_cov_rows = cov_lower->GetNrows();
      if (n_cov_rows != fTotalSystParameters) {
        std::cout << "Error: " << n_cov_rows << " " << fTotalSystParameters <<
                     std::endl;
      }

      std::vector<double> throw_limits, throw_limits_up;
      for (auto it = fSystParameters.begin();
           it != fSystParameters.end(); ++it) {
        throw_limits.push_back(it->second.GetGenThrowLimit());
        throw_limits_up.push_back(it->second.GetGenThrowLimitUp());
      }


      bool rethrow = true;
      while (rethrow) {
        bool all_pos = true;
        TVectorD rand(fTotalSystParameters);
        for (int i = 0; i < (int)fTotalSystParameters; ++i) {
          rand[fCovarianceBinsSimple[i]] = fRNG.Gaus(); 
          std::cout << "Throwing " << i << " " << rand[i] << std::endl;
        }
        TVectorD rand_times_chol = (*cov_lower)*rand;

        for (size_t i = 0; i < fTotalSystParameters; ++i) {
          double val = rand_times_chol[fCovarianceBinsSimple[i]] +
                       nominal_vals[nonsyst_pars + i];
          if (val < throw_limits[i]) {
            all_pos = false;
            std::cout << "Rethrowing " << i << " " << val << " " <<
                         throw_limits[i] <<
                         std::endl;
          }
          if (val > throw_limits_up[i]) {
            all_pos = false;
            std::cout << "Rethrowing " << i << " " << val << " " <<
                         throw_limits_up[i] <<
                         std::endl;
          }
          vals[nonsyst_pars + i] = val;
          std::cout << "Set " << nonsyst_pars + i << " " << vals[nonsyst_pars + i] << std::endl;
        }

        rethrow = !all_pos;
      }

      std::cout << "DOne throwing" << std::endl;
      for (size_t i = 0; i < vals.size(); ++i) {
        std::cout << i << " " << vals[i] << std::endl;
      }
    }

    fFillIncidentInFunction = true;
    fUseFakeSamples = true;
    fThinSliceDriver->SetFillFakeInMain(true);

    std::cout << "Vary stats: " << fVaryMCStatsForFakeData << std::endl;
    std::cout << fDataBeamFluxes.size() << std::endl;
    fThinSliceDriver->SetStatVar(fVaryMCStatsForFakeData); //Force on stat vars in driver
    fThinSliceDriver->TurnOnFakeData();

    std::cout << fDataBeamFluxes.size() << std::endl;
    fFitFunction(&vals[0]);
    fDataSet.FillHistsFromSamples(fFakeSamples, fDataFlux, fDataBeamFluxes,
                                  (fDoFluctuateStats && fFluctuateInSamples));
    std::cout << "Fake Data fluxes: " << std::endl;
    for (const auto & f: fDataBeamFluxes) std::cout << f << " ";
    std::cout << std::endl;
    fThinSliceDriver->TurnOffFakeData();
    BuildFakeDataXSecs(false);

    fUseFakeSamples = false;
    fThinSliceDriver->SetFillFakeInMain(false);
    fThinSliceDriver->SetStatVar(false); //Force off stat vars in driver
    //Refill the hists for comparisons
    fFitFunction(&nominal_vals[0]); 
    fFillIncidentInFunction = false;
  }
  else {
    if (!fNoFakeReset) {
      for (auto it = fFakeSamples.begin(); it != fFakeSamples.end(); ++it) {
        for (size_t i = 0; i < it->second.size(); ++i) {
          for (size_t j = 0; j < it->second[i].size(); ++j) {
            it->second[i][j].Reset();
          }
        }
      }
    }
    for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
      std::vector<ThinSliceSample> & samples_vec = it->second[0];
      fFakeDataScales[it->first] = std::vector<double>(samples_vec.size(), 0.);
    }
    fThinSliceDriver->BuildFakeData(tree, fFakeDataEvents, fFakeSamples,
                                    fIsSignalSample, fDataSet,
                                    fDataFlux, fFakeDataScales, fBeamEnergyBins,
                                    fDataBeamFluxes,
                                    fSplitVal, fScaleToDataBeamProfile);

    BuildFakeDataXSecs(false);

    std::cout << std::endl << " Data by beam bin: ";
    for (const auto & f : fDataBeamFluxes) {
      std::cout << f << " ";
    }

    if (fDoFluctuateStats && fFluctuateInSamples) {
      fDataSet.FillHistsFromSamples(fFakeSamples, fDataFlux, fDataBeamFluxes,
                                    (fDoFluctuateStats && fFluctuateInSamples));
    }

    std::cout << std::endl << " Data by beam bin: ";
    for (const auto & f : fDataBeamFluxes) {
      std::cout << f << " ";
    }
    std::cout << std::endl;
  }

  fDataSet.MakeRebinnedHists();
  fDataSet.SetDirectory();

  inputFile->Close();
}

void protoana::PDSPThinSliceFitter::SaveDataSet() {
  fOutputFile.cd();
  TDirectory * out = (TDirectory *)fOutputFile.mkdir("Data");
  out->cd();
  TVectorD out_fluxes(fDataBeamFluxes.size());
  for (size_t i = 0; i < fDataBeamFluxes.size(); ++i) {
    out_fluxes[i] = fDataBeamFluxes[i];
  }
  out_fluxes.Write("data_fluxes");
  std::map<int, TH1 *> & selected_hists = fDataSet.GetSelectionHists();
  for (auto it = selected_hists.begin(); it != selected_hists.end(); ++it) {
    it->second->Write();
  }
}


void protoana::PDSPThinSliceFitter::SaveMCSamples() {
  fOutputFile.cd();
  fOutputFile.mkdir("MC_Samples");
  fOutputFile.cd("MC_Samples");

  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      auto vec = it->second.at(i);
      for (size_t j = 0; j < vec.size(); ++j) {
        //const std::map<int, TH1 *> & hists = vec[j].GetSelectionHists();
        const auto & hists = vec[j].GetSelectionHists();
        for (auto it2 = hists.begin(); it2 != hists.end(); ++it2) {
          it2->second->Write();
        }
      }
    }
  }


  fOutputFile.cd();
}


double protoana::PDSPThinSliceFitter::GetNMuons() {
  double total = 0.;
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        if (it->second[i][j].GetFluxType() == 2) {
          total += it->second[i][j].GetVariedFlux();
        }
      }
    }
  }
  return total;

}

void protoana::PDSPThinSliceFitter::CompareDataMC(
    std::string extra_name, TDirectory * xsec_dir, TDirectory * plot_dir,
    bool post_fit) {
      fThinSliceStrategy->CompareDataMC(fMCDists, fDataSet, fOutputFile);
      fThinSliceStrategy->CompareSelections(fMCDists, fDataSet, fOutputFile, fPlotStyle, false, post_fit, plot_dir);
  // fThinSliceDriver->CompareDataMC(
  //     fSamples, fDataSet, fOutputFile,
  //     fPlotStyle, (fTotalFluxParameters + fTotalSignalParameters +
  //                  fTotalSystParameters + fTotalG4RWParameters),
  //     plot_dir, fPlotRebinned, post_fit);

  // xsec_dir->cd();

  // //Get the incident histogram from all of the relevant samples
  // for (size_t i = 0; i < fMeasurementSamples.size(); ++i) {
  //   int sample_ID = fMeasurementSamples[i];
  //   auto & samples_vec_2D
  //       = fSamples[sample_ID]; 

  //   std::vector<double> signal_bins = fSignalBins[sample_ID];
  //   if (fDrawXSecUnderflow) {
  //     signal_bins.insert(signal_bins.begin(), 0.);
  //     signal_bins.push_back(2000.);//Make this like inf or something big
  //   }

  //   std::string signal_name = extra_name;
  //   signal_name += samples_vec_2D[0][1].GetName() + "Hist";
  //   TH1D signal_hist(signal_name.c_str(), "", signal_bins.size() - 1,
  //                    &signal_bins[0]);

  //   for (size_t j = 0; j < samples_vec_2D.size(); ++j) {
  //     auto & samples_vec = samples_vec_2D[j];

  //     size_t start = (fDrawXSecUnderflow ? 0 : 1);
  //     size_t end = samples_vec.size() - (fDrawXSecUnderflow ? 0 : 1);
  //     for (size_t k = start; k < end; ++k) {
  //     //for (size_t k = (fDrawXSecUnderflow ? 0 : 1);
  //     //     k < samples_vec.size()-1; ++k) {
  //       ThinSliceSample & sample = samples_vec[k];

  //       int the_bin = (fDrawXSecUnderflow ? k+1 : k);
  //       signal_hist.AddBinContent(the_bin, sample.GetVariedFlux());

  //     }
  //   }
  //   signal_hist.Write();

  //   //Get the incident histogram from all of the relevant samples
  //   std::string inc_name = extra_name;
  //   inc_name += "TotalIncident" + samples_vec_2D[0][0].GetName();

  //   TH1D * total_incident_hist = new TH1D(inc_name.c_str(), "",
  //                                         signal_bins.size() - 1,
  //                                         &signal_bins[0]);
  //   std::map<int, std::vector<TH1D*>> temp_hists;
  //   for (size_t i = 0; i < fIncidentSamples.size(); ++i) {
  //     auto & vec_2D = fSamples[fIncidentSamples[i]];
  //     temp_hists[fIncidentSamples[i]] = std::vector<TH1D*>();
  //     for (size_t j = 0; j < vec_2D.size(); ++j) {
  //       auto & samples_vec = vec_2D[j];
  //       for (size_t k = 0; k < samples_vec.size(); ++k) {
  //         if (j == 0) {
  //           std::string name = extra_name;
  //           name += "Incident" + samples_vec_2D[0][1].GetName() +
  //                    samples_vec[k].GetName() + std::to_string(k);
  //           std::string title = samples_vec[k].GetSelectionHists().begin()->second->GetTitle();
  //           temp_hists[fIncidentSamples[i]].push_back(
  //               new TH1D(name.c_str(),
  //                        title.c_str(),
  //                        signal_bins.size() - 1,
  //                        &signal_bins[0]));
  //         }
  //         samples_vec[k].FillHistFromIncidentEnergies(*total_incident_hist);
  //         samples_vec[k].FillHistFromIncidentEnergies(
  //             *(temp_hists[fIncidentSamples[i]][k]));
  //       }
  //     }
  //   }
  //   total_incident_hist->Write();
  //   if (post_fit) {
  //     fBestFitIncs[sample_ID] = total_incident_hist;
  //   }
  //   else {
  //     fNominalIncs[sample_ID] = total_incident_hist;
  //   }

  //   std::string xsec_name = extra_name + samples_vec_2D[0][1].GetName() +
  //                           "XSec";
  //   TH1D * xsec_hist = (TH1D*)signal_hist.Clone(xsec_name.c_str());
  //   xsec_hist->Sumw2();
  //   xsec_hist->Divide(total_incident_hist);
  //   CalculateCrossSection(xsec_hist);

  //   xsec_hist->Write();
  //   if (post_fit) {
  //     fBestFitXSecs[sample_ID] = xsec_hist;
  //   }
  //   else {
  //     fNominalXSecs[sample_ID] = xsec_hist;
  //   }

  //   std::string stack_name = extra_name;
  //   stack_name += "IncidentStack" + samples_vec_2D[0][1].GetName();
  //   THStack inc_stack(stack_name.c_str(), ";Reconstructed KE (MeV)");
  //   int iColor = 0;
  //   for (auto it = temp_hists.begin(); it != temp_hists.end(); ++it) {
  //     for (size_t i = 0; i < it->second.size(); ++i) {
  //       std::pair<int, int> color_fill =
  //           fThinSliceDriver->GetColorAndStyle(iColor, fPlotStyle);
  //       it->second[i]->SetFillColor(color_fill.first);
  //       it->second[i]->SetFillStyle(color_fill.second);
  //       it->second[i]->SetLineColor(1);
  //       inc_stack.Add(it->second[i]);
  //       ++iColor;
  //     }
  //   }
  //   inc_stack.Write();
  // }

  // fOutputFile.cd(); 
  // TVector muons(1);
  // muons[0] = GetNMuons();

  // auto muon_names = extra_name + "Muons";
  // muons.Write(muon_names.c_str());


  // //Get the Best fit selection hists
  // std::map<int, std::vector<TH1*>> temp_map;
  // fThinSliceDriver->GetCurrentHists(fSamples, fDataSet, temp_map,
  //                                   fPlotRebinned);
  // for (auto it = temp_map.begin(); it != temp_map.end(); ++it) {
  //   fBestFitSelectionHists[it->first] = (TH1D*)it->second[0];
  // }

  // //Get the best fit truth vals
  // for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
  //   fBestFitTruthVals[it->first]
  //       = std::vector<double>(it->second[0].size(), 0.);
   
  //   for (size_t i = 0; i < it->second[0].size(); ++i) {
  //     double best_fit_val_i = 0.;
  //     for (size_t j = 0; j < it->second.size(); ++j) {
  //       best_fit_val_i += it->second[j][i].GetVariedFlux();
  //     }

  //     fBestFitTruthVals[it->first][i] = best_fit_val_i;
  //   }
  // }
}

void protoana::PDSPThinSliceFitter::AdjustInitConds(size_t fit_attempt) {
  //Force Random start
  //if (fit_attempt > 4) {
    fRandomStart = true;
  //}

  size_t remainder = fit_attempt % 4;
  switch (remainder) {
    case 1: 
      fMinimizer->SetTolerance(10.*fTolerance);
      break;
    case 2:
      fMinimizer->SetTolerance(1.e-1*fTolerance);
      break;
    case 3: 
      fMinimizer->SetTolerance(100.*fTolerance);
      break;
    default:
      fMinimizer->SetTolerance(1.e-2*fTolerance);
  }
}

void protoana::PDSPThinSliceFitter::NormalFit() {
  std::cout << "Running Fit" << std::endl;
  bool found_minimum = false;

  SetupTree();

  size_t n_fit_attempts = 0;
  size_t max_attempts = (fRerunFit ? fFitAttempts : 1);

  TVector output_fit_status(max_attempts);
  TVector output_cov_status(max_attempts);
  TVector output_hesse_cov_status(max_attempts);

  while (n_fit_attempts < max_attempts) {
    if (n_fit_attempts > 0) {
      AdjustInitConds(n_fit_attempts);
      fOutputTree->Reset();
      fMinimizer->Clear();
      ResetParameters();
    }
    std::cout << "Fit attempt " << n_fit_attempts << std::endl;
    try {
      found_minimum = fMinimizer->Minimize();
    }
    catch (const std::exception & e) {
      std::cerr << e.what() << std::endl;
      std::cerr << "exiting safely" << std::endl;
      found_minimum = 0;
    }

    output_fit_status[n_fit_attempts] = fMinimizer->Status();
    output_cov_status[n_fit_attempts] = fMinimizer->CovMatrixStatus();

    if (found_minimum && fRunHesse) {
      found_minimum = DoHesse();
      if (fRequireGoodHesse)
        found_minimum = (fMinimizer->CovMatrixStatus() == 3);
      output_hesse_cov_status[n_fit_attempts] = fMinimizer->CovMatrixStatus();
    }

    if (found_minimum) break;
    ++n_fit_attempts;

  }

  WrapUpTree();

  if (!found_minimum) {
    std::cout << "Failed to find minimum" << std::endl;
    fOutputFile.cd();
    output_fit_status.Write("fit_status");
    output_cov_status.Write("cov_status");
    output_hesse_cov_status.Write("post_hesse_cov_status");

    if (fMinimizer->Status() == 3) {
      std::cout << "Edm too large: " << fMinimizer->Edm() << std::endl;
    }

    if (fDoScans) ParameterScans();
  }
  else {
    std::vector<double> vals;
    std::cout << "Found minimimum. Fit status: " << fMinimizer->Status() << std::endl;
    double chi2_syst = (fAddSystTerm ? CalcChi2SystTerm() : 0.);
    double chi2_stat = fThinSliceDriver->CalculateChi2(fSamples, fDataSet).first;
    double chi2_reg = (fAddRegTerm ? CalcRegTerm() : 0.);
    double chi2_sig = (fAddSignalConstraint ? CalcSignalConstraint() : 0.);
    std::cout << "chi2 Syst: " << chi2_syst << std::endl;
    std::cout << "chi2 Stat: " << chi2_stat << std::endl;
    std::cout << "chi2 Reg: " << chi2_reg << std::endl;
    std::cout << "chi2 Sig: " << chi2_sig << std::endl;
    output_cov_status.Write("cov_status");

    GetCovarianceVals("Migrad");

    //if (fRunHesse) SaveHesse();
    output_hesse_cov_status.Write("post_hesse_cov_status");
    if (fRunMinos1D) DoMinos1D();
    if (fRunMinosConts) DoMinosConts();
    if (fRunMultiConts) DoMultiConts();


    // fOutputFile.cd();
    // auto * dir = fOutputFile.mkdir("ExtraHistsMCPostFit");

    // SaveExtraHists(dir);
    fOutputFile.cd();

    size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters + fTotalG4RWParameters;
    std::cout << total_parameters << std::endl;
    for (size_t i = 0; i < total_parameters; ++i) {
      std::cout << fMinimizer->VariableName(i) << " " << fMinimizer->X()[i] <<
                   std::endl;
      vals.push_back(fMinimizer->X()[i]);
    }

    TH2D covHist("covHist", "", total_parameters, 0,
                 total_parameters, total_parameters, 0,
                 total_parameters);
    TH2D corrHist("corrHist", "", total_parameters, 0,
                  total_parameters, total_parameters, 0,
                  total_parameters);

    TH1D parsHist("postFitPars", "", total_parameters, 0,
                  total_parameters);
    TH1D parsHist_normal_val("postFitParsNormal", "", 
                             total_parameters, 0, total_parameters);
    TH1D * toyParsHist = 0x0,
         * toyParsHist_normal = 0x0;
    if (fFakeDataRoutine == "Toy") {
      toyParsHist = (TH1D*)parsHist.Clone("toyPars");
      toyParsHist->Reset();
      toyParsHist_normal = (TH1D*)parsHist.Clone("toyPars_normal");
      toyParsHist_normal->Reset();
    }



    int n_par = 0;
    for (auto it = fSignalParameters.begin();
         it != fSignalParameters.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
        parsHist.SetBinError(n_par+1,
                             sqrt(fMinimizer->CovMatrix(n_par, n_par)));
        parsHist.GetXaxis()->SetBinLabel(
            n_par+1, fSignalParameterNames[it->first][i].c_str());

        parsHist_normal_val.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
        parsHist_normal_val.SetBinError(
            n_par+1, sqrt(fMinimizer->CovMatrix(n_par, n_par)));
        parsHist_normal_val.GetXaxis()->SetBinLabel(
            n_par+1, fSignalParameterNames[it->first][i].c_str());

        if (fFakeDataRoutine == "Toy") 
          toyParsHist->SetBinContent(n_par+1, -999.);


        ++n_par;
      }
    }
    for (auto it = fFluxParameters.begin();
         it != fFluxParameters.end(); ++it) {
      parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
      parsHist.SetBinError(n_par+1,
                           sqrt(fMinimizer->CovMatrix(n_par, n_par)));
      parsHist.GetXaxis()->SetBinLabel(
          n_par+1, fFluxParameterNames[it->first].c_str());

      parsHist_normal_val.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
      parsHist_normal_val.SetBinError(
          n_par+1, sqrt(fMinimizer->CovMatrix(n_par, n_par)));
      parsHist_normal_val.GetXaxis()->SetBinLabel(
          n_par+1, fFluxParameterNames[it->first].c_str());


      if (fFakeDataRoutine == "Toy") 
        toyParsHist->SetBinContent(n_par+1, -999.);


      ++n_par;
    }

    int cov_bin_counter = 0;
    for (auto it = fSystParameters.begin();
         it != fSystParameters.end(); ++it) {

      if (abs(it->second.GetCentral()) < 1.e-5) {
        parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par] + 1.);
        parsHist.SetBinError(n_par+1,
                             sqrt(fMinimizer->CovMatrix(n_par, n_par)));
        if (fFakeDataRoutine == "Toy") 
          toyParsHist->SetBinContent(n_par+1, fToyValues[it->first] + 1.);

      }
      else {
        parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par]/it->second.GetCentral());
        parsHist.SetBinError(n_par+1,
                             sqrt(fMinimizer->CovMatrix(n_par, n_par))/it->second.GetCentral());
        if (fFakeDataRoutine == "Toy")
          toyParsHist->SetBinContent(n_par+1, fToyValues[it->first]/it->second.GetCentral());
      }

      if (fFakeDataRoutine == "Toy") {
        toyParsHist_normal->SetBinContent(n_par+1, fToyValues[it->first]);
        //size_t cov_bin = fCovarianceBins[it->first];
        size_t cov_bin = fCovarianceBinsSimple[cov_bin_counter];
        ++cov_bin_counter;

        toyParsHist_normal->SetBinError(
            n_par+1, sqrt((*fCovMatrixDisplay)[cov_bin][cov_bin]));
      }

      parsHist.GetXaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());

      parsHist_normal_val.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
      parsHist_normal_val.SetBinError(
          n_par+1, sqrt(fMinimizer->CovMatrix(n_par, n_par)));
      parsHist_normal_val.GetXaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());


      ++n_par;
    }


    for (auto it = fG4RWParameters.begin();
         it != fG4RWParameters.end(); ++it) {

      if (abs(it->second.GetCentral()) < 1.e-5) {
        parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par] + 1.);
        parsHist.SetBinError(n_par+1,
                             sqrt(fMinimizer->CovMatrix(n_par, n_par)));
      }
      else {
        parsHist.SetBinContent(n_par+1, fMinimizer->X()[n_par]/it->second.GetCentral());
        parsHist.SetBinError(n_par+1,
                             sqrt(fMinimizer->CovMatrix(n_par, n_par))/it->second.GetCentral());
      }

      parsHist.GetXaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());

      parsHist_normal_val.SetBinContent(n_par+1, fMinimizer->X()[n_par]);
      parsHist_normal_val.SetBinError(
          n_par+1, sqrt(fMinimizer->CovMatrix(n_par, n_par)));
      parsHist_normal_val.GetXaxis()->SetBinLabel(
          n_par+1, it->second.GetName().c_str());


      ++n_par;
    }

    TMatrixD * cov = new TMatrixD(total_parameters, total_parameters);

    //std::cout << "First cov: " << std::endl;
    for (size_t i = 0; i < total_parameters; ++i) {
      for (size_t j = 0; j < total_parameters; ++j) {
        double cov_val = fMinimizer->CovMatrix(i, j);
        covHist.SetBinContent(i+1, j+1, cov_val);
        corrHist.SetBinContent(i+1, j+1, fMinimizer->Correlation(i, j));
        //std::cout << fMinimizer->Correlation(i, j) << " ";
        (*cov)(i, j) = fMinimizer->CovMatrix(i, j);
      }
      //std::cout << std::endl;
    }

    if (fFixVariables && fAddSystTerm) {
      size_t n_par = fTotalSignalParameters + fTotalFluxParameters;
      int cov_bin_counter = 0;
      for (auto it = fSystParameters.begin(); it != fSystParameters.end();
           ++it) {
        if (fSystsToFix.find(it->first) != fSystsToFix.end()) {
          //size_t cov_bin = fCovarianceBins[it->first];
          size_t cov_bin = fCovarianceBinsSimple[cov_bin_counter];
          covHist.SetBinContent(n_par+1, n_par+1,
                                (*fCovMatrixDisplay)[cov_bin][cov_bin]);
          (*cov)(n_par, n_par) = (*fCovMatrixDisplay)[cov_bin][cov_bin];
        }
        ++n_par;
        ++cov_bin_counter;
      }

    }

    fOutputFile.cd();
    output_fit_status.Write("fit_status");
    parsHist.SetFillColor(kRed);
    parsHist.SetFillStyle(3144);
    parsHist.SetMarkerColor(kRed);
    parsHist.SetMarkerStyle(20);
    parsHist.Write();
    parsHist_normal_val.Write();

    TVectorD chi2_syst_out(1);
    chi2_syst_out[0] = chi2_syst;
    chi2_syst_out.Write("chi2_syst");
    TVectorD chi2_stat_out(1);
    chi2_stat_out[0] = chi2_stat;
    chi2_stat_out.Write("chi2_stat");
    TVectorD chi2_reg_out(1);
    chi2_reg_out[0] = chi2_reg;
    chi2_reg_out.Write("chi2_reg");
    TVectorD chi2_sig_out(1);
    chi2_sig_out[0] = chi2_sig;
    chi2_sig_out.Write("chi2_sig");

    fBestFitSignalPars = fSignalParameters; 
    fBestFitFluxPars = fFluxParameters; 
    fBestFitSystPars = fSystParameters;
    fBestFitG4RWPars = fG4RWParameters;

    //Drawing pre + post fit pars
    TCanvas cPars("cParameters", "");
    cPars.SetTicks();
    parsHist.GetYaxis()->SetTitle("Parameter Values");
    parsHist.SetTitleSize(.05, "Y");
    //parsHist.Draw("e2");
    TH1D * preFitParsHist = (TH1D*)fOutputFile.Get("preFitPars");
    preFitParsHist->SetTitle(";Parameter;Parameter Value");
    preFitParsHist->GetXaxis()->SetTitleOffset(.8);
    preFitParsHist->GetYaxis()->SetTitleOffset(.8);
    preFitParsHist->SetTitleSize(.05, "XY");
    preFitParsHist->SetFillColor(kBlue);
    preFitParsHist->SetFillStyle(3001);
    //preFitParsHist->Draw("pe2 same");
    preFitParsHist->SetMaximum(3.);
    preFitParsHist->SetMinimum(-1.);
    preFitParsHist->Draw("pe2");
    parsHist.Draw("e2 same");
    TLine l(0., 1., parsHist.GetXaxis()->GetBinUpEdge(parsHist.GetNbinsX()), 1.);
    l.SetLineColor(kBlack);
    l.Draw("same");
    TLegend leg(.15, .65, .45, .85);
    leg.AddEntry(preFitParsHist, "Pre-Fit #pm 1 #sigma", "pf");
    leg.AddEntry(&parsHist, "Post-Fit #pm 1 #sigma", "pf");

    if (fFakeDataRoutine == "Toy") {
      toyParsHist->SetMarkerStyle(20);
      toyParsHist->SetMarkerColor(kBlack);
      toyParsHist->Draw("same p");
      leg.AddEntry(toyParsHist, "Toy Value", "p");
      toyParsHist->Write();
      toyParsHist_normal->Write();
    }
    preFitParsHist->SetLineWidth(0);
    preFitParsHist->Draw("p same");
    parsHist.SetLineWidth(0);
    parsHist.Draw("p same");


    leg.Draw();
    cPars.Write();

    covHist.SetTitle(";Parameter;Parameter");
    covHist.SetTitleSize(.05, "XY");
    covHist.GetXaxis()->SetTitleOffset(.8);
    covHist.GetYaxis()->SetTitleOffset(.8);
    covHist.Write();
    cov->Write("CovMatrix");
    corrHist.SetTitle(";Parameter;Parameter");
    corrHist.SetTitleSize(.05, "XY");
    corrHist.GetXaxis()->SetTitleOffset(.8);
    corrHist.GetYaxis()->SetTitleOffset(.8);
    corrHist.SetMaximum(1.);
    corrHist.SetMinimum(-1.);
    corrHist.Write();

    TCanvas cCorr("cCorr", "");
    corrHist.Draw("colz");
    cCorr.Write();

    if (fDoScans)
      ParameterScans();

    fFillIncidentInFunction = true;
    fFitFunction(&vals[0]);

    //SetBestFit();  //Deprecated

    //save post fit stacks
    TDirectory * top_dir = fOutputFile.GetDirectory("");
    TDirectory * xsec_dir = fOutputFile.mkdir("PostFitXSec");
    std::string extra_name = "PostFit";
    CompareDataMC(extra_name, xsec_dir, top_dir, true);


    //Make 'fixed' 
    if (fFixPostFit) {
      GetFixFactors();
      fFixSamplesInFunction = true;
      fFitFunction(&vals[0]);
      TDirectory * fixed_plot_dir = fOutputFile.mkdir("FixedPlots");
      TDirectory * fixed_xsec_dir = fixed_plot_dir->mkdir("FixedXSec");
      extra_name = "Fixed";
      CompareDataMC(extra_name, fixed_xsec_dir, fixed_plot_dir, true);
      fFixSamplesInFunction = false;
    }


    if (fDoThrows) 
      DoThrows(parsHist_normal_val, cov);

    if (fDo1DShifts) {
      if (fAddSystTerm)
        Do1DShifts(fPreFitParsNormal, true);
      Do1DShifts(parsHist_normal_val);
    }
  }
}


void protoana::PDSPThinSliceFitter::RunFitAndSave() {
  DefineFitFunction();
  fMinimizer->SetFunction(fFitFunction);

  fThinSliceDriver->SetUseMCStatVarWeight(fUseMCStatVarWeightFakeData);
  BuildDataHists();
  if (fDoFluctuateStats && !fFluctuateInSamples) {
    //fDataSet.GetCumulatives();
    fDataSet.GenerateStatFluctuation(fDataBeamFluxes);
  }
  SaveDataSet();


  fThinSliceDriver->SetUseMCStatVarWeight(fUseMCStatVarWeight);
  fThinSliceDriver->SetStatVar(fVaryMCStats); //Force on stat vars in driver
  if (fSetValsPreFit) {
    size_t par_position = 0;
    for (auto it = fSignalParameters.begin();
         it != fSignalParameters.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        it->second.at(i) = fPreFitVals[par_position];
        ++par_position;
      }
    }
    for (auto it = fFluxParameters.begin();
         it != fFluxParameters.end(); ++it) {
      it->second = fPreFitVals[par_position];
      ++par_position;
    }

    for (auto it = fSystParameters.begin();
         it != fSystParameters.end(); ++it) {
      it->second.SetValue(fPreFitVals[par_position]);
      //std::cout << it->first << " " << it->second.GetValue() << std::endl;
      ++par_position;
    }


    for (auto it = fG4RWParameters.begin();
         it != fG4RWParameters.end(); ++it) {
      it->second.SetValue(fPreFitVals[par_position]);
      //std::cout << it->first << " " << it->second.GetValue() << std::endl;
      ++par_position;
    }

    std::cout << fPreFitVals.size() << std::endl;
    for (auto p : fPreFitVals) {
      std::cout << p << std::endl;
    }
    fFillIncidentInFunction = true;
    fFitFunction(&fPreFitVals[0]);
    fFillIncidentInFunction = false;
  }
  else {
    std::vector<double> vals;
    for (auto it = fSignalParameters.begin();
         it != fSignalParameters.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        vals.push_back(it->second.at(i));
      }
    }
    for (auto it = fFluxParameters.begin();
         it != fFluxParameters.end(); ++it) {
      vals.push_back(it->second);
    }

    for (auto it = fSystParameters.begin();
         it != fSystParameters.end(); ++it) {
      std::cout << "syst: " << it->second.GetValue() << std::endl;
      vals.push_back(it->second.GetValue());
    }

    for (auto it = fG4RWParameters.begin();
         it != fG4RWParameters.end(); ++it) {
      vals.push_back(it->second.GetValue());
    }

    std::cout << vals.size() << std::endl;
    for (auto p : vals) {
      std::cout << p << std::endl;
    }

    fFillIncidentInFunction = true;
    fFitFunction(&vals[0]);
    fFillIncidentInFunction = false;
  }

  SaveMCSamples();


  TDirectory * top_dir = fOutputFile.GetDirectory("");
  TDirectory * xsec_dir = fOutputFile.mkdir("PreFitXSec");
  std::string extra_name = "PreFit";
  CompareDataMC(extra_name, xsec_dir, top_dir);

  if (fFitType == "Normal") {
    NormalFit();
  }

  else if (fFitType == "None") {
    if (fDoScans) ParameterScans();
    if (fDo1DShifts) {
      if (fAddSystTerm)
        Do1DShifts(fPreFitParsNormal, true);
    }
  }
  else if (fFitType == "ThrowsOnly") {
    //Load up parsHist and covariance
    TFile * prior_fit_file = TFile::Open(fFitResultsFile.c_str());
    TH1D * parsHist_normal_val = (TH1D*)prior_fit_file->Get("postFitParsNormal");
    TMatrixD * cov = (TMatrixD*)prior_fit_file->Get("CovMatrix");

    //Also load up best fit xsecs
    for (auto & s : fMeasurementSamples) {
      auto & samples_vec_2D = fSamples[s]; 

      std::string xsec_name = "PostFitXSec/PostFit" + samples_vec_2D[0][1].GetName() +
                              "XSec";
      fBestFitXSecs[s]
          = (TH1D*)prior_fit_file->Get(xsec_name.c_str());

      std::string inc_name = "PostFitXSec/PostFitTotalIncident" +
                             samples_vec_2D[0][0].GetName();
      fBestFitIncs[s]
          = (TH1D*)prior_fit_file->Get(inc_name.c_str());
      std::cout << fBestFitXSecs[s] << " " << fBestFitIncs[s] << std::endl;

      //Nominals
      xsec_name = "PreFitXSec/PreFit" + samples_vec_2D[0][1].GetName() +
                              "XSec";
      fNominalXSecs[s]
          = (TH1D*)prior_fit_file->Get(xsec_name.c_str());
    }

    DoThrows(*parsHist_normal_val, cov);

    //Copy some hists/canvases/graphs to new output
    TDirectory * postfit_input = (TDirectory*)prior_fit_file->Get("PostFitXSec");
    TList * postfit_keys = postfit_input->GetListOfKeys();
    TList * top_keys = prior_fit_file->GetListOfKeys();

    fOutputFile.cd();
    TDirectory * postfit_dir = (TDirectory*)fOutputFile.mkdir("PostFitXSec");
    postfit_dir->cd();
    for (int i = 0; i < postfit_keys->GetSize(); ++i) {
      postfit_input->Get(postfit_keys->At(i)->GetName())->Write();
    }

    std::vector<std::string> to_copy = {
      "cCorr", "covHist", "corrHist", "CovMatrix", "chi2_syst",
      "chi2_stat", "chi2_reg", "chi2_sig"
    };

    fOutputFile.cd();
    for (int i = 0; i < top_keys->GetSize(); ++i) {
      std::string name = top_keys->At(i)->GetName();

      if (name.find("ostFit") != std::string::npos) {
        std::cout << name << " " << (name.find("ostFit") != std::string::npos) << std::endl;
        prior_fit_file->Get(name.c_str())->Write(name.c_str());
      }

      if (std::find(to_copy.begin(), to_copy.end(), name) != to_copy.end()) {
        std::cout << name << std::endl;
        prior_fit_file->Get(name.c_str())->Write(name.c_str());
      }
    }
  }
  else {
    std::string message = "Error: invalid fit type given -- ";
    message += fFitType;
    throw std::runtime_error(message);
  }


}

void protoana::PDSPThinSliceFitter::BuildDataFromToy() {
  if (!fAddSystTerm) {
    std::string message
        = "Error: Need to include systematic term and covariance matrix";
    throw std::runtime_error(message);
  }

  std::vector<double> vals, init_vals;
  for (auto it = fSignalParameters.begin();
       it != fSignalParameters.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      vals.push_back(it->second.at(i));
      init_vals.push_back(it->second.at(i));
    }
  }
  for (auto it = fFluxParameters.begin();
       it != fFluxParameters.end(); ++it) {
    vals.push_back(it->second);
    init_vals.push_back(it->second);
  }

  std::vector<size_t> bins;
  std::vector<std::string> names;
  std::vector<double> throw_limits, throw_limits_up;
  for (auto it = fSystParameters.begin();
       it != fSystParameters.end(); ++it) {
    bins.push_back(fCovarianceBins[it->first]);
    names.push_back(it->first);
    init_vals.push_back(it->second.GetValue());
    throw_limits.push_back(it->second.GetGenThrowLimit());
    throw_limits_up.push_back(it->second.GetGenThrowLimitUp());
  }

  size_t base_par = fTotalSignalParameters + fTotalFluxParameters;

  if (fAnalysisOptions.get<bool>("SetToy", false)) {
    std::vector<double> toy_vals
        = fAnalysisOptions.get<std::vector<double>>("ToyVals");
    for (double & v : toy_vals) {
      vals.push_back(v);
    }

    for (size_t i = 0; i < vals.size(); ++i) {
      std::cout << vals[i] << " ";
      if (i >= base_par) {
        std::cout << names[i - base_par];
        fToyValues[names[i - base_par]] = vals[i];
      }
      std::cout << std::endl;
    }
  }
  else {

    TMatrixD * cov_lower = (TMatrixD*)fInputChol->GetU().Clone();
    cov_lower->Transpose(fInputChol->GetU());

    size_t n_cov_rows = cov_lower->GetNrows();
    if (n_cov_rows != fTotalSystParameters) {
      std::cout << "Error: " << n_cov_rows << " " << fTotalSystParameters <<
                   std::endl;
    }

    bool rethrow = true;
    while (rethrow) {
      bool all_pos = true;
      TVectorD rand(fTotalSystParameters);
      for (size_t i = 0; i < fTotalSystParameters; ++i) {
        rand[i] = fRNG.Gaus();
      }
      //TVectorD rand_times_chol = fInputChol->GetU()*rand;
      TVectorD rand_times_chol = (*cov_lower)*rand;

      for (size_t i = 0; i < fTotalSystParameters; ++i) {
        //double val = rand_times_chol[bins[i]] + init_vals[base_par + i];
        double val = rand_times_chol[fCovarianceBinsSimple[i]] +
                     init_vals[base_par + i];
        //if (val < fSystParameters[names[i]].GetGenThrowLimit()) {
        if (val < throw_limits[i]) {
          all_pos = false;
          std::cout << "Rethrowing " << i << " " << val << " " <<
                       throw_limits[i] <<
                       std::endl;
        }
        //if (val > fSystParameters[names[i]].GetGenThrowLimitUp()) {
        if (val > throw_limits_up[i]) {
          all_pos = false;
          std::cout << "Rethrowing " << i << " " << val << " " <<
                       throw_limits_up[i] <<
                       std::endl;
        }

      }

      if (all_pos) {
        for (size_t i = 0; i < fTotalSystParameters; ++i) {
          if (fSystsToFix.find(names[i]) != fSystsToFix.end()) {

          }
          double val = (fSystsToFix.find(names[i]) == fSystsToFix.end() ?
                        //rand_times_chol[bins[i]] + init_vals[base_par + i] :
                        (rand_times_chol[fCovarianceBinsSimple[i]] +
                         init_vals[base_par + i]) : fSystsToFix[names[i]]);

          vals.push_back(val);
          fToyValues[names[i]] = val;
        }
      }
      rethrow = !all_pos;
    }

    std::cout << "Vals: " << vals.size() << " total " <<
                 base_par +
                 fTotalSystParameters << std::endl;
    for (size_t i = 0; i < vals.size(); ++i) {
      std::cout << vals[i] << " ";
      if (i >= base_par)
        std::cout << names[i - base_par];
      std::cout << std::endl;
    }
  }

  for (auto it = fG4RWParameters.begin();
       it != fG4RWParameters.end(); ++it) {
    vals.push_back(it->second.GetValue());
  }

  fFillIncidentInFunction = true;
  fUseFakeSamples = true;
  fThinSliceDriver->SetStatVar(fVaryMCStatsForFakeData); //Force on stat vars in driver
  fFitFunction(&vals[0]);
  fDataSet.FillHistsFromSamples(fFakeSamples, fDataFlux, fDataBeamFluxes,
                                (fDoFluctuateStats && fFluctuateInSamples));
  BuildFakeDataXSecs(false);

  fUseFakeSamples = false;
  fThinSliceDriver->SetStatVar(false); //Turn it back off
  //Refill the hists for comparisons
  fFitFunction(&init_vals[0]);
  fFillIncidentInFunction = false;
  

}

std::vector<double> protoana::PDSPThinSliceFitter::GetBestFitParsVec() {
  std::vector<double> results;

  for (auto it = fSignalParameters.begin();
       it != fSignalParameters.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      results.push_back(it->second.at(i));
    }
  }
  for (auto it = fFluxParameters.begin();
       it != fFluxParameters.end(); ++it) {
    results.push_back(it->second);
  }

  for (auto it = fSystParameters.begin();
       it != fSystParameters.end(); ++it) {
    results.push_back(it->second.GetValue());
  }
  for (auto it = fG4RWParameters.begin();
       it != fG4RWParameters.end(); ++it) {
    results.push_back(it->second.GetValue());
  }

  return results;
}

void protoana::PDSPThinSliceFitter::ParameterScans() {
  fOutputFile.cd();
  TDirectory * out = (TDirectory *)fOutputFile.mkdir("Scans");
  out->cd();

  size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters + fTotalG4RWParameters;
  std::cout << "Total parameters " << total_parameters << std::endl;

  size_t start = 0;
  if (fOnlySystScans) {
    start = fTotalSignalParameters + fTotalFluxParameters;
  }
  else if (fOnlyG4RWScans) {
    start = fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters;
  }

  double * x = new double[fNScanSteps] {};
  double * y = new double[fNScanSteps] {};
  for (size_t i = start; i < total_parameters; ++i) {
    std::cout << "\tParameter " << fMinimizer->VariableName(i) << std::endl;
    bool scanned = fMinimizer->Scan(i, fNScanSteps, x, y);
    if (scanned) {
      TGraph gr(fNScanSteps - 1, x, y);
      gr.Write(fMinimizer->VariableName(i).c_str());
    }
  }

  delete[] x;
  delete[] y;
}


void protoana::PDSPThinSliceFitter::DoThrows(const TH1D & pars, const TMatrixD * cov) {
  std::vector<std::vector<double>> vals;
  TVectorD rand(pars.GetNbinsX());

  TMatrixD * cov_copy = 0x0;
  //TODO -- comment this
  if (fThrowType == "SingleBinThrow") {
    std::cout << "Uncorrelating" << std::endl;
    cov_copy = new TMatrixD(cov->GetNrows(), cov->GetNcols());
    for (int i = 0; i < cov_copy->GetNrows(); ++i) {
      if (i == fSingleThrowBin) {
        (*cov_copy)[i][i] = (*cov)[i][i];
      }
      else {
        (*cov_copy)[i][i] = 0.;
      }

      std::cout << i << " " << (*cov_copy)[i][i] << std::endl;
    }
  }
  else if (fThrowType == "CorrelateRangeThrow") {
    std::cout << "Uncorrelating" << std::endl;
    cov_copy = new TMatrixD(cov->GetNrows(), cov->GetNcols());
    for (int i = 0; i < cov_copy->GetNrows(); ++i) {
      for (int j = 0; j < cov_copy->GetNrows(); ++j) {
        if (i == j) {
          (*cov_copy)[i][j] = (*cov)[i][j];
        }
        else if ((i >= fRemainCorrRange.first) &&
            (i <= fRemainCorrRange.second) &&
            (j >= fRemainCorrRange.first) &&
            (j <= fRemainCorrRange.second)) {
          (*cov_copy)[i][j] = (*cov)[i][j];
        }
        else {
          (*cov_copy)[i][j] = 0.;
        }

        std::cout << i << " " << j << " " << (*cov_copy)[i][j] << std::endl;
      }
    }
  }

  TDecompChol chol((fThrowType == "CorrelateRangeThrow" ? *cov_copy : *cov));
  bool success = chol.Decompose();
  if (!success) {
    std::cout << "Error" << std::endl;
    return;
  }

  TMatrixD * cov_lower = (TMatrixD*)chol.GetU().Clone();
  cov_lower->Transpose(chol.GetU());

  (TDirectory*)fOutputFile.mkdir("Throws");
  fOutputFile.cd("Throws");

  std::vector<TVectorD *> throws_arrays;
  std::vector<TVectorD *> throws_int_arrays;
  std::vector<TVectorD *> throws_xsec_arrays;

  for (auto it = fSignalParameters.begin(); it != fSignalParameters.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      throws_arrays.push_back(new TVectorD(fNThrows));
      throws_xsec_arrays.push_back(new TVectorD(fNThrows));
      throws_int_arrays.push_back(new TVectorD(fNThrows));
    }
  }

  for (auto it = fFluxParameters.begin(); it != fFluxParameters.end(); ++it) {
    throws_arrays.push_back(new TVectorD(fNThrows));
  }

  for (auto it = fSystParameters.begin(); it != fSystParameters.end(); ++it) {
    throws_arrays.push_back(new TVectorD(fNThrows));
  }

  for (auto it = fG4RWParameters.begin(); it != fG4RWParameters.end(); ++it) {
    throws_arrays.push_back(new TVectorD(fNThrows));
  }


  TCanvas cPars("cParametersThrown", "");
  cPars.SetTicks();

  //Build up the map of selection hists
  std::map<int, std::vector<TH1*>> throw_hists;
  for (auto it = fDataSet.GetSelectionHists().begin();
       it != fDataSet.GetSelectionHists().end(); ++it) {
    throw_hists[it->first] = std::vector<TH1*>();
  }

  std::map<int, std::vector<TH1*>> truth_throw_hists;
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    truth_throw_hists[it->first] = std::vector<TH1*>();
  }

  std::map<int, std::vector<TH1*>> truth_inc_throw_hists,
                                   truth_xsec_throw_hists;
  for (auto it = fMeasurementSamples.begin();
       it != fMeasurementSamples.end(); ++it) {
    truth_inc_throw_hists[*it] = std::vector<TH1*>();
    truth_xsec_throw_hists[*it] = std::vector<TH1*>();
  }

  // SetupExtraHistsThrows();

  std::vector<double> throw_means(pars.GetNbinsX());
  //int n_signal_flux_pars = fTotalSignalParameters + fTotalFluxParameters;
  fCalcChi2InFCN = false;
  for (size_t i = 0; i < fNThrows; ++i) {
    if (! (i%100) ) {
      std::cout << "Throw " << i << "/" << fNThrows << std::endl;
    //auto start = std::chrono::high_resolution_clock::now();
    }
    vals.push_back(std::vector<double>(pars.GetNbinsX(), 1.));
    
    if (fThrowType == "SingleBinThrow") {
      GenerateUncorrelatedThrow(pars, cov_copy, vals.back());
    }
    else {
      GenerateCorrelatedThrow(pars, cov_lower, vals.back());
    }

    for (size_t j = 0; j < vals.back().size(); ++j) {
      (*throws_arrays[j])[i] = vals.back()[j];
      throw_means[j] += vals.back()[j];
    }

    //Applies the variations according to the thrown parameters
    fFillIncidentInFunction = true;
    fFitFunction(&(vals.back()[0]));
    fThinSliceDriver->GetCurrentHists(fSamples, fDataSet, throw_hists, fPlotRebinned);
    GetCurrentTruthHists(truth_throw_hists,
                         truth_inc_throw_hists,
                         truth_xsec_throw_hists);
    int xsec_bin = 0;
    for(auto it = truth_xsec_throw_hists.begin();
        it != truth_xsec_throw_hists.end(); ++it) {
      TH1D * xsec_hist = static_cast<TH1D*>(it->second.back());
      CalculateCrossSection(xsec_hist);

      TH1D * int_hist = (TH1D*)(truth_throw_hists.at(it->first).back());
      for (int j = 1; j <= xsec_hist->GetNbinsX(); ++j) {
        (*throws_xsec_arrays[xsec_bin])[i] = xsec_hist->GetBinContent(j);
        (*throws_int_arrays[xsec_bin])[i] = int_hist->GetBinContent(j);
        ++xsec_bin;
      }
    }


  }
  fCalcChi2InFCN = true;

  //throws_tree.Write();
  for (size_t i = 0; i < throws_arrays.size(); ++i) {
    std::string name = "array" + std::to_string(i);
    throws_arrays[i]->Write(name.c_str());
  }
  for (size_t i = 0; i < throws_xsec_arrays.size(); ++i) {
    std::string name = "xsec_array" + std::to_string(i);
    throws_xsec_arrays[i]->Write(name.c_str());
  }
  for (size_t i = 0; i < throws_int_arrays.size(); ++i) {
    std::string name = "int_array" + std::to_string(i);
    throws_int_arrays[i]->Write(name.c_str());
  }

  for (size_t i = 0; i < throw_means.size(); ++i) {
    throw_means[i] /= fNThrows;
  }

  fOutputFile.cd("Throws");

  //std::cout << "Calculating" << std::endl;
  TMatrixD recalculated_cov(pars.GetNbinsX(), pars.GetNbinsX());
  for (size_t i = 0; i < fNThrows; ++i) {
    for (int j = 0; j < pars.GetNbinsX(); ++j) {
      for (int k = 0; k < pars.GetNbinsX(); ++k) {
        //std::cout << i << " " << j << " " << k << std::endl;
        recalculated_cov[j][k]
            += (throw_means[j] - vals[i][j])*(throw_means[k] - vals[i][k]) /
               fNThrows;
      }
    }
  }

  recalculated_cov.Write("recalculated_cov");
  //std::cout << "Wrote" << std::endl;

  PlotThrows(throw_hists, truth_throw_hists,
             truth_inc_throw_hists, truth_xsec_throw_hists);

  fOutputFile.cd();
}

void protoana::PDSPThinSliceFitter::Do1DShifts(const TH1D & pars, bool prefit) {
  //std::string dir_name = (prefit ? "PreFit1DShifts" : "1DShifts");
  TDirectory * shift_dir = fOutputFile.mkdir((prefit ? "PreFit1DShifts" : "1DShifts"));
  fFillIncidentInFunction = true;

  std::cout << "1D shifts" << std::endl;
  //Iterate over the bins in the parameter hist
  for (int i = 1; i <= pars.GetNbinsX(); ++i) {
    
    std::string par_name = pars.GetXaxis()->GetBinLabel(i);
    std::cout << i << " " << par_name << std::endl;

    std::vector<double> vals;
    if (pars.GetBinError(i) < 1.e-9) continue;
    //Set to +1 sigma for parameter i
    for (int j = 1; j <= pars.GetNbinsX(); ++j) {
      if (j == i) {
        vals.push_back(pars.GetBinContent(j) + pars.GetBinError(i));
      }
      else {
        vals.push_back(pars.GetBinContent(j));
      }
    }
    fFitFunction(&vals[0]);

    shift_dir->cd();
    std::string dir_name = par_name + "_PlusSigma";
    TDirectory * plus_dir = shift_dir->mkdir(dir_name.c_str()); 
    TDirectory * xsec_dir = plus_dir->mkdir("XSec");
    CompareDataMC(dir_name, xsec_dir, plus_dir);

    vals.clear();

    //Set to -1 sigma for parameter i
    for (int j = 1; j <= pars.GetNbinsX(); ++j) {
      if (j == i) {
        vals.push_back(pars.GetBinContent(j) - pars.GetBinError(i));
      }
      else {
        vals.push_back(pars.GetBinContent(j));
      }
    }
    fFitFunction(&vals[0]);

    shift_dir->cd();
    dir_name = par_name + "_MinusSigma";
    TDirectory * minus_dir = shift_dir->mkdir(dir_name.c_str()); 
    xsec_dir = minus_dir->mkdir("XSec");
    CompareDataMC(dir_name, xsec_dir, minus_dir);
  }
}

void protoana::PDSPThinSliceFitter::NewRefillLoop(
  const std::vector<ThinSliceEvent> & events,
  std::vector<double> & mc_beam_fluxes,
  size_t worker_id, std::vector<size_t> n_events
) {
  
  size_t start_event = 0;
  for (size_t i = 0; i < worker_id; ++i) {
    start_event += n_events[i];
  }


  size_t end_event = start_event + n_events[worker_id];
  for (size_t i = start_event; i < end_event; ++i) {
    const ThinSliceEvent & event = events.at(i);
    int sample_ID = event.GetSampleID();
    double weight = 1.;
    if (fCoutLevel > 1) {
      std::cout << "Event " << i << event.GetSelectionID() << std::endl;
    }

    bool is_signal = (std::find(
      fMCDists.GetSignalIDs().begin(), fMCDists.GetSignalIDs().end(),
      sample_ID
    ) != fMCDists.GetSignalIDs().end());

    if (is_signal) {
      int signal_bin = fThinSliceStrategy->GetSignalBin(event, fMCDists);
      weight *= fSignalParameters.at(sample_ID).at(signal_bin);
    }


    int beam_bin = fThinSliceStrategy->GetBeamBin(fBeamEnergyBins, event, true);

    fThinSliceStrategy->FillHistsFromEvent(
      event, fMCDists, beam_bin, std::ref(fRefillMutex),
      weight);

    std::lock_guard<std::mutex> lock(fFluxMutex);
    mc_beam_fluxes[beam_bin] += weight;
  }


}


// TODO -- Call this something new
void protoana::PDSPThinSliceFitter::NewRefill(
  const std::vector<ThinSliceEvent> & events
) {

  fMCDists.Reset();
  // fFakeDataDists.Reset();

  size_t fNWorkers = 4;
  std::vector<double> mc_beam_fluxes(fBeamEnergyBins.size()-1, 0.);

  //Calculate events per workers
  size_t nevents = events.size() / fNWorkers;
  size_t remainder = events.size() % fNWorkers;
  std::vector<size_t> events_split(fNWorkers, nevents);
  while (remainder > fNWorkers) {
    nevents = nevents/fNWorkers;
    remainder = nevents % fNWorkers;

    for (size_t i = 0; i < fNWorkers; ++i) {
      events_split[i] += nevents;
    }
  }
  for (size_t i = 0; i < remainder; ++i) {
    events_split[i] += 1;
  }


  std::vector<std::thread> workers;
  for (auto workerid : ROOT::TSeqI(fNWorkers)) {
    std::thread worker(&PDSPThinSliceFitter::NewRefillLoop, this,
        std::ref(events), std::ref(mc_beam_fluxes), workerid, events_split);
    workers.emplace_back(std::move(worker));
  }

  for (auto &&worker : workers) { worker.join();}
  // if (fFakeResolution > 0. && !fUseStoredEnergies) fUseStoredEnergies = true; //turn this on so that every time after the first, the same energy is used
  // Go through and rescale according to the data flux
  for (size_t i = 0; i < mc_beam_fluxes.size(); ++i) {
    double ratio = ((fDataBeamFluxes[i] > 0.) ? fDataBeamFluxes[i] / mc_beam_fluxes[i] : 1.);
    if (fDataBeamFluxes[i] > 0.) {
      ratio = fDataBeamFluxes[i] / mc_beam_fluxes[i];
      std::cout << "Ratio " << i << " " << ratio << std::endl;
    }
    fMCDists.ScaleInBeamBin(i, ratio);
  }

  fThinSliceStrategy->CalcXSecs(fMCDists, 1.E27/ (fPitch * 1.4 * 6.022E23 / 39.948 ));
  fMCDists.CalcTotalDists();
}

void protoana::PDSPThinSliceFitter::DefineFitFunction() {
  ///use samples
  fFitFunction = ROOT::Math::Functor(
      [&](double const * coeffs) {
        auto new_time = std::chrono::high_resolution_clock::now();
        auto delta =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                new_time - fTime).count();
        fTime = new_time;
        if (fCoutLevel > 0)
          std::cout << "Fit step " << fNFitSteps << " took " << delta << std::endl;

        //Set all the parameters
        //coeffs are ordered according to
        //keys of signal par map
        //----vector for given signal
        //then keys of flux par map
        //then values of syst parameters
        size_t par_position = 0;
        //if (fDebugPars) {
        //  std::cout << "Pars: " << std::endl;
        //}
        for (auto it = fSignalParameters.begin();
             it != fSignalParameters.end(); ++it) {
          for (size_t i = 0; i < it->second.size(); ++i) {
            it->second.at(i) = coeffs[par_position];
            ++par_position;
            //if (fDebugPars) 
            //  std::cout << "\t" << par_position << " " << coeffs[par_position] << std::endl;
          }
        }
        for (auto it = fFluxParameters.begin();
             it != fFluxParameters.end(); ++it) {
          it->second = coeffs[par_position];
          ++par_position;
          //if (fDebugPars) 
          //  std::cout << "\t" << par_position << " " << coeffs[par_position] << std::endl;
        }

        for (auto it = fSystParameters.begin();
             it != fSystParameters.end(); ++it) {
          it->second.SetValue(coeffs[par_position]);
          //std::cout << it->first << " " << it->second.GetValue() << std::endl;
          //if (fDebugPars) 
          //  std::cout << "\t" << par_position << " " << coeffs[par_position] << std::endl;
          ++par_position;
        }
        for (auto it = fG4RWParameters.begin();
             it != fG4RWParameters.end(); ++it) {
          it->second.SetValue(coeffs[par_position]);
          //std::cout << it->first << " " << it->second.GetValue() << std::endl;
          //if (fDebugPars) 
          //  std::cout << "\t" << par_position << " " << coeffs[par_position] << std::endl;
          ++par_position;
        }

        //Refill the hists 
        // if (!fUseFakeSamples) {
        //   auto begin_time = std::chrono::high_resolution_clock::now();
        //   fThinSliceDriver->RefillMCSamples(
        //       fEvents, fSamples, fIsSignalSample,
        //       fBeamEnergyBins, fSignalParameters,
        //       fFluxParameters, fSystParameters,
        //       fG4RWParameters, 
        //       fFitUnderOverflow, fTieUnderOver, fScaleToDataBeamProfile,
        //       fFillIncidentInFunction,
        //       (fFixSamplesInFunction ? &fFixFactorHists : 0x0));
        //   auto end_time = std::chrono::high_resolution_clock::now();
        //   auto delta =
        //       std::chrono::duration_cast<std::chrono::microseconds>(
        //           end_time - begin_time).count();
        //   if (fCoutLevel > 0)
        //     std::cout << "Refilling took " << delta << " us" << std::endl;
        // }
        // else {
        //   fThinSliceDriver->RefillMCSamples(
        //       fFakeDataEvents, fFakeSamples, fIsSignalSample,
        //       fBeamEnergyBins, fSignalParameters,
        //       fFluxParameters, fSystParameters,
        //       fG4RWParameters, 
        //       fFitUnderOverflow, fTieUnderOver, fScaleToDataBeamProfile,
        //       fFillIncidentInFunction);
        // }
        NewRefill(fEvents);

        //Simpler renormalization...
        // if (fDataBeamFluxes.size()) {
        //   std::vector<double> fluxes(fDataBeamFluxes.size());
        //   for (auto it = (!fUseFakeSamples ? fSamples.begin() : fFakeSamples.begin());
        //        it != (!fUseFakeSamples ? fSamples.end() : fFakeSamples.end()); ++it) {
        //     for (size_t i = 0; i < it->second.size(); ++i) {
        //       for (size_t j = 0; j < it->second[i].size(); ++j) {
        //         fluxes[i] += it->second[i][j].GetVariedFlux();
        //       }
        //     }
        //   }


        //   std::vector<double> new_fluxes(fluxes.size());
        //   for (auto it = (!fUseFakeSamples ? fSamples.begin() : fFakeSamples.begin());
        //        it != (!fUseFakeSamples ? fSamples.end() : fFakeSamples.end()); ++it) {
        //     for (size_t i = 0; i < it->second.size(); ++i) {
        //       for (size_t j = 0; j < it->second[i].size(); ++j) {
        //         it->second[i][j].SetFactorAndScale(fDataBeamFluxes[i]/fluxes[i]);
        //         new_fluxes[i] += it->second[i][j].GetVariedFlux();
        //       }
        //     }
        //   }
        // }

        double chi2 = (
          fCalcChi2InFCN ?
          fThinSliceStrategy->CalcChi2(fMCDists, fDataSet) :
          0.
        );
        ++fNFitSteps;

        auto begin_time = std::chrono::high_resolution_clock::now();
        double syst_chi2 = (fAddSystTerm ? CalcChi2SystTerm() : 0.);
        auto end_time = std::chrono::high_resolution_clock::now();
        delta =
            std::chrono::duration_cast<std::chrono::microseconds>(
                end_time - begin_time).count();
        if (fCoutLevel > 0)
          std::cout << "Syst chi2 took " << delta << " us" << std::endl;

        double reg_term  = (fAddRegTerm ? CalcRegTerm() : 0.);
        double sig_term = (fAddSignalConstraint ? CalcSignalConstraint() : 0.);

        if (fDebugChi2)
          std::cout << chi2 << " " << syst_chi2 << std::endl;
        //if (chi2 < 0.) {
        if (chi2 < -1.*FLT_EPSILON ||
            std::isnan(chi2)) {
          std::string message = "Chi2 went negative " +
                                std::to_string(chi2) + "\n";


          size_t a = 0;
          for (auto it = fSignalParameters.begin();
               it != fSignalParameters.end(); ++it) {
            for (size_t i = 0; i < it->second.size(); ++i) {
               message += std::to_string(a) + " " +
                          std::to_string(it->second.at(i)) + "\n";
              ++a;
            }
          }
          for (auto it = fFluxParameters.begin();
               it != fFluxParameters.end(); ++it) {
            if (fCoutLevel > 0)
              std::cout << a << " " << it->second << std::endl;
            message += std::to_string(a) + " " + std::to_string(it->second) +
                       "\n";
            ++a;
          }

          for (auto it = fSystParameters.begin();
               it != fSystParameters.end(); ++it) {
            //it->second.SetValue(coeffs[a]);
            message += std::to_string(a) + " " +
                       std::to_string(it->second.GetValue()) + " " +
                       it->second.GetName() + "\n";
            ++a;
          }

          for (auto it = fG4RWParameters.begin();
               it != fG4RWParameters.end(); ++it) {
            //it->second.SetValue(coeffs[a]);
            message += std::to_string(a) + " " +
                       std::to_string(it->second.GetValue()) + " " +
                       it->second.GetName() + "\n";
            ++a;
          }
          throw std::runtime_error(message);
        }


        if (fSaveFitTree) {
          fOutputChi2Stat = chi2;
          fOutputChi2Syst = syst_chi2;
          fOutputParVals.clear();
          for (size_t ipar = 0; ipar < fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters + fTotalG4RWParameters; ++ipar) {
            fOutputParVals.push_back(coeffs[ipar]);
          }
          fOutputTree->Fill();
        }


        if (fCoutLevel > 0)
          std::cout << (chi2 + syst_chi2 + reg_term + sig_term) << std::endl;
        return (chi2 + syst_chi2 + reg_term + sig_term);
      },
      fTotalSignalParameters + fTotalFluxParameters + fTotalSystParameters + fTotalG4RWParameters);
      std::cout << "Done F2" << std::endl;

}

void protoana::PDSPThinSliceFitter::Configure(std::string fcl_file) {
  fhicl::ParameterSet pset;

  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;

  if (!fhicl_env) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing " <<
                 "or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};
  pset = fhicl::ParameterSet::make(fcl_file, lookupPolicy);

  //Getting the configurable parameters 
  if (fMCFileName == "")
    fMCFileName = pset.get<std::string>("MCFileName");
  if (fDataFileName == "")
    fDataFileName = pset.get<std::string>("DataFileName");

  fTreeName = pset.get<std::string>("TreeName");

  fSelectionSets = pset.get<std::vector<fhicl::ParameterSet>>("Selections");
  std::cout << "Selection Sets" << std::endl;

  fTrueIncidentBins = pset.get<std::vector<double>>("TrueIncidentBins");
  fBeamEnergyBins = pset.get<std::vector<double>>("BeamEnergyBins");
  fIncidentSamples = pset.get<std::vector<int>>("IncidentSamples");
  fMeasurementSamples = pset.get<std::vector<int>>("MeasurementSamples");

  fDrawXSecUnderflow = pset.get<bool>("DrawXSecUnderflow");
  
  fSampleSets = pset.get<std::vector<fhicl::ParameterSet>>("Samples");
  std::vector<std::pair<int, std::string>> temp_vec
      = pset.get<std::vector<std::pair<int, std::string>>>("FluxTypes");
  fFluxTypes = std::map<int, std::string>(temp_vec.begin(), temp_vec.end());

  fFitFlux = pset.get<bool>("FitFlux");
  auto temp_vec_fix = pset.get<std::vector<std::pair<int, double>>>(
      "FixFluxes", {});
  fFixFluxes = std::map<int, double>(temp_vec_fix.begin(), temp_vec_fix.end());

  if (fFitFlux) {
    for (size_t i = 0; i < temp_vec.size() - 1; ++i) {
      auto flux_id = temp_vec[i].first;
      fFluxParameterNames[flux_id] = "par_" + temp_vec[i].second +
                                               "_flux";
      fFluxParameters[flux_id] = 1.;

      if (fFixFluxes.find(flux_id) != fFixFluxes.end())  {
        fFluxParameters[flux_id] = fFixFluxes[flux_id];
      }
    }
    fTotalFluxParameters = temp_vec.size() - 1;
  }

  fMaxCalls = pset.get<int>("MaxCalls", 1e9);
  fMaxIterations = pset.get<int>("MaxIterations", 1e6);
  fPrintLevel = pset.get<int>("PrintLevel", 0);
  fCoutLevel = pset.get<int>("CoutLevel", 0);
  fNScanSteps = pset.get<unsigned int>("NScanSteps") + 1;
  fTolerance = pset.get<double>("Tolerance");
  fLowerLimit = pset.get<double>("LowerLimit");
  fUpperLimit = pset.get<double>("UpperLimit");
  fPlotStyle = pset.get<std::vector<std::pair<int,int>>>("PlotStyle");
  fPlotRebinned = pset.get<bool>("PlotRebinned");
  fRandomStart = pset.get<bool>("RandomStart", false);

  fDriverName = pset.get<std::string>("DriverName");
  fStrategyName = pset.get<std::string>("StrategyName");
  fAnalysisOptions = pset.get<fhicl::ParameterSet>("AnalysisOptions");
  fAltFluxVar = fAnalysisOptions.get<bool>("AltFluxVar", false);
  std::cout << "ALT FLUX VAR " << fAltFluxVar << std::endl;
  fPitch = fAnalysisOptions.get<double>("WirePitch");
  fPitchCorrection = fAnalysisOptions.get<double>("PitchCorrection", 0.9628);
  fPitch /= fPitchCorrection;
  fSliceMethod = fAnalysisOptions.get<std::string>("SliceMethod");
  fMultinomial = fAnalysisOptions.get<bool>("Multinomial", true);
  fDoFakeData = pset.get<bool>("DoFakeData");
  if (fDoFakeData) {
    fFakeDataRoutine = fAnalysisOptions.get<std::string>("FakeDataRoutine");
  }

  if (fNJobs > 0) {
    fAnalysisOptions.put_or_replace_compatible<size_t>("NWorkers", fNJobs);
  }
  fTrajZStart = fAnalysisOptions.get<double>("TrajZStart", 0.);
  fNoFakeReset = pset.get<bool>("NoFakeReset", false);
  fMaxEntries = pset.get<int>("MaxEntries", -1);
  fMaxDataEntries = pset.get<int>("MaxDataEntries", -1);
  fSplitMC = pset.get<bool>("SplitMC");
  fShuffle = pset.get<bool>("Shuffle", false);
  fDoThrows = pset.get<bool>("DoThrows");
  fDoThrowSystCentrals = pset.get<bool>("DoThrowSystCentrals", false);
  fSystsToThrow = pset.get<std::vector<size_t>>("SystsToThrow", {});

  fThrowType = pset.get<std::string>("ThrowType", "NormalThrow");
  fRemainCorrRange = pset.get<std::pair<int, int>>("RemainCorrRange", {-1, -1});
  fSingleThrowBin = pset.get<int>("SingleThrowBin", 0);

  fDoScans = pset.get<bool>("DoScans");
  fOnlySystScans = pset.get<bool>("OnlySystScans", false);
  fOnlyG4RWScans = pset.get<bool>("OnlyG4RWScans", false);
  fRunHesse = pset.get<bool>("RunHesse");
  fRerunFit = pset.get<bool>("RerunFit", false);
  fRequireGoodHesse = pset.get<bool>("RequireGoodHesse", true);
  fFitAttempts = pset.get<size_t>("FitAttempts", 1);
  fRunMinos1D = pset.get<bool>("RunMinos1D", false);
  fRunMinosConts = pset.get<bool>("RunMinosConts", false);
  fRunMultiConts = pset.get<bool>("RunMultiConts", false);
  fSignalContoursOnly = pset.get<bool>("SignalContoursOnly", true);
  fNContourPoints = pset.get<unsigned int>("NContourPoints", 10);
  fDo1DShifts = pset.get<bool>("Do1DShifts");
  fDoSysts = pset.get<bool>("DoSysts");
  fSetSystLimits = pset.get<bool>("SetSystLimits", true);
  // fSetSelVarLimits = pset.get<bool>("SetSelVarLimits", true);
  fSetSigLimits = pset.get<bool>("SetSigLimits", true);
  fRemoveSigThrowLimits = pset.get<bool>("RemoveSigThrowLimits", true);
  fFixVariables = pset.get<bool>("FixVariables", false);
  fFitUnderOverflow = pset.get<bool>("FitUnderOverflow", false);
  fGetMeanXSec = pset.get<bool>("GetMeanXSec", false);
  fTieUnderOver = pset.get<bool>("TieUnderOver", false);
  if (fFixVariables) {
    std::vector<std::pair<std::string, double>> temp_vec
        = pset.get<std::vector<std::pair<std::string, double>>>("SystsToFix");
    fSystsToFix = std::map<std::string, double>(temp_vec.begin(), temp_vec.end());
    temp_vec
        = pset.get<std::vector<std::pair<std::string, double>>>(
            "FixSystsPostFit", std::vector<std::pair<std::string, double>>());
    fFixSystsPostFit = std::map<std::string, double>(temp_vec.begin(), temp_vec.end());
  }

  fSetValsPreFit = pset.get<bool>("SetValsPreFit", false);
  if (fSetValsPreFit)
    fPreFitVals = pset.get<std::vector<double>>("PreFitVals");

  fDoFluctuateStats = pset.get<bool>("FluctuateStats");
  fFluctuateInSamples = pset.get<bool>("FluctuateInSamples", false);
  fVaryMCStats = pset.get<bool>("VaryMCStats", false);
  fVaryMCStatsForFakeData = pset.get<bool>("VaryMCStatsForFakeData", false);
  fUseMCStatVarWeight = pset.get<bool>("UseMCStatVarWeight", false);
  fUseMCStatVarWeightFakeData = pset.get<bool>("UseMCStatVarWeightFakeData", false);

  if (fVaryMCStats && fUseMCStatVarWeight) {
    std::string message = "PDSPThinSliceFitter::Configure: ";
    message += "Error. Can not use both 'VaryMCStats' and 'UseMCStatVarWeight' ";
    message += "at the same time.";
    throw std::runtime_error(message);
  }

  fDebugMCDataScale = pset.get<bool>("DebugMCDataScale", false);
  fDebugChi2 = pset.get<bool>("DebugChi2", false);
  fScaleToDataBeamProfile = pset.get<bool>("ScaleToDataBeamProfile", false);
  fFixPostFit = pset.get<bool>("FixPostFit", false);

  fNThrows = pset.get<size_t>("NThrows");
  fMaxRethrows = pset.get<size_t>("MaxRethrows");

  fFitType = pset.get<std::string>("FitType");

  // fXSecCalcStyle = pset.get<std::string>("XSecCalcStyle", "Approx");

  if (fFitType == "ThrowsOnly") {
    fFitResultsFile = pset.get<std::string>("PriorFitResults");
  }


  //fDoScaleDataToNorm = pset.get<bool>("ScaleDataToNorm", false);

  fRNG = TRandom3(pset.get<int>("RNGSeed", 0));

  //fMCXSecFileName = pset.get<std::string>("MCXSecFile", "");

  fAddRegTerm = pset.get<bool>("AddRegTerm", false);
  fRegFactor = pset.get<double>("RegFactor", 0.);

  // fExtraHistSets = fAnalysisOptions.get<std::vector<fhicl::ParameterSet>>("ExtraHists", {});

  fAddDiffInQuadrature = pset.get<bool>("AddDiffInQuadrature", false);
  if (fAddDiffInQuadrature) {
    //std::vector<std::pair<int, std::string>> temp_vec
    //    = pset.get<std::vector<std::pair<int, std::string>>>("Diffs");
    //fDiffGraphs = std::map<int, std::string>(temp_vec.begin(), temp_vec.end());
    fDiffCovName = pset.get<std::string>("DiffCovName");
    fDiffGraphFile = pset.get<std::string>("DiffGraphFile");
  }

  fAddSignalConstraint = pset.get<bool>("AddSignalConstraint", false);
  fSignalConstraintSize = pset.get<double>("SignalConstraintSize", 10.);
  
  fMCDists.SetIDs(pset);
  fFakeDataDists.SetIDs(pset);
  try {
    fThinSliceStrategy = protoana::ThinSliceStrategyRegistry::Instance()->GetStrategy(
      fStrategyName, fAnalysisOptions);
  }
  catch (const std::runtime_error & e) {
    protoana::ThinSliceStrategyRegistry::Instance()->PrintAvailableStrategies();
    throw e;
  }
  fThinSliceStrategy->BuildDists(fMCDists, pset);
  fThinSliceStrategy->BuildDists(fFakeDataDists, pset, "Fake");
  
  //Initialize systematics
  if (fDoSysts) {
    fAddSystTerm = pset.get<bool>("AddSystTerm");
    std::vector<fhicl::ParameterSet> par_vec
        = pset.get<std::vector<fhicl::ParameterSet>>("Systematics");
    for (size_t i = 0; i < par_vec.size(); ++i) {
      ThinSliceSystematic syst(par_vec[i]);
      fSystParameters[par_vec[i].get<std::string>("Name")] = syst;
      ++fTotalSystParameters;
    }
    fThrowAndFixSysts = pset.get<bool>("ThrowAndFixSysts", false);
    fSystsToThrowAndFix
        = pset.get<std::vector<std::string>>("SystsToThrowAndFix", {});
    if (fThrowAndFixSysts && fSystsToThrowAndFix.empty()) {
      std::string message = "PDSPThinSliceFitter::Configure: ";
      message += "Error. Requesting to Throw and fix systs, ";
      message += "but provided empty vector of systs to throw and fix";
      throw std::runtime_error(message);
    }


    //Moving this
    //if (fAddSystTerm) {
      std::vector<std::pair<std::string, int>> temp_vec
          = pset.get<std::vector<std::pair<std::string, int>>>("CovarianceBins");
      fCovarianceBins
          = std::map<std::string, size_t>(temp_vec.begin(), temp_vec.end());

      for (auto it = fSystParameters.begin(); it != fSystParameters.end();
           ++it) {
        fCovarianceBinsSimple.push_back(fCovarianceBins[it->first]);
      }

      TFile * cov_file = TFile::Open(pset.get<std::string>("CovarianceFile").c_str());
      std::cout << "Getting covariance from " <<
                   pset.get<std::string>("CovarianceFile").c_str() << std::endl;

      //Turn this into a temporary TMatrixD
      TMatrixD* base_matrix = (TMatrixD*)cov_file->Get(
          pset.get<std::string>("CovarianceMatrix").c_str());
      int n_base_rows = base_matrix->GetNrows();

      fCovMatrix = new TMatrixD(n_base_rows, n_base_rows);
      for (int i = 0; i < n_base_rows; ++i) {
        for (int j = 0; j < n_base_rows; ++j) {
          (*fCovMatrix)[i][j] = (*base_matrix)[i][j];
        }
      }



      fCovMatrixDisplay = (TMatrixD*)fCovMatrix->Clone();
      if (!fCovMatrix->IsSymmetric()) {
        std::string message = "PDSPThinSliceFitter::Configure: ";
        message += "Error. Input covariance matrix ";
        message += "is not symmetric";
        throw std::runtime_error(message);
      }

      fInputChol = new TDecompChol(*fCovMatrix);
      fInputChol->Decompose();

      size_t nrows = static_cast<size_t>(fCovMatrix->GetNrows());
      if (fTotalSystParameters != nrows) {
        std::string message = "PDSPThinSliceFitter::Configure: ";
        message += "Error. Input covariance matrix and nsyst differ\n";
        message += std::to_string(fTotalSystParameters);
        message += " " + std::to_string(fCovMatrix->GetNrows());
        throw std::runtime_error(message);
      }


      fCovMatrix->Invert();
      cov_file->Close();


      //Here, throw a new central according 
      //Adding syst term check here
      if (fAddSystTerm && fDoThrowSystCentrals) ThrowSystCentrals();

    //}
  }

  //New initialize g4rw systematics
  std::vector<fhicl::ParameterSet> par_vec
      = pset.get<std::vector<fhicl::ParameterSet>>("G4RWParameters", {});
  fTuneG4RWPars = pset.get<bool>("TuneG4RWPars", false);
  fNominalG4RWFakeData = pset.get<bool>("NominalG4RWFakeData", false);
  for (size_t i = 0; i < par_vec.size(); ++i) {
    ThinSliceSystematic syst(par_vec[i]);
    fG4RWParameters[par_vec[i].get<std::string>("Name")] = syst;
    ++fTotalG4RWParameters;
  }
}

void protoana::PDSPThinSliceFitter::ThrowSystCentrals() {
  std::cout << "Throwing Syst Centrals" << std::endl;
  TMatrixD * cov_lower = (TMatrixD*)fInputChol->GetU().Clone();
  cov_lower->Transpose(fInputChol->GetU());

  size_t n_cov_rows = cov_lower->GetNrows();
  if (n_cov_rows != fTotalSystParameters) {
    std::cout << "Error: " << n_cov_rows << " " << fTotalSystParameters <<
                 std::endl;
  }

  std::vector<double> init_vals, new_centrals;
  std::vector<double> throw_limits, throw_limits_up;
  for (auto & it : fSystParameters) {
    init_vals.push_back(it.second.GetCentral());
    throw_limits.push_back(it.second.GetGenThrowLimit());
    throw_limits_up.push_back(it.second.GetGenThrowLimitUp());
  }


  bool rethrow = true;
  while (rethrow) {
    new_centrals.clear();
    bool all_pos = true;
    TVectorD rand(fTotalSystParameters);
    for (size_t i = 0; i < fTotalSystParameters; ++i) {

      bool throw_this = (std::find(fSystsToThrow.begin(),
                                   fSystsToThrow.end(), i) !=
                         fSystsToThrow.end());
      rand[i] = (throw_this ? fRNG.Gaus() : 0.);
    }
    TVectorD rand_times_chol = (*cov_lower)*rand;

    for (size_t i = 0; i < fTotalSystParameters; ++i) {
      double val = rand_times_chol[fCovarianceBinsSimple[i]] +
                   init_vals[i];
      bool below = (val < throw_limits[i]);
      bool above = (val > throw_limits_up[i]);
      if (below || above) {
        all_pos = false;
        std::cout << "Rethrowing " << i << " " << val << " " <<
                     (below ? throw_limits[i] : throw_limits_up[i]) <<
                     std::endl;
        break;
      }
      new_centrals.push_back(val);
    }

    if (all_pos) {
      size_t i_par = 0;
      for (auto & it : fSystParameters) {
        it.second.SetCentral(new_centrals[i_par]);
        it.second.SetValue(new_centrals[i_par]);
        ++i_par;
      }
    }
    rethrow = !all_pos;
  }
}

void protoana::PDSPThinSliceFitter::SetBestFit() {
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        it->second[i][j].SetBestFit();
      }
    }
  }
}

void protoana::PDSPThinSliceFitter::PlotThrows(
    std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<TH1*>> & truth_throw_hists,
    std::map<int, std::vector<TH1*>> & truth_inc_hists,
    std::map<int, std::vector<TH1*>> & truth_xsec_hists) {
  std::map<int, TH1*> data_hists
      = (fPlotRebinned ?
         fDataSet.GetRebinnedSelectionHists() :
         fDataSet.GetSelectionHists());

  //Build best fit hists and get bins for covariance 
  //std::map<int, TH1D*> best_fit_selection_hists;
  int nBins = 0;
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it ) {
    //TH1D * best_fit_hist = (TH1D*)it->second->Clone();
    //best_fit_hist->Reset();
    //for (auto it2 = fSamples.begin(); it2 != fSamples.end(); ++it2) {
    //  for (size_t i = 0; i < it2->second.size(); ++i) {
    //    for (size_t j = 0; j < it2->second[i].size(); ++j) {
    //      it2->second[i][j].SetFactorToBestFit();
    //      best_fit_hist->Add(
    //          (TH1D*)(fPlotRebinned ?
    //                  it2->second[i][j].GetRebinnedSelectionHist(it->first) :
    //                  it2->second[i][j].GetSelectionHist(it->first)));
    //    }
    //  }
    //}
    //best_fit_selection_hists[it->first] = best_fit_hist;
    //nBins += best_fit_hist->GetNbinsX();
    nBins += it->second->GetNbinsX();
  }

  TH2D selection_cov("SelectionCov", "", nBins, 0, nBins, nBins, 0, nBins);

  nBins = 0;
  std::map<int, size_t> sample_bins;
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    nBins += it->second[0].size();
    sample_bins[it->first] = it->second[0].size();
  }

  //std::map<int, std::vector<double>> best_fit_truth;
  std::map<int, std::vector<double>> best_fit_errs;

  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    //best_fit_truth[it->first]
    //    = std::vector<double>(sample_bins[it->first], 0.);
    best_fit_errs[it->first]
        = std::vector<double>(sample_bins[it->first], 0.);
   
    //for (size_t i = 0; i < sample_bins[it->first]; ++i) {
    //  double best_fit_val_i = 0.;
    //  for (size_t j = 0; j < it->second.size(); ++j) {
    //    best_fit_val_i += it->second[j][i].GetVariedFlux();
    //  }

    //  best_fit_truth[it->first][i] = best_fit_val_i;
    //}
  }

  TH2D interaction_cov("interaction_cov", "", nBins, 0, nBins, nBins, 0, nBins);
  std::map<int, std::vector<double>> best_fit_inc_truth;
  std::map<int, std::vector<double>> best_fit_xsec_truth;
  std::map<int, std::vector<double>> best_fit_inc_errs;
  std::map<int, std::vector<double>> best_fit_xsec_errs;

  nBins = 0;
  std::map<int, size_t> xsec_bins;
  for (auto it = fBestFitIncs.begin(); it != fBestFitIncs.end(); ++it) {
    int s = it->first;
    nBins += it->second->GetNbinsX();
    xsec_bins[s] = it->second->GetNbinsX();

    best_fit_inc_truth[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_xsec_truth[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_inc_errs[s] = std::vector<double>(xsec_bins[s], 0.);
    best_fit_xsec_errs[s] = std::vector<double>(xsec_bins[s], 0.);
    
    for (size_t i = 0; i < xsec_bins[s]; ++i) {
      best_fit_inc_truth[s][i] = it->second->GetBinContent(i+1);
      best_fit_xsec_truth[s][i] = fBestFitXSecs[s]->GetBinContent(i+1);
    }
  }

  //TH2D incident_cov("incident_cov", "", nBins, 0, nBins, nBins, 0, nBins);
  TH2D xsec_cov("xsec_cov", "", nBins, 0, nBins, nBins, 0, nBins);
  TH2D xsec_corr("xsec_corr", "", nBins, 0, nBins, nBins, 0, nBins);
  xsec_cov.SetTitle(";Cross Section Bin;Cross section Bin");
  xsec_cov.SetTitleSize(.05, "XY");
  xsec_cov.GetXaxis()->SetTitleOffset(.8);
  xsec_corr.SetTitle(";Cross Section Bin;Cross section Bin");
  xsec_corr.SetTitleSize(.05, "XY");
  xsec_corr.GetXaxis()->SetTitleOffset(.8);
  xsec_corr.SetMaximum(1.);
  xsec_corr.SetMinimum(-1.);
  TMatrixD xsec_cov_matrix(nBins, nBins);

  std::map<int, std::vector<double>> mean_xsecs;
  for (size_t z = 0; z < fNThrows; ++z) {
    for (auto it = truth_xsec_hists.begin(); it != truth_xsec_hists.end(); ++it) {
      std::vector<TH1 *> xsec_hists_i = it->second;
      if (z == 0)
        mean_xsecs[it->first] = std::vector<double>(xsec_bins[it->first], 0.);
      for (size_t i = 0; i < xsec_bins[it->first]; ++i) {
        mean_xsecs[it->first][i] += xsec_hists_i[z]->GetBinContent(i+1)/fNThrows;
      }
    }
  }

  for (size_t z = 0; z < fNThrows; ++z) {
    int bin_i = 1;
    for (auto it = fBestFitSelectionHists.begin();
         it != fBestFitSelectionHists.end(); ++it) {
      TH1D * best_fit = it->second;
      int selection_ID = it->first;
      std::vector<TH1*> & temp_throws = throw_hists[selection_ID];
      for (int i = 1; i <= best_fit->GetNbinsX(); ++i) {
        double best_fit_val_i = best_fit->GetBinContent(i);
        int bin_j = 1;
        for (auto it2 = fBestFitSelectionHists.begin();
             it2 != fBestFitSelectionHists.end(); ++it2) {

          TH1D * best_fit_2 = it2->second;
          int selection_ID_2 = it2->first;
          std::vector<TH1*> & temp_throws_2 = throw_hists[selection_ID_2];
          for (int j = 1; j <= best_fit_2->GetNbinsX(); ++j) {
            double best_fit_val_j = best_fit_2->GetBinContent(j);
            double val = (best_fit_val_i - temp_throws[z]->GetBinContent(i))*
                         (best_fit_val_j - temp_throws_2[z]->GetBinContent(j));
            selection_cov.SetBinContent(
                bin_i, bin_j, (val/temp_throws.size() +
                               selection_cov.GetBinContent(bin_i, bin_j)));
            ++bin_j;
          }
        }
        ++bin_i;
      }
    }

    bin_i = 1;
    for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
      std::vector<TH1 *> throw_hists_i = truth_throw_hists[it->first];
     
      for (size_t i = 0; i < sample_bins[it->first]; ++i) {
        double best_fit_val_i = fBestFitTruthVals[it->first][i];

        int bin_j = 1;
        for (auto it2 = fSamples.begin(); it2 != fSamples.end(); ++it2) {
          std::vector<TH1 *> throw_hists_j = truth_throw_hists[it2->first];
          for (size_t j = 0; j < sample_bins[it2->first]; ++j) {
            double best_fit_val_j = fBestFitTruthVals[it2->first][j];

            double val
                = (throw_hists_i[z]->GetBinContent(i+1) - best_fit_val_i)*
                  (throw_hists_j[z]->GetBinContent(j+1) - best_fit_val_j);
            interaction_cov.SetBinContent(
                bin_i, bin_j,
                (interaction_cov.GetBinContent(bin_i, bin_j) +
                 val/throw_hists_i.size()));
            if (bin_i == bin_j && (z == fNThrows - 1)) {
              best_fit_errs[it->first][i]
                  = sqrt(interaction_cov.GetBinContent(bin_i, bin_j));
            }
            ++bin_j;
          }
        }

        ++bin_i;
      }
    }




    bin_i = 1;
    for (auto it = truth_inc_hists.begin(); it != truth_inc_hists.end(); ++it) {
      std::vector<TH1 *> xsec_hists_i = truth_xsec_hists[it->first];
      for (size_t i = 0; i < xsec_bins[it->first]; ++i) {
        //double best_fit_xsec_i = mean_xsecs[it->first][i];
        double best_fit_xsec_i = (fGetMeanXSec ?
                                  mean_xsecs[it->first][i] :
                                  best_fit_xsec_truth[it->first][i]);

        int bin_j = 1;
        for (auto it2 = truth_inc_hists.begin(); it2 != truth_inc_hists.end();
             ++it2) {
          std::vector<TH1 *> xsec_hists_j = truth_xsec_hists[it2->first];
          for (size_t j = 0; j < xsec_bins[it2->first]; ++j) {
            //double best_fit_xsec_j = mean_xsecs[it2->first][j];
            double best_fit_xsec_j = (fGetMeanXSec ?
                                      mean_xsecs[it2->first][j] :
                                      best_fit_xsec_truth[it2->first][j]);

            double val
                = (xsec_hists_i[z]->GetBinContent(i+1) - best_fit_xsec_i)*
                  (xsec_hists_j[z]->GetBinContent(j+1) - best_fit_xsec_j);
            xsec_cov.SetBinContent(
                bin_i, bin_j,
                (xsec_cov.GetBinContent(bin_i, bin_j) +
                 val/fNThrows));
            if (bin_i == bin_j && (z == fNThrows - 1)) {
              best_fit_xsec_errs[it->first][i]
                  = sqrt(xsec_cov.GetBinContent(bin_i, bin_j));
            }
            ++bin_j;
          }
        }
        ++bin_i;
      }
    }
  }


  //Add in quadrature here   
  if (fAddDiffInQuadrature) {
    std::cout << "Adding diff in quad" << std::endl;

    std::map<int, TGraph*> diff_graph_map;
    TFile diff_file(fDiffGraphFile.c_str(), "open");

    fDiffCov = (TH2D*)diff_file.Get(fDiffCovName.c_str());

    /*
    for (auto it = fDiffGraphs.begin(); it != fDiffGraphs.end(); ++it) {
      diff_graph_map[it->first] = (TGraph*)diff_file.Get(it->second.c_str());
      std::cout << it->first << " " << diff_graph_map[it->first] << std::endl;
    }*/
    
    xsec_cov.Add(fDiffCov);

    int diag_bin = 1;
    for (auto it = diff_graph_map.begin(); it != diff_graph_map.end(); ++it) {
      std::cout << it->first << " " << it->second->GetN() << " " <<
                   best_fit_xsec_errs[it->first].size() << std::endl;
      for (int i = 0; i < it->second->GetN(); ++i) {
        //double diff_err = std::pow(it->second->GetY()[i], 2);
        //double best_fit_err = std::pow(best_fit_xsec_errs[it->first][i], 2);
        //std::cout << "Adding " << diff_err << " " << best_fit_err << " " <<
        //             sqrt(best_fit_err + diff_err) << " " <<
        //             sqrt(best_fit_err) << std::endl;
        //xsec_cov.SetBinContent(diag_bin, diag_bin,
        //                       (xsec_cov.GetBinContent(diag_bin, diag_bin)
        //                        + diff_err));
        //best_fit_xsec_errs[it->first][i]
        //    = sqrt(best_fit_err + diff_err);
        best_fit_xsec_errs[it->first][i] = xsec_cov.GetBinContent(diag_bin,  diag_bin);

        diag_bin++;
      }
    }
  }

  for (int i = 0; i < nBins; ++i) {
    for (int j = 0; j < nBins; ++j) {
      xsec_cov_matrix[i][j] = xsec_cov.GetBinContent(i+1, j+1);
      double corr_val =
          xsec_cov.GetBinContent(i+1, j+1)/
              sqrt(xsec_cov.GetBinContent(i+1, i+1)*
                   xsec_cov.GetBinContent(j+1, j+1));
      xsec_corr.SetBinContent(i+1, j+1, corr_val);
    }
  }
  xsec_cov_matrix.Invert();


  fOutputFile.cd("Throws");
  selection_cov.Write();
  interaction_cov.Write();
  xsec_cov.Write();
  xsec_corr.Write();

  TCanvas cXSecCorr("cXSecCorr", "");
  xsec_corr.Draw("colz");
  cXSecCorr.Write();



  double nominal_xsec_chi2 = 0.;
  double fake_xsec_chi2 = 0.;
  int bin_i = 0;

  auto & measured_xsec = (fGetMeanXSec ? mean_xsecs : best_fit_xsec_truth);

  for (auto it = measured_xsec.begin(); it != measured_xsec.end();
       ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      double measured_val_i = it->second[i];
      double mc_val_i = fNominalXSecs[it->first]->GetBinContent(i+1);
      double fake_val_i = ((fDoFakeData /*&& fFakeDataRoutine != "Toy"*/) ?
                           fFakeDataXSecs[it->first]->GetBinContent(i+1) :
                           0.);
      int bin_j = 0;
      for (auto it2 = measured_xsec.begin(); it2 != measured_xsec.end(); ++it2) {
        for (size_t j = 0; j < it2->second.size(); ++j) {
          double measured_val_j = it2->second[j];
          double mc_val_j = fNominalXSecs[it2->first]->GetBinContent(j+1);
          double fake_val_j = ((fDoFakeData /*&& fFakeDataRoutine != "Toy"*/) ?
                               fFakeDataXSecs[it2->first]->GetBinContent(j+1) :
                               0.);
          nominal_xsec_chi2 += ((measured_val_i - mc_val_i)*
                                xsec_cov_matrix[bin_i][bin_j]*
                                (measured_val_j - mc_val_j));
          /*std::cout << bin_i << " " << bin_j << " " << mc_val_i << " " << mc_val_j <<
                       " " << nominal_xsec_chi2 << " " <<
                       xsec_cov_matrix[bin_i][bin_j] << std::endl;*/
          if (fDoFakeData /*&& fFakeDataRoutine != "Toy"*/) {
            fake_xsec_chi2 += ((measured_val_i - fake_val_i)*
                               xsec_cov_matrix[bin_i][bin_j]*
                               (measured_val_j - fake_val_j));
          }
          ++bin_j;
        }
        //++bin_j;
      }
      ++bin_i;
    }
    //++bin_i;
  }

  int bin_count = 0;
  for (auto it = data_hists.begin(); it != data_hists.end(); ++it) {
    int selection_ID = it->first;
    std::vector<TH1*> hists = throw_hists.at(selection_ID);

    std::string canvas_name = "cThrow" +
                              fDataSet.GetSelectionName(selection_ID);
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();

    std::string name = "Throw" + fDataSet.GetSelectionName(selection_ID);
    auto data_hist = it->second;
    std::vector<double> xs, xs_width;
    std::vector<double> ys, errs;
    for (int i = 1;
         i <= fBestFitSelectionHists[it->first]->GetNbinsX(); ++i) {
      ys.push_back(
          fBestFitSelectionHists[it->first]->GetBinContent(i));
      errs.push_back(
          sqrt(selection_cov.GetBinContent(bin_count+i, bin_count+i)));
      xs.push_back(data_hist->GetBinCenter(i));
      xs_width.push_back(data_hist->GetBinWidth(i)/2.);
    } 

    TGraphAsymmErrors throw_gr(data_hist->GetNbinsX(),
                               &xs[0], &ys[0], 
                               &xs_width[0], &xs_width[0], &errs[0], &errs[0]);

    throw_gr.SetTitle(fDataSet.GetSelectionName(selection_ID).c_str());
    throw_gr.GetXaxis()->SetTitle(data_hist->GetXaxis()->GetTitle());
    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    data_hist->Draw();
    throw_gr.Draw("same a2");
    data_hist->Draw("same e1");
    fOutputFile.cd("Throws");
    cThrow.Write();

    bin_count += data_hist->GetNbinsX();
  }

  bin_count = 0;
  for (auto it = truth_throw_hists.begin(); it != truth_throw_hists.end(); ++it) {
    int sample_ID = it->first;

    std::vector<double> xs, xs_width;
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      xs.push_back(i + 0.5);
      xs_width.push_back(.5);
    }

    std::string name = "hNominal" + fSamples[sample_ID][0][0].GetName();
    TH1D temp_nominal(name.c_str(), "", xs.size(), 0, xs.size());
    std::string dummy_name = "hDummy" + fSamples[sample_ID][0][0].GetName();
    TH1D dummy(dummy_name.c_str(), "", xs.size(), 0, xs.size());
    if (fIsSignalSample[sample_ID]) {
      dummy.GetXaxis()->SetBinLabel(1, "Underflow");
      dummy.GetXaxis()->SetBinLabel(dummy.GetNbinsX(), "Overflow");
      for (int i = 2; i < dummy.GetNbinsX(); ++i) {
        std::ostringstream low_stream, high_stream;
        low_stream.precision(2);
        high_stream.precision(2);
        low_stream << std::fixed <<
                      fSamples[sample_ID][0][i-1].RangeLowEnd();
        high_stream
            << std::fixed <<
               fSamples[sample_ID][0][i-1].RangeHighEnd();
        std::string label = low_stream.str();
        label += " - ";
        label += high_stream.str();
        //std::string label
        //    = std::to_string(fSamples[sample_ID][0][i-1].RangeLowEnd());
        //label += " - ";
        //label += std::to_string(fSamples[sample_ID][0][i-1].RangeHighEnd());
        dummy.GetXaxis()->SetBinLabel(i, label.c_str());
      }
    }
    else {
      dummy.GetXaxis()->SetBinLabel(1, "");
    }


    std::vector<std::vector<ThinSliceSample>> & samples_vec_2D
        = fSamples[sample_ID];
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {
      for (size_t j = 0; j < samples_vec_2D[i].size(); ++j) {
        temp_nominal.AddBinContent(j+1, samples_vec_2D[i][j].GetNominalFlux());
      }
    }

    double max = -999.;
    for (size_t i = 0; i < sample_bins[it->first]; ++i) {
      if ((fBestFitTruthVals[sample_ID][i] + best_fit_errs[sample_ID][i]) > max)
        max = (fBestFitTruthVals[sample_ID][i] + best_fit_errs[sample_ID][i]);

      if (temp_nominal.GetBinContent(i+1) > max)
        max = temp_nominal.GetBinContent(i+1);
    }

    fOutputFile.cd("Throws");
    std::string canvas_name = "cTruthThrow" + fSamples[sample_ID][0][0].GetName();
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();
    TGraphAsymmErrors throw_gr(xs.size(),
                                &xs[0], &fBestFitTruthVals[it->first][0], 
                                &xs_width[0], &xs_width[0],
                                &best_fit_errs[it->first][0],
                                &best_fit_errs[it->first][0]);
    throw_gr.SetFillStyle(3144);
    throw_gr.SetFillColor(kRed);
    throw_gr.SetMinimum(0.);
    throw_gr.SetMaximum(1.5*max);
    dummy.SetMaximum(1.5*max);
    dummy.Draw();
    //throw_gr.Draw("a2");
    throw_gr.Draw("same 2");
    throw_gr.Draw("p");

    temp_nominal.SetMarkerColor(kBlue);
    temp_nominal.SetMarkerStyle(20);
    temp_nominal.Draw("same p");

    TLegend leg;
    leg.AddEntry(&throw_gr, "Throws", "lpf");
    leg.AddEntry(&temp_nominal, "Nominal", "p");

    if (fDoFakeData && fFakeDataRoutine != "Toy" && fFakeDataRoutine != "G4RWGrid") {
      name = "hVaried" + fSamples[sample_ID][0][0].GetName();
      TH1D * temp_varied = (TH1D*)temp_nominal.Clone(name.c_str());
      for (size_t i = 0; i < xs.size(); ++i) {
        temp_varied->SetBinContent(
            i+1, temp_varied->GetBinContent(i+1)*fFakeDataScales[sample_ID][i]);
      }
      temp_varied->SetMarkerColor(kBlack);
      temp_varied->SetMarkerStyle(20);
      temp_varied->Draw("same p");
      leg.AddEntry(temp_varied, "Fake Data", "p");
    }

    leg.Draw();
    cThrow.Write();

    bin_count += xs.size();
  }

  //for (auto it = best_fit_xsec_truth.begin(); it != best_fit_xsec_truth.end();
  for (auto it = measured_xsec.begin(); it != measured_xsec.end(); ++it) {
    int sample_ID = it->first;

    std::vector<double> xs, xs_width;
    for (size_t i = 0; i < xsec_bins[sample_ID]; ++i) {
      xs.push_back(i + 0.5);
      xs_width.push_back(0.);
    }

    fOutputFile.cd("Throws");
    std::string canvas_name = "cXSecThrow" + fSamples[sample_ID][0][0].GetName();
    std::string gr_name = "grXSecThrow" + fSamples[sample_ID][0][0].GetName();
    TCanvas cThrow(canvas_name.c_str(), "");
    cThrow.SetTicks();

    std::string dummy_name = "hDummyXSec" + fSamples[sample_ID][0][0].GetName();
    TH1D dummy(dummy_name.c_str(), "", xs.size(), 0, xs.size());
    for (int i = 1; i <= dummy.GetNbinsX(); ++i) {
      std::ostringstream low_stream, high_stream;
      low_stream.precision(2);
      high_stream.precision(2);
      low_stream << std::fixed <<
                    fSamples[sample_ID][0][i].RangeLowEnd();
      high_stream
          << std::fixed <<
             fSamples[sample_ID][0][i].RangeHighEnd();
      std::string label = low_stream.str();
      label += " - ";
      label += high_stream.str();
      //std::string label
      //    = std::to_string(fSamples[sample_ID][0][i-1].RangeLowEnd());
      //label += " - ";
      //label += std::to_string(fSamples[sample_ID][0][i-1].RangeHighEnd());
      dummy.GetXaxis()->SetBinLabel(i, label.c_str());
    }

    TGraphAsymmErrors throw_gr(xs.size(),
                               &xs[0], &measured_xsec[it->first][0],//&best_fit_xsec_truth/*mean_xsecs*/[it->first][0], 
                               &xs_width[0], &xs_width[0],
                               &best_fit_xsec_errs[it->first][0],
                               &best_fit_xsec_errs[it->first][0]);
    TGraph mean_gr(xs.size(), &xs[0], &mean_xsecs[it->first][0]);
    double max = -999.;
    for (size_t i = 0; i < measured_xsec/*best_fit_xsec_truth*/[it->first].size(); ++i) {
      if ((measured_xsec/*best_fit_xsec_truth*/[it->first][i] +
           best_fit_xsec_errs[it->first][0]) > max) {
        max = (measured_xsec/*best_fit_xsec_truth*/[it->first][i] +
               best_fit_xsec_errs[it->first][0]);
      }
    }
    //throw_gr.SetFillStyle(3144);
    //throw_gr.SetFillColor(kRed);
    throw_gr.SetMarkerStyle(20);
    throw_gr.SetMinimum(0.);
    dummy.SetMaximum(1.5*max);
    dummy.SetMinimum(0.);
    dummy.Draw();
    dummy.GetXaxis()->SetTitle("Kinetic Energy (MeV)");
    dummy.GetYaxis()->SetTitle("#sigma (mb)");
    //throw_gr.Draw("same 2");
    throw_gr.Draw("same 2");
    //mean_gr.Draw("same 2");
    throw_gr.Draw("p");

    std::vector<double> nominal_xsec_vals;
    for (size_t i = 0; i < xs.size(); ++i) {
      nominal_xsec_vals.push_back(
          fNominalXSecs[it->first]->GetBinContent(i+1));
    }
    TGraph nominal_gr(xs.size(), &xs[0], &nominal_xsec_vals[0]);
    nominal_gr.SetMarkerColor(kBlue);
    nominal_gr.SetMarkerStyle(20);
    nominal_gr.Draw("same p");

    TLegend leg;
    leg.AddEntry(&throw_gr, "Measured", "lf");
    //leg.AddEntry(&mean_gr, "Mean", "lpf");
    leg.AddEntry(&nominal_gr, "Nominal", "p");

    if (fDoFakeData /*&& fFakeDataRoutine != "Toy"*/) {
      //std::cout << "Plotting fake data" << std::endl;
      std::vector<double> fake_xsec_vals;
      //std::cout << xs.size() << std::endl;
      for (size_t i = 0; i < xs.size(); ++i) {

        //std::cout << fFakeDataXSecs[it->first] << std::endl;
        //std::cout << "Testing" << std::endl;
        //std::cout << fFakeDataXSecs[it->first]->GetNbinsX() << std::endl;
        fake_xsec_vals.push_back(
            fFakeDataXSecs[it->first]->GetBinContent(i+1));
        //std::cout << "Adding " << fake_xsec_vals.back() << std::endl;
      }
      TGraph * fake_gr = new TGraph(xs.size(), &xs[0], &fake_xsec_vals[0]);
      fake_gr->Write();
      fake_gr->SetMarkerColor(kRed);
      fake_gr->SetMarkerStyle(20);
      fake_gr->Draw("same p");
      leg.AddEntry(fake_gr, "Fake Data", "p");
    }

    TString chi2_str;
    chi2_str.Form("Nominal #chi^{2} = %.2f", nominal_xsec_chi2);
    leg.AddEntry((TObject*)0x0, chi2_str, "");
    //std::string chi2_str = "Nominal #chi^{2} = " +
    //                       std::to_string(nominal_xsec_chi2);
    //leg.AddEntry((TObject*)0x0, chi2_str.c_str(), "");
    if (fDoFakeData /*&& fFakeDataRoutine != "Toy"*/) {
      //std::string fake_chi2_str = "Fake Data #chi^{2} = " +
      //                            std::to_string(fake_xsec_chi2);
      //leg.AddEntry((TObject*)0x0, fake_chi2_str.c_str(), "");
      TString fake_chi2_str;
      fake_chi2_str.Form("Fake Data #chi^{2} = %.2f", fake_xsec_chi2);
      leg.AddEntry((TObject*)0x0, fake_chi2_str, "");
    }
    leg.Draw("same");
    gPad->RedrawAxis();
    cThrow.Write();
    throw_gr.Write(gr_name.c_str());
    mean_gr.Write((gr_name + "_mean").c_str());
  }

  TVectorD fake_chi2_out(1);
  fake_chi2_out[0] = fake_xsec_chi2;
  fake_chi2_out.Write("FakeDataChi2");

  TVectorD nominal_chi2_out(1);
  nominal_chi2_out[0] = nominal_xsec_chi2;
  nominal_chi2_out.Write("NominalChi2");
}

void protoana::PDSPThinSliceFitter::GetCurrentTruthHists(
    std::map<int, std::vector<TH1*>> & throw_hists,
    std::map<int, std::vector<TH1*>> & throw_inc_hists,
    std::map<int, std::vector<TH1*>> & throw_xsec_hists) {
  //Loop over the samples
  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    //Get the number of bins from the first entry of the beam energy bins
    std::vector<std::vector<ThinSliceSample>> & samples_vec_2D = it->second;
    size_t nBins = samples_vec_2D[0].size();
    std::string name = it->second[0][0].GetName() + "Throw" +
                       std::to_string(throw_hists[it->first].size());
    TH1D * temp_hist = new TH1D(name.c_str(), "", nBins, 0, nBins); 
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {
      std::vector<ThinSliceSample> & samples_vec = samples_vec_2D[i];
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        temp_hist->AddBinContent(j+1, samples_vec[j].GetVariedFlux());
      }
    }
    throw_hists[it->first].push_back(temp_hist);
  }

  for (auto it = throw_inc_hists.begin(); it != throw_inc_hists.end(); ++it) {
    int s = it->first;
    auto & samples_vec_2D = fSamples[s];
    const std::vector<double> & bins = fSignalBins[s];
    std::string name = samples_vec_2D[0][0].GetName();
    name += "IncidentThrow" +
             std::to_string(throw_inc_hists[it->first].size());
    TH1D * temp_inc_hist = new TH1D(name.c_str(), "", bins.size() - 1, &bins[0]); 
    

    name = samples_vec_2D[0][0].GetName();
    name += "XSecThrow" +
             std::to_string(throw_inc_hists[it->first].size());
    TH1D * temp_xsec_hist = new TH1D(name.c_str(), "", bins.size() - 1,
                                     &bins[0]);
    for (auto i_s : fIncidentSamples) {
      auto & incident_vec_2D = fSamples[i_s];
      for (size_t i = 0; i < incident_vec_2D.size(); ++i) {
        for (size_t j = 0; j < incident_vec_2D[i].size(); ++j) {

            incident_vec_2D[i][j].FillHistFromIncidentEnergies(*temp_inc_hist);
        }
      }
    }
    throw_inc_hists[s].push_back(temp_inc_hist);

    for (int i = 1; i <= temp_xsec_hist->GetNbinsX(); ++i) {
      temp_xsec_hist->SetBinContent(
          i, throw_hists[s].back()->GetBinContent(i+1));
    }
    temp_xsec_hist->Divide(temp_inc_hist);
    throw_xsec_hists[s].push_back(temp_xsec_hist);
  }  
}




void protoana::PDSPThinSliceFitter::BuildFakeDataXSecs(bool use_scales) {
  //First, set all samples to fake data scales
  if (use_scales) {
    for (auto it = fFakeSamples.begin(); it != fFakeSamples.end(); ++it) {
      for (size_t i = 0; i < it->second.size(); ++i) {
        for (size_t j = 0; j < it->second[i].size(); ++j) {
        std::cout << "Fake scale " << it->first << " " << i << " " << j << " " << fFakeDataScales[it->first][j] << std::endl;
          it->second[i][j].SetFactorAndScale(fFakeDataScales[it->first][j]);
        }
      }
    }
  }

  std::map<int, TH1D*> fake_data_totals;

  for (auto s : fMeasurementSamples) {
    //auto & samples_vec_2D = fSamples[s];
    auto & samples_vec_2D = fFakeSamples[s];
    std::vector<double> bins = fSignalBins[s];
    std::string xsec_name = "FakeData" +
                             samples_vec_2D[0][1].GetName() + "XSec";
    if (fDrawXSecUnderflow) {
      bins.insert(bins.begin(), 0.);
      bins.push_back(2000.);
    }
    TH1D * temp_xsec = new TH1D(xsec_name.c_str(), "",
                                bins.size() - 1, &bins[0]);
    for (size_t i = 0; i < samples_vec_2D.size(); ++i) {

      size_t start = (fDrawXSecUnderflow ? 0 : 1);
      size_t end = samples_vec_2D[i].size() - (fDrawXSecUnderflow ? 0 : 1);
      for (size_t j = start; j < end; ++j) {
        int the_bin = (fDrawXSecUnderflow ? j+1 : j);
        temp_xsec->AddBinContent(the_bin, samples_vec_2D[i][j].GetVariedFlux());
      }
    }

    std::string inc_name = "FakeData" +
                             samples_vec_2D[0][1].GetName() + "Inc";
    TH1D * temp_inc = new TH1D(inc_name.c_str(), "", bins.size() - 1,
                               &bins[0]);
    for (size_t i = 0; i < fIncidentSamples.size(); ++i) {
      auto & vec_2D = fFakeSamples[fIncidentSamples[i]];
      for (size_t j = 0; j < vec_2D.size(); ++j) {
        auto & samples_vec = vec_2D[j];
        for (size_t k = 0; k < samples_vec.size(); ++k) {
            samples_vec[k].FillHistFromIncidentEnergies(*temp_inc);
        }
      }
    }

    std::string tot_name = "FakeData" +
                             samples_vec_2D[0][1].GetName() + "Tot";
    fake_data_totals[s] = (TH1D*)temp_xsec->Clone(tot_name.c_str());
    temp_xsec->Divide(temp_inc);

    CalculateCrossSection(temp_xsec);

    fFakeDataXSecs[s] = temp_xsec;
    fFakeDataIncs[s] = temp_inc;

    fFakeDataXSecs[s]->SetDirectory(0);
    fFakeDataIncs[s]->SetDirectory(0);
    fake_data_totals[s]->SetDirectory(0);

  }
  fOutputFile.cd();
  TDirectory * out = (TDirectory *)fOutputFile.mkdir("FakeDataXSecs");
  out->cd();
  for (auto it = fFakeDataXSecs.begin(); it != fFakeDataXSecs.end(); ++it) {
    it->second->Write();
    fFakeDataIncs[it->first]->Write();
    fake_data_totals[it->first]->Write();
  }

}

double protoana::PDSPThinSliceFitter::CalcChi2SystTerm() {
  std::vector<double> syst_vals, central_vals;
  for (auto it = fSystParameters.begin(); it != fSystParameters.end(); ++it) {
    syst_vals.push_back(it->second.GetValue());
    central_vals.push_back(it->second.GetCentral());
  }


  double new_result = 0.;
  for (size_t i = 0; i < syst_vals.size(); ++i) {
    int bin_1 = fCovarianceBinsSimple[i];
    for (size_t j = 0; j < syst_vals.size(); ++j) {
      int bin_2 = fCovarianceBinsSimple[j];
      new_result += ((syst_vals[i] - central_vals[i])*
                     (syst_vals[j] - central_vals[j])*
                     (*fCovMatrix)[bin_1][bin_2]);
    }
  }
  //std::cout << result << " " << new_result << std::endl;
  return new_result;
}

double protoana::PDSPThinSliceFitter::CalcSignalConstraint() {
  double result = 0.;
  for (auto & pars : fSignalParameters) {
    for (auto & par : pars.second) {
      result += std::pow((par - 1.)/fSignalConstraintSize, 2);
    }
  }

  return result;
}

//Trim TODO
void protoana::PDSPThinSliceFitter::MakeThrowsTree(TTree & tree, std::vector<double> & branches) {
  for (auto it = fSignalParameters.begin(); it != fSignalParameters.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      branches.push_back(0.);
      tree.Branch(fSignalParameterNames[it->first][i].c_str(),
                  &branches.back());
    }
 }

  for (auto it = fFluxParameters.begin(); it != fFluxParameters.end(); ++it) {
    branches.push_back(0.);
    tree.Branch(fFluxParameterNames[it->first].c_str(),
                &branches.back());
  }

  for (auto it = fSystParameters.begin(); it != fSystParameters.end(); ++it) {
    branches.push_back(0.);
    tree.Branch(it->second.GetName().c_str(), &branches.back());
  }

}


double protoana::PDSPThinSliceFitter::CalcRegTerm() {
  double result = 0.;
  for (auto it = fSignalParameters.begin(); it != fSignalParameters.end(); ++it) {
    for (size_t i = 0; i < it->second.size()-1; ++i) {
      result += std::pow((it->second[i] - it->second[i+1]), 2);
    }
  }

  return fRegFactor*result;
}

void protoana::PDSPThinSliceFitter::SetupTree() {
  fOutputTree = new TTree("tree", "");
  fOutputTree->Branch("chi2_stat", &fOutputChi2Stat);
  fOutputTree->Branch("chi2_syst", &fOutputChi2Syst);
  fOutputTree->Branch("par_vals", &fOutputParVals);
  fSaveFitTree = true;
}

void protoana::PDSPThinSliceFitter::WrapUpTree() {
  fOutputFile.cd();  
  fOutputTree->Write();
  fSaveFitTree = false;
}

void protoana::PDSPThinSliceFitter::GetFixFactors(
    /*std::map<int, std::vector<std::vector<ThinSliceSample>>> & samples,
    ThinSliceDataSet & data_set*/) {
  auto & sel_hists = fDataSet.GetSelectionHists();
//  std::map<int, TH1 *> fFixFactorHists;
  std::map<int, TH1 *> total_mc_hists;
  for (auto it = sel_hists.begin(); it != sel_hists.end(); ++it) {
    fFixFactorHists[it->first]
        = (TH1D *)it->second->Clone(("FixFactor" + std::to_string(it->first)).c_str());
    total_mc_hists[it->first]
        = (TH1D *)it->second->Clone(("TotalMC" + std::to_string(it->first)).c_str());
    total_mc_hists[it->first]->Reset();
  }

  for (auto it = fSamples.begin(); it != fSamples.end(); ++it) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      for (size_t j = 0; j < it->second[i].size(); ++j) {
        const auto & hists = it->second[i][j].GetSelectionHists();
        for (auto it2 = hists.begin(); it2 != hists.end(); ++it2) {
          total_mc_hists[it2->first]->Add(it2->second);
        }
      }
    }
  }

  std::cout << "Factors:" << std::endl;
  for (auto it = fFixFactorHists.begin(); it != fFixFactorHists.end(); ++it) {
    it->second->Divide(total_mc_hists[it->first]);
    for (int i = 1; i <= it->second->GetNbinsX(); ++i) {
      std::cout << it->second->GetBinContent(i) << " ";
    }
    std::cout << std::endl;
  }
}

void protoana::PDSPThinSliceFitter::GenerateCorrelatedThrow(
    const TH1D & pars, const TMatrixD * cov_lower/*TDecompChol & chol*/, std::vector<double> & vals) {
  bool rethrow = true;
  size_t nRethrows = 0;
  TVectorD rand(pars.GetNbinsX());
  while (rethrow && nRethrows < fMaxRethrows) {
    bool all_pos = true;


    //std::cout << "Throw Pars" << std::endl;
    for (int j = 1; j <= pars.GetNbinsX(); ++j) {
      if ((fThrowType == "NormalThrow") ||
          (j-1 >= fRemainCorrRange.first && j-1 <= fRemainCorrRange.second)) {
        rand[j-1] = fRNG.Gaus();
      }
      else {
        rand[j-1] = 0.;
      }
      //std::cout << j-1 << " " << rand[j-1] << std::endl;
    }

    //TVectorD rand_times_chol = chol.GetU()*rand;
    TVectorD rand_times_chol = (*cov_lower)*rand;

    for (int j = 1; j <= pars.GetNbinsX(); ++j) {
      //std::cout << rand_times_chol[j-1] << " ";
      vals[j-1] = pars.GetBinContent(j) + rand_times_chol[j-1];

      //Configure this to truncate at 0 or rethrow
      if (vals[j-1] < fParLimits[j-1]) {
        all_pos = false;
        std::cout << "Rethrowing " << j-1 << " " << vals[j-1] << " " <<
                     fParLimits[j-1] << std::endl;
      }
      if (vals[j-1] > fParLimitsUp[j-1]) {
        all_pos = false;
        std::cout << "Rethrowing " << j-1 << " " << vals[j-1] << " " <<
                     fParLimitsUp[j-1] << std::endl;
      }
      rethrow = !all_pos;
    }

    for (auto it = fFixSystsPostFit.begin(); it != fFixSystsPostFit.end();
         ++it) {
      vals[fSystParameterIndices[it->first]] = it->second;
    }

    //std::cout << std::endl;
    ++nRethrows;
  }
}

void protoana::PDSPThinSliceFitter::GenerateUncorrelatedThrow(
    const TH1D & pars, const TMatrixD * cov, std::vector<double> & vals) {

  for (int j = 1; j <= pars.GetNbinsX(); ++j) {

    bool rethrow = true;
    while (rethrow) {
      vals[j-1] = fRNG.Gaus(pars.GetBinContent(j), (*cov)[j-1][j-1]);

      //Configure this to truncate at 0 or rethrow
      if (vals[j-1] < fParLimits[j-1]) {
        std::cout << "Rethrowing " << j-1 << " " << vals[j-1] << " " <<
                     fParLimits[j-1] << std::endl;
      }
      else if (vals[j-1] > fParLimitsUp[j-1]) {
        std::cout << "Rethrowing " << j-1 << " " << vals[j-1] << " " <<
                     fParLimitsUp[j-1] << std::endl;
      }
      else {
        rethrow = false;
      }
    }
    std::cout << j << " " << vals[j-1] << " " << pars.GetBinContent(j) << " " <<
                 (*cov)[j-1][j-1] << std::endl;
  }

}

void protoana::PDSPThinSliceFitter::CalcFullCrossSection(TH1D * xsec_hist) {
  for (int i = 1; i <= xsec_hist->GetNbinsX(); ++i) {
    xsec_hist->SetBinContent(i, -1.*log(1. - xsec_hist->GetBinContent(i)));
  }
}

void protoana::PDSPThinSliceFitter::CalculateCrossSection(TH1D * xsec_hist) {
  CalcFullCrossSection(xsec_hist);

  xsec_hist->Scale(1.E27/ (fPitch * 1.4 * 6.022E23 / 39.948 ));
}

void protoana::PDSPThinSliceFitter::DoMinos1D() {
  fOutputFile.cd();
  fOutputFile.mkdir("Minos");
  fOutputFile.cd("Minos");

  size_t total_parameters = fTotalSignalParameters +
                            fTotalFluxParameters +
                            fTotalSystParameters +
                            fTotalG4RWParameters;
  TH1D errors_low("hMinosErrsLow", "",
                  total_parameters, 0, total_parameters),
       errors_up("hMinosErrsUp", "",
                  total_parameters, 0, total_parameters);
  for (size_t i = 0; i < total_parameters; ++i) {
    std::cout << "Running MINOS " << i << std::endl;
    double err_low = 0., err_up = 0.;
    fMinimizer->GetMinosError(i, err_low, err_up);
    std::cout << "\t" << err_low << " " << err_up << std::endl;
    errors_low.SetBinContent(i+1, err_low);
    errors_up.SetBinContent(i+1, err_up);
  }
  errors_low.Write();
  errors_up.Write();
  fOutputFile.cd();
}

void protoana::PDSPThinSliceFitter::GetCovarianceVals(TString dir) {
  fOutputFile.cd();
  fOutputFile.mkdir(dir);
  fOutputFile.cd(dir);

  size_t total_parameters = fTotalSignalParameters +
                            fTotalFluxParameters +
                            fTotalSystParameters +
                            fTotalG4RWParameters;
  TH2D covHist("covHist", "", total_parameters, 0,
               total_parameters, total_parameters, 0,
               total_parameters);
  TH2D corrHist("corrHist", "", total_parameters, 0,
                total_parameters, total_parameters, 0,
                total_parameters);
  TH1D varsHist("varsHist", "", total_parameters, 0,
                total_parameters);

  for (size_t i = 0; i < total_parameters; ++i) {
    varsHist.SetBinContent(i+1, fMinimizer->CovMatrix(i, i));
    for (size_t j = 0; j < total_parameters; ++j) {
      covHist.SetBinContent(i+1, j+1, fMinimizer->CovMatrix(i, j));
      corrHist.SetBinContent(i+1, j+1, fMinimizer->Correlation(i, j));
    }
  }
  covHist.Write();
  corrHist.Write();
  varsHist.Write();

  fOutputFile.cd();
}


void protoana::PDSPThinSliceFitter::DoMultiConts() {
  fOutputFile.cd();
  fOutputFile.mkdir("MultiConts");
  fOutputFile.cd("MultiConts");

  size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters;
  if (!fSignalContoursOnly) total_parameters += fTotalSystParameters;

  std::cout << "Running Contours" << std::endl;
  //Contour every other pair of parameters 
  std::vector<double> errors = {.5, 1., 2., 3.};
  for (size_t j = 0; j < 3; ++j) {
    fMinimizer->SetErrorDef(errors[j]);
    for (size_t i = 0; i < total_parameters; i += 2) {
      RunOneContour(i, i+1, fNContourPoints);
    }

    //If odd, contour last with second-to-last
    if (total_parameters % 2) {
      RunOneContour(total_parameters-2, total_parameters-1, fNContourPoints);
    }
  }

  fOutputFile.cd();
}

void protoana::PDSPThinSliceFitter::DoMinosConts() {
  fOutputFile.cd();
  fOutputFile.mkdir("MinosConts");
  fOutputFile.cd("MinosConts");

  size_t total_parameters = fTotalSignalParameters + fTotalFluxParameters;
  if (!fSignalContoursOnly) total_parameters += fTotalSystParameters;
      

  std::cout << "Running Contours" << std::endl;
  //Contour every other pair of parameters 
  /*unsigned int npoints = 10;*/
  for (size_t i = 0; i < total_parameters; i += 2) {
    RunOneContour(i, i+1, fNContourPoints);
  }

  //If odd, contour last with second-to-last
  if (total_parameters % 2) {
    RunOneContour(total_parameters-2, total_parameters-1, fNContourPoints);
  }

  fOutputFile.cd();
}

void protoana::PDSPThinSliceFitter::RunOneContour(size_t i, size_t j,
                                                  unsigned int & npoints) {
  std::cout << "\tPars: " << i << " " << j << std::endl;
  double * vals_i = new double[npoints], * vals_j = new double[npoints];

  fMinimizer->Contour(i, j, npoints, &vals_i[0], &vals_j[0]);
  TGraph gr_contour(npoints, &vals_i[0], &vals_j[0]);
  gr_contour.Write(TString::Format("grContour_%zu_%zu", i, j));

  delete[] vals_i;
  delete[] vals_j;
}

void protoana::PDSPThinSliceFitter::SaveHesse() {
  fOutputFile.cd();
  TVector output_hesse_status(1),
          output_post_hesse_status(1),
          output_cov_status(1);

  //output_hesse_status[0] = (hesse_good ? 1 : 0);
  output_hesse_status.Write("hesse_status");

  output_post_hesse_status[0] = fMinimizer->Status();
  output_post_hesse_status.Write("post_hesse_status");

  output_cov_status[0] = fMinimizer->CovMatrixStatus();
  output_cov_status.Write("post_hesse_cov_status");
}

bool protoana::PDSPThinSliceFitter::DoHesse() {
  fMinimizer->SetPrintLevel(10);
  std::cout << "Running Hesse" << std::endl;
  bool hesse_good = fMinimizer->Hesse();
  std::cout << "hesse good? " << hesse_good << std::endl;

  fMinimizer->SetPrintLevel(fPrintLevel);
  return hesse_good;
}