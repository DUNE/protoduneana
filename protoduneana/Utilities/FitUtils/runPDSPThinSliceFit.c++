#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "cetlib/filepath_maker.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"

#include "PDSPThinSliceFitter.h"

#include "TFile.h"

int main(int argc, char ** argv){

  std::string fcl_file;
  std::string output_file;
  std::string mc_file;
  std::string data_file;
  std::string refit_file = "";
  std::string tune_file = "";
  int njobs = -1;
  bool retune = false;
  // Options to run
  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-o")) {
      output_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-m")) {
      mc_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-d")) {
      data_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"--refit")) {
      refit_file = argv[++iArg]; 
    }
    if (!strcasecmp(argv[iArg],"--tune")) {
      tune_file = argv[++iArg]; 
    }
    if (!strcasecmp(argv[iArg],"--retune")) {
      retune = true; 
      ++iArg;
    }
    if (!strcasecmp(argv[iArg],"--njobs")) {
      njobs = std::atoi(argv[++iArg]);
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: runPDSPThinSliceFit -c fclfile.fcl " << 
                    "-o outputfile.root " << std::endl;
      return 1;
    }
  }

  protoana::PDSPThinSliceFitter * fit
      = new protoana::PDSPThinSliceFitter(fcl_file, output_file, mc_file,
                                          data_file, refit_file, tune_file,
                                          retune, njobs);
  //fit->InitializeMCSamples();
  //fit->Tune(refit_file);
  fit->BuildMCSamples();
  fit->RunFitAndSave();

  return 0;
}
