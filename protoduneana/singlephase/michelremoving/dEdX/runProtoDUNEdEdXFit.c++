#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "cetlib/filepath_maker.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSet.h"

#include "ProtoDUNEdEdXFitter.h"

#include "TFile.h"
#include "TROOT.h"


int main(int argc, char ** argv) {


  std::string fcl_file = "";
  std::string output_file = "";
  std::string input_file = "";
  // Options to run
  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-o")) {
      output_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-i")) {
      input_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: runPDSPThinSliceFit -c fclfile.fcl " << 
                    "-o outputfile.root -i input_file.root" << std::endl;
      return 1;
    }
  }


  gROOT->SetBatch();
  ProtoDUNEdEdXFitter fit(fcl_file, input_file, output_file);
  //fit.RunFit();
  fit.ParameterScans();

}
