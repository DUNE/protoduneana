# Electromagnetic shower analysis (version : v10_04_04d01_e26_prof)
### using the branch : 
You will need to checkout ```dunesw```, ```dunecore```, ```duneprototypes```, ```duneopdet```, ```dunereco``` as well as this branch ```protoduneana/feature/chalamet_PDVD_EM_shower```, using the following command in ```$MRB_SOURCE``` folder : ```mrb g -t version dunesw```, version should be >v10_04_04d01_e26_prof.

### Goals :
- optimize the reconstructed energy of electromagnetic shower using a fiducial volume
- study the different backgrounds sources (cosmological and radiological)

### Reading the ouput of the LArSoft simulation and keeping what we want in another ROOT tree
```lar -c run_analyseEvents.fcl path/to/simulation/output ``` will use ```ParticleEvents_module.cc``` to create and populate a ROOT tree with the data we want to keep for the analysis. The output will be located in the current folder and named ```analysisOutput.root```.

### Using ROOT macros
In the Macros folder, run ```root -l -e '.L AnalysisModule/libAnalysisModule.so' macro_name.C``` specifying the path to ```analysisOutput.root``` in the macro you want to use.
