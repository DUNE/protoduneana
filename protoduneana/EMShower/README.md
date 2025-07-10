# Electromagnetic shower analysis (version : v10_04_04d01)
### Goals :
- optimize the reconstructed energy of electromagnetic shower using a fiducial volume
- study the different backgrounds sources (cosmological and radiological)

### Reading the ouput of the LArSoft simulation and keeping what we want in another ROOT tree
```lar -c run_analyseEvents.fcl path/to/simulation/output ``` will use ```ParticleEvents_module.cc``` to create and populate a ROOT tree with the data we want to keep for the analysis. The file is will be named ```analysisOutput.root``` in the current folder.

### Using ROOT macros
In the Macros folder, run ```root -l -e '.L AnalysisModule/libAnalysisModule.so' macro_name.C``` specifying the path to ```analysisOutput.root``` in the macro you want to use.
