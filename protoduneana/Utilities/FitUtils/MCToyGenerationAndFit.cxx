#include "MCToyGenerationAndFit.h"
#include "ProtoDUNEFitUtils.h"

// ROOT
#include "TROOT.h"
#include "TList.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include <Math/MinimizerOptions.h>

#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "RooAbsData.h"
#include "RooWorkspace.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooMultiVarGaussian.h"
#include "RooArgList.h"
#include "RooArgSet.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooFitLegacy/RooCatTypeLegacy.h"

//********************************************************************
protoana::MCToyGenerationAndFit::MCToyGenerationAndFit(){
  //********************************************************************

  _minimizer = "Minuit2";
  _algorithm = TString(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  _fitstrategy = 1;
  _minoserror = false;
  _conflevel = 0.95;

}

//********************************************************************
protoana::MCToyGenerationAndFit::MCToyGenerationAndFit(std::string minimizer, int fitstrategy, bool minoserror, double cl){
  //********************************************************************

  _minimizer = TString(minimizer.c_str());
  _algorithm = TString(ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  _fitstrategy = fitstrategy;
  _minoserror = minoserror;
  _conflevel = cl;

}

//********************************************************************
protoana::MCToyGenerationAndFit::~MCToyGenerationAndFit(){
  //********************************************************************
 
}

//********************************************************************
RooAbsData* protoana::MCToyGenerationAndFit::GenerateToyMC(RooWorkspace* ws, bool datanorm){
  //********************************************************************
  
  if(!ws){
    std::cerr << "ERROR::NULL workspace. No dataset generated!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace!" << std::endl;
    return NULL;
  }

  // Get the pdf
  RooAbsPdf* pdf = combined_config->GetPdf();
  if(!pdf){
    std::cerr << "ERROR::No pdf found in ModelConfig. No dataset generated!" << std::endl;
    return NULL;
  }

  RooAbsData* data = ws->data("obsData");
  if(!data){
    std::cerr << "ERROR::No dataset found in workspace. No dataset generated!" << std::endl;
    return NULL;
  }

  const RooArgSet* obs = combined_config->GetObservables();
  if(!obs){
    std::cerr << "ERROR::No observables found in ModelConfig. No dataset generated!" << std::endl;
    return NULL;
  }

  // Silence the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  RooAbsData* mctoy;
  if(datanorm){
    Int_t NEventsToGenerate = (Int_t)(data->sumEntries());
    mctoy = pdf->generate(*obs, RooFit::NumEvents(NEventsToGenerate), RooFit::AutoBinned(false));
  }
  else{
    mctoy = pdf->generate(*obs, RooFit::AutoBinned(false), RooFit::Extended());
  }

  // Reset verbosity
  RooMsgService::instance().reset();
  
  return mctoy;
  
}

//********************************************************************
TTree* protoana::MCToyGenerationAndFit::GenerateAndFit(RooWorkspace* ws, int nexp){
  //********************************************************************
  
  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }
  
  // Silence the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  RooAbsPdf* pdf = combined_config->GetPdf();
  const RooArgSet* obs_set = combined_config->GetObservables();

  // Create mc study. Only used to store the fit result
  RooMCStudy* mcstudy = new RooMCStudy( *pdf, *obs_set, RooFit::Save());

  // Create test statistics object
  RooStats::ProfileLikelihoodTestStat ts(*combined_config->GetPdf());

  // Create toy mc sampler
  RooStats::ToyMCSampler sampler(ts,nexp);
  sampler.SetPdf(*combined_config->GetPdf());
  sampler.SetObservables(*combined_config->GetObservables());
  sampler.SetGlobalObservables(*combined_config->GetGlobalObservables());
  sampler.SetParametersForTestStat(*combined_config->GetParametersOfInterest());
  
  RooArgSet poiAndNuisance;
  poiAndNuisance.add(*combined_config->GetParametersOfInterest());
  if(combined_config->GetNuisanceParameters())
    poiAndNuisance.add(*combined_config->GetNuisanceParameters());
  RooArgSet* nullParams = (RooArgSet*) poiAndNuisance.snapshot(); 
  //nullParams->Print();

  // Will be used as start values of fit
  ws->saveSnapshot("paramsToFit",poiAndNuisance);

  for(int i=0; i<nexp; ++i) {
    if(i%10==0)
      std::cout << "INFO::Running on toy " << i << std::endl;

    // Reset starting values of fit
    ws->loadSnapshot("paramsToFit");

    RooAbsData* toyMC = sampler.GenerateToyData( *nullParams );
    toyMC->Print();

    RooFitResult* fresult = combined_config->GetPdf()->fitTo(*toyMC, RooFit::Minimizer(_minimizer.Data(), _algorithm.Data()), RooFit::Save(), RooFit::Strategy(_fitstrategy), RooFit::Minos(_minoserror), RooFit::PrintLevel(-1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*toyMC)), RooFit::Extended(), RooFit::Offset(true));

    mcstudy->addFitResult(*fresult);
    
    delete toyMC;
  }

  // Reset verbosity
  RooMsgService::instance().reset();

  TTree* mcstree = RooMCStudyToTTree(mcstudy);

  delete mcstudy;

  return mcstree;

}

//********************************************************************
RooFitResult * protoana::MCToyGenerationAndFit::GenerateAndFitOneToy(
    RooWorkspace* ws) {
//********************************************************************
  
  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }
  
  // Silence the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  //RooAbsPdf* pdf = combined_config->GetPdf();
 // const RooArgSet* obs_set = combined_config->GetObservables();

  // Create mc study. Only used to store the fit result
  //RooMCStudy* mcstudy = new RooMCStudy( *pdf, *obs_set, RooFit::FitOptions("r"));

  // Create test statistics object
  RooStats::ProfileLikelihoodTestStat ts(*combined_config->GetPdf());

  // Create toy mc sampler
  RooStats::ToyMCSampler sampler(ts, 1);
  sampler.SetPdf(*combined_config->GetPdf());
  sampler.SetObservables(*combined_config->GetObservables());
  sampler.SetGlobalObservables(*combined_config->GetGlobalObservables());
  sampler.SetParametersForTestStat(*combined_config->GetParametersOfInterest());
  
  RooArgSet poiAndNuisance;
  poiAndNuisance.add(*combined_config->GetParametersOfInterest());
  if(combined_config->GetNuisanceParameters())
    poiAndNuisance.add(*combined_config->GetNuisanceParameters());
  RooArgSet* nullParams = (RooArgSet*) poiAndNuisance.snapshot(); 
  //nullParams->Print();

  // Will be used as start values of fit
  ws->saveSnapshot("paramsToFit",poiAndNuisance);

//  for(int i=0; i<nexp; ++i) {
//    if(i%10==0)
//      std::cout << "INFO::Running on toy " << i << std::endl;

    // Reset starting values of fit
    ws->loadSnapshot("paramsToFit");

    RooAbsData* toyMC = sampler.GenerateToyData( *nullParams );
    toyMC->Print();

    RooFitResult* fresult = combined_config->GetPdf()->fitTo(*toyMC, RooFit::Minimizer(_minimizer.Data(), _algorithm.Data()), RooFit::Save(), RooFit::Strategy(_fitstrategy), RooFit::Minos(_minoserror), RooFit::PrintLevel(-1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*toyMC)), RooFit::Extended(), RooFit::Offset(true));

    return fresult;

//    mcstudy->addFitResult(*fresult);
//    
//    delete toyMC;
//  }

}

//********************************************************************
RooFitResult* protoana::MCToyGenerationAndFit::FitData(RooWorkspace* ws, bool isWeighted){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }
  
  // Silence the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  // Dataset
  RooAbsData *obsdata = ws->data("obsData");
  
  RooFitResult* fresult = NULL;
  if(!isWeighted)
    fresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer(_minimizer.Data(), _algorithm.Data()), RooFit::Save(), RooFit::Strategy(_fitstrategy), RooFit::Minos(_minoserror), RooFit::PrintLevel(1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true) );
  else
    fresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer(_minimizer.Data(), _algorithm.Data()), RooFit::Save(), RooFit::Strategy(_fitstrategy), RooFit::Minos(_minoserror), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true), RooFit::SumW2Error(true) );

  // Reset the verbosity
  RooMsgService::instance().reset();

  return fresult;

}

//********************************************************************
RooFitResult* protoana::MCToyGenerationAndFit::FitAsimovData(RooWorkspace* ws){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }
 
  // Dataset
  RooAbsData *obsdata = ws->data("asimovData");

  RooFitResult* fresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer(_minimizer.Data(), _algorithm.Data()), RooFit::Save(), RooFit::Strategy(_fitstrategy), RooFit::Minos(_minoserror), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true));
  
  return fresult;

}

//********************************************************************
RooFitResult* protoana::MCToyGenerationAndFit::FitToyData(RooWorkspace* ws, RooAbsData* obsdata){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  RooFitResult* fresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer(_minimizer.Data(), _algorithm.Data()), RooFit::Save(), RooFit::Strategy(_fitstrategy), RooFit::Minos(_minoserror), RooFit::PrintLevel(1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true) );
  
  return fresult;

}

//********************************************************************
TTree* protoana::MCToyGenerationAndFit::FitToyMCFromWorkspace(std::vector<RooWorkspace*> wsvec, bool fitdata){
  //********************************************************************

  // List to store all trees
  TList *list = new TList;

  for(UInt_t j = 0 ; j < wsvec.size(); j++){
    RooWorkspace* ws = wsvec.at(j);
    if(fitdata){
      RooFitResult* result = FitData(ws);
      TTree *tree = RooFitResultToTTree(ws, result);
      if(tree)
	list->Add(tree);
    }
    else{
      RooDataSet* dataset = (RooDataSet*)GenerateToyMC(ws);
      RooFitResult* result = FitToyData(ws, dataset);
      TTree *tree = RooFitResultToTTree(ws, result);
      if(tree)
	list->Add(tree);
    }
  }

  TTree *newtree = TTree::MergeTrees(list); 
  newtree->SetName("MCToysFromWorkspace");

  return newtree;

}

//********************************************************************
TTree* protoana::MCToyGenerationAndFit::RooFitResultToTTree(RooWorkspace* ws, RooFitResult* res){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }
  
  RooAbsPdf* pdf = combined_config->GetPdf();
  const RooArgSet* obs_set = combined_config->GetObservables();

  // Create mc study. Only used to store the fit result
  RooMCStudy* mcstudy = new RooMCStudy( *pdf, *obs_set, RooFit::Save());
  mcstudy->addFitResult(*res);

  TTree* myTree = RooMCStudyToTTree(mcstudy);

  delete mcstudy;

  return myTree;

}

//********************************************************************
TTree* protoana::MCToyGenerationAndFit::RooMCStudyToTTree(RooMCStudy* mc){
  //********************************************************************

  TTree* myTree = new TTree("mcstree","mcstree");

  const RooDataSet& toydata = mc->fitParDataSet();

  // fitStatus = migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult
  // Minuit will return 0 if fit is good and 4 if fit failed (for example migrad did not converge)

  // Fit status for Migrad / improveResult (Migrad mininimization followed by an improved search for global minima)
  //status = 0 : OK
  //status = 1 : Covariance was made pos defined
  //status = 2 : Hesse is invalid
  //status = 3 : Edm is above max
  //status = 4 : Reached call limit
  //status = 5 : Any other failure

  // Fit status for Hesse
  // status = 0 : OK
  // status = 1 : hesse failed 
  // status = 2 : matrix inversion failed 
  // status = 3 : matrix is not pos defined

  // Fit status for Minos
  // status = 0 : OK
  // status = 1 : maximum number of function calls exceeded when running for lower error 
  // status = 2 : maximum number of function calls exceeded when running for upper error 
  // status = 3 : new minimum found when running for lower error 
  // status = 4 : new minimum found when running for upper error 
  // status = 5 : any other failure

  // Covariance matrix status (covQual)
  // -1 :  not available (inversion failed or Hesse failed)
  // 0  : available but not positive defined
  // 1  : covariance only approximate
  // 2  : full matrix but forced pos def
  // 3  : full accurate matrix

  Int_t covQual=0, status=0, numInvalidNLL=0;
  Float_t minNll=0, edm=0;
  myTree->Branch("covQual",        &covQual,        "covQual/I" );
  myTree->Branch("status",         &status,         "status/I" );
  myTree->Branch("numInvalidNLL",  &numInvalidNLL,  "numInvalidNLL/I" );
  myTree->Branch("minNll",         &minNll,         "minNll/F" );
  myTree->Branch("edm",            &edm,            "edm/F" );

  std::vector<Float_t> varVals;
  const RooArgSet* args = toydata.get();
  varVals.resize( args->getSize(), -999. );
  
  RooRealVar* var(0);
  TIterator* Itr = args->createIterator();
  for(Int_t i=0; (var = (RooRealVar*)Itr->Next()); ++i){
    TString varName = var->GetName();
    TString varNameF = TString(var->GetName()) + "/F";
    myTree->Branch( varName.Data(), &varVals[i], varNameF.Data() ); 
  }
  delete Itr;

  // Save correlations between parameters
  std::vector<Float_t> corrVals;
  const RooArgList& fitparams = mc->fitResult(0)->floatParsFinal();
  corrVals.resize( fitparams.getSize()*fitparams.getSize(), -999. );

  Int_t counter = 0;
  for(Int_t i=0; i < fitparams.getSize(); ++i){
    RooRealVar* fvar = (RooRealVar*)fitparams.at(i);
    TString varName = fvar->GetName();
    for(Int_t j=0; j < fitparams.getSize(); ++j){
      if(j == i) continue;
      RooRealVar* fvar2 = (RooRealVar*)fitparams.at(j);
      TString varName2 = varName + TString("_corr_") + TString(fvar2->GetName());
      TString varNameF = varName2 + "/F";
      myTree->Branch( varName2.Data(), &corrVals[counter], varNameF.Data() );
      counter++;
    }
  }  
  
  // Fill the tree by looping over mcstudy
  for(int iToy=0; iToy<toydata.numEntries(); iToy++){
    const RooFitResult* result = mc->fitResult(iToy);
    edm     = result->edm();
    covQual = result->covQual();
    status  = result->status();
    minNll  = result->minNll();
    numInvalidNLL = result->numInvalidNLL();
    
    toydata.get(iToy); // reset args to new value
    Itr = args->createIterator();
    for (Int_t i=0; (var=(RooRealVar*)Itr->Next()); ++i) { varVals[i] = var->getVal(); }
    delete Itr;

    const RooArgList& fparams = result->floatParsFinal();
    counter = 0;
    for(Int_t i=0; i < fparams.getSize(); ++i){
      RooRealVar* fvar = (RooRealVar*)fparams.at(i);
      for(Int_t j=0; j < fparams.getSize(); ++j){
	if(j == i) continue;
	RooRealVar* fvar2 = (RooRealVar*)fparams.at(j);
	corrVals[counter] = result->correlation(fvar->GetName(), fvar2->GetName());
	counter++;
      }
    }
    
    myTree->Fill();  
  }
  
  return myTree;

}

//********************************************************************
TTree* protoana::MCToyGenerationAndFit::RooDataSetToTTree(RooAbsData* data, TString treename){
  //********************************************************************

  if(!data){
    std::cerr << "ERROR::No RooAbsData. Will return NULL TTree!" << std::endl;
    return NULL;
  }

  TTree* myTree = new TTree(treename, treename);

  std::vector<Float_t> varVals;
  const RooArgSet* args = data->get();
  varVals.resize( args->getSize(), -999. );

  RooRealVar* var(0);
  TIterator* Itr = args->createIterator();
  for(Int_t i=0; (var = (RooRealVar*)Itr->Next()); ++i) {
    TString varName = var->GetName();
    TString varNameF = TString(var->GetName()) + "/F";
    myTree->Branch( varName.Data(), &varVals[i], varNameF.Data() );
  }
  delete Itr;

  // Fill the tree by looping over mcstudy
  for(int iToy=0; iToy<data->numEntries(); iToy++){
    data->get(iToy); // reset args to new value
    TIterator* Itr2 = args->createIterator();
    for (Int_t i=0; (var=(RooRealVar*)Itr2->Next()); ++i) {varVals[i] = var->getVal();}
    delete Itr2;

    myTree->Fill();
  }

  return myTree;

}

//********************************************************************
RooArgSet* protoana::MCToyGenerationAndFit::GetChi2Set(RooWorkspace *work, RooAbsData* data){
  //********************************************************************

  if(!work){
    std::cerr << "ERROR::NULL workspace. Will return empty vector!" << std::endl;
    return NULL;
  }

  // Get pdf from workspace
  RooSimultaneous* pdf = (RooSimultaneous*)work->pdf("simPdf");
  if(!pdf){
    std::cerr << "ERROR::No pdf found in workspace. Will return empty vector!" << std::endl;
    return NULL;
  }
  
   // If not provide, get data from workspace
  if(!data)
    data = work->data("obsData");

  // Get category components
  RooCategory* categories = work->cat("channelCat");

  // Silence output
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);

  std::vector<TString> categoriesName;  
  TIterator* iter = categories->typeIterator();
  RooCatType* catType;
  while( (catType = (RooCatType*) iter->Next())) {
    TString catname = catType->GetName();
    categoriesName.push_back(catname);
  }

  RooArgSet* chi2varset = new RooArgSet("chi2varset");
  for(UInt_t i = 0; i < categoriesName.size(); i++){
    TString catname = categoriesName[i];

    if(categories->setLabel(catname)){
      std::cout << "WARNING::Category " << catname.Data() << " is not a member of channelCat" << std::endl;
      continue;
    }
    
    RooAbsPdf* subpdf = (RooAbsPdf*)pdf->getPdf(catname.Data());
    if(!subpdf){
      std::cout << "WARNING::Can't find sub-pdf for region " << catname.Data() << std::endl;
      continue;
    }

    TString subdataset_str = Form("channelCat==channelCat::%s",catname.Data());
    RooAbsData* subdataset = (RooAbsData*) data->reduce(subdataset_str.Data());
    if(!subdataset){
      std::cout << "WARNING::Can't find sub-dataset for region " << catname.Data() << std::endl;
      continue;
    }

    RooRealVar* var =(RooRealVar*) ((RooArgSet*) subpdf->getObservables(*subdataset))->find(Form("obs_x_%s", catname.Data()));
 
    RooPlot* frame_dummy = var->frame();
    subdataset->plotOn(frame_dummy,RooFit::Cut(subdataset_str),RooFit::DataError(RooAbsData::Poisson));
    subpdf->plotOn(frame_dummy,RooFit::Normalization(1,RooAbsReal::RelativeExpected),RooFit::Precision(1e-5));

    RooRealVar* chi2var = new RooRealVar(Form("%s_var",subpdf->GetName()), Form("%s_var",subpdf->GetTitle()), frame_dummy->chiSquare());
    chi2varset->add(*chi2var);
  }

  return chi2varset;

}

//********************************************************************
TTree* protoana::MCToyGenerationAndFit::GenerateChi2Tree(RooWorkspace* ws, int nexp, int seed){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  RooAbsPdf* pdf = combined_config->GetPdf();
  if(!pdf){
    std::cerr << "ERROR::No pdf found in ModelConfig. No dataset generated!" << std::endl;
    return NULL;
  }

  const RooArgSet* obs = combined_config->GetObservables();
  const RooArgSet* globalobs = combined_config->GetGlobalObservables();
  if(!obs){
    std::cerr << "ERROR::No observables found in ModelConfig. No dataset generated!" << std::endl;
    return NULL;
  }

  const RooArgSet* pars = pdf->getParameters(obs);
  RooArgList floatParList;
  if(pars){
    TIterator* iter = pars->createIterator();
    RooAbsArg* arg;
    while((arg=(RooAbsArg*)iter->Next())){
      if(arg->InheritsFrom("RooRealVar") && !arg->isConstant()){
	floatParList.add(*arg);
      }
    }
    delete iter;
  }

  // Reset values and errors to nominal
  protoana::ProtoDUNEFitUtils::ResetValues(ws, floatParList);
  protoana::ProtoDUNEFitUtils::ResetError(ws, floatParList);
  protoana::ProtoDUNEFitUtils::ResetValuesToNominal(ws, *globalobs);
  
  // Silence the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  RooDataSet* chi2dataset = NULL;

  for(int i=0; i<nexp; ++i) {
    if(i%10==0)
      std::cout << "INFO::Generating toy " << i << std::endl;

    RooRandom::randomGenerator()->SetSeed(seed + i);

    RooAbsData* toyMC = pdf->generate(*obs, RooFit::AutoBinned(false), RooFit::Extended());
    toyMC->Print();

    RooFitResult* toyfitresult = FitToyData(ws,toyMC);
    
    RooArgSet* chi2set = GetChi2Set(ws,toyMC);
    if(i==0){
      chi2dataset = new RooDataSet("chi2dataset","chi2dataset",*chi2set);
    }

    // Only consider good quality fits
    if(toyfitresult->status() == 0)
      chi2dataset->add(*chi2set);

    protoana::ProtoDUNEFitUtils::ResetValues(ws, floatParList);
    protoana::ProtoDUNEFitUtils::ResetError(ws, floatParList);
    protoana::ProtoDUNEFitUtils::ResetValuesToNominal(ws, *globalobs);
  }

  // Reset verbosity
  RooMsgService::instance().reset();

  TTree* chi2_tree = RooDataSetToTTree(chi2dataset, "chi2tree");

  return chi2_tree;

}
