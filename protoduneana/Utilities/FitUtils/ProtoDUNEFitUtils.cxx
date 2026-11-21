#include "ProtoDUNEFitUtils.h"

#include <TIterator.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TF1.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TKey.h>

#include <RooDataSet.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooRealIntegral.h>
#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooAbsArg.h>
#include <Roo1DTable.h>
#include <RooRealSumPdf.h>
#include <RooProduct.h>
#include <RooProdPdf.h>
#include <RooHist.h>

#include <RooStats/HistFactory/PiecewiseInterpolation.h>
#include <RooStats/ModelConfig.h>
#include "RooFitLegacy/RooCatTypeLegacy.h"

//********************************************************************
double GetArgonNumberDensity(double argon_density, double argon_molecularmass){
  //********************************************************************

  // argon_density: g/cm3, argon_molecularmass: g/mol
  double avogadro_number = 6.0221*10e23;
  return (avogadro_number*argon_density/argon_molecularmass);

}

//********************************************************************
std::vector<TH1*> protoana::ProtoDUNEFitUtils::GetSystHistograms(std::string name){
  //********************************************************************

  std::vector<TH1*> vec;

  TFile* f = new TFile(name.c_str(), "READ");
  if(!f) return vec;

  TIter next(f->GetListOfKeys());
  TKey *key;

  while ((key = (TKey*)next())){
    TString classname(key->GetClassName());
    if(!classname.Contains("TH1")) continue;

    TH1* h = (TH1*)key->ReadObj();
    h->SetDirectory(0);
    vec.push_back(h);
  }

  f->Close();

  return vec;

}

//********************************************************************
TH1* protoana::ProtoDUNEFitUtils::GetStatsSystHistogram(TH1* nominal){
  //********************************************************************

  TString name = TString(nominal->GetName()) + TString("_mcshapestats");
  TH1* histo = (TH1*)nominal->Clone(name.Data());

  for(Int_t i=1; i <= nominal->GetNbinsX(); i++){
    Double_t mcerror = sqrt(nominal->GetBinContent(i));

    Double_t ratio = 0.0;
    if(nominal->GetBinContent(i) > 0.)
      ratio = mcerror/nominal->GetBinContent(i);

    histo->SetBinContent(i, ratio);
    histo->SetBinError(i, 0.0);
  }

  return histo;

}

//********************************************************************
TH1* protoana::ProtoDUNEFitUtils::GetSystematicHistoFromNominal(TH1* nominal, TH1* syst, TString name){
  //********************************************************************

  TString histoname = TString(syst->GetName()) + TString("_") + name;
  TH1D* histo = new TH1D(histoname.Data(), histoname.Data(), nominal->GetNbinsX(), nominal->GetBinLowEdge(1), nominal->GetBinLowEdge(nominal->GetNbinsX()+1));
  if(nominal->GetNbinsX() != syst->GetNbinsX()){
    std::cout << "WARNING::Nominal and systematic histograms does not have the same number of bins. No systematic is applied!" << std::endl;
    return histo;
  }
 
  for(Int_t i=0; i < nominal->GetNbinsX(); i++){
    double nom = nominal->GetBinContent(i+1);
    double sys = syst->GetBinContent(i+1);

    if(nom <= 0.){
      std::cout << "WARNING::No events in histogram " << nominal->GetName() << " for bin " << i << ". Skip systematics propagation for systematic histogram " << syst->GetName() << std::endl;
      histo->SetBinContent(i+1, 0.0);
      continue;
    }

    if(sys < 0.)
      std::cout << "WARNING::Negative fractional error for histogram " << syst->GetName() << " and bin " << i << std::endl;

    // Avoid large fluctuation in the systematics
    if(sys > 1.0) sys = 1.0;

    double w = 1.0;
    if(name == "up" || name == "Up" || name == "UP" || name == "UP2" || name == "HIGH" || name == "High" || name == "high"){
      w = nom + nom*sys;
    }
    else if(name == "down" || name == "Down" || name == "DOWN" || name == "DOWN2" || name == "LOW" || name == "Low" || name == "low"){
      w = nom - nom*sys;
    }
 
    if(w < 0.) w = 0.;

    histo->SetBinContent(i+1, w);
  }

  return histo;

}

//********************************************************************
bool protoana::ProtoDUNEFitUtils::IsSingleBinHisto(TH1* histo){
  //********************************************************************

  int nbinswithevents = 0;
  for(int j = 0; j < histo->GetNbinsX(); j++){
    if(histo->GetBinContent(j) > 0.0)
      nbinswithevents++;
  }

  if(nbinswithevents == 1) return true;

  return false;

}

//********************************************************************
TH2* protoana::ProtoDUNEFitUtils::GetFitCovariance(RooFitResult* result){
  //********************************************************************

  RooArgList floatParsFinal = result->floatParsFinal();
  const TMatrixDSym& covmatrix = result->covarianceMatrix();
  Int_t n = covmatrix.GetNcols();
  TH2D* CovarianceHisto = new TH2D("FitCovariance", "FitCovariance", n,0,n, n,0,n);
  CovarianceHisto->GetXaxis()->SetLabelSize(0.01);
  CovarianceHisto->GetYaxis()->SetLabelSize(0.01);
  for(Int_t j = 0 ; j<n ; j++) {
    for(Int_t k = 0 ; k<n; k++) {
      CovarianceHisto->Fill(j+0.5,n-k-0.5, covmatrix(j,k));
    }
    CovarianceHisto->GetXaxis()->SetBinLabel(j+1, floatParsFinal.at(j)->GetName());
    CovarianceHisto->GetYaxis()->SetBinLabel(n-j, floatParsFinal.at(j)->GetName());
  }
  CovarianceHisto->SetMinimum(-1);
  CovarianceHisto->SetMaximum(+1);

  return CovarianceHisto;

}

//********************************************************************
TH2* protoana::ProtoDUNEFitUtils::GetFitCorrelationMatrix(RooFitResult* result){
  //********************************************************************

  TH2* corrhisto = result->correlationHist();
  corrhisto->GetXaxis()->SetLabelSize(0.01);
  corrhisto->GetYaxis()->SetLabelSize(0.01);

  return corrhisto;

}

//********************************************************************
void protoana::ProtoDUNEFitUtils::SaveSnapshot(RooWorkspace* ws, TString snapshotname){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::Workspace not found. Will not save snapshot." << std::endl;
    return;
  }

  RooSimultaneous* simpdf = (RooSimultaneous*)ws->pdf("simPdf");
  if(!simpdf)
    simpdf = (RooSimultaneous*)ws->pdf("combPdf");
  if(!simpdf){
    std::cerr << "ERROR::Pdf not found in workspace. Will not save snapshot." << std::endl;
    return;
  }

  RooAbsData* obsdata = (RooAbsData*)ws->data("obsData");
  RooArgSet* parameters = (RooArgSet*)simpdf->getParameters(*obsdata);
  if(!ws->loadSnapshot(snapshotname.Data())){
    ws->saveSnapshot(snapshotname.Data(),*parameters);
  } 
  else{
    std::cout << "WARNING::Snapshot " << snapshotname.Data() << " found in workspace. Will not overwrite it." << std::endl;
  }

  //gDirectory->Add(ws);

}

//********************************************************************
bool protoana::ProtoDUNEFitUtils::LoadSnapshot(RooWorkspace* ws, TString snapshotname){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::Workspace not found. Can't find snapshot." << std::endl;
    return false;
  }

  if(!ws->loadSnapshot(snapshotname)){
    std::cerr << "ERROR::Snapshot not found in workspace." << std::endl;
    return false;
  }

   ws->loadSnapshot(snapshotname);
   std::cout << "INFO::Workspace snapshot " << snapshotname.Data() << " loaded." << std::endl;

   return true;

}

//********************************************************************
void protoana::ProtoDUNEFitUtils::SaveWorkspace(RooWorkspace* ws, TString outFileName){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::Workspace not found. Won't save." << std::endl;
    return;
  }

  ws->writeToFile(outFileName.Data());
  std::cout << "INFO::Workspace written to file " << outFileName.Data() << std::endl;

  return;

}

//********************************************************************
void protoana::ProtoDUNEFitUtils::SetInterpolationCode(RooWorkspace* ws, Int_t code){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL workspace. No interpolation code added." << std::endl;
    return;
  }

  RooArgSet fun = ws->allFunctions();
  TIterator* iter = fun.createIterator();

  RooAbsArg* arg(0);
  while( (arg=(RooAbsArg*)iter->Next()) ) {
    if(arg->ClassName()!=TString("PiecewiseInterpolation") ) continue;
    PiecewiseInterpolation* p = (PiecewiseInterpolation*)ws->function(arg->GetName() );
    p->setAllInterpCodes(code);
  }

  delete iter;

}

//********************************************************************
void protoana::ProtoDUNEFitUtils::RemoveEmptyBins(RooPlot* frame){
  //********************************************************************

  const char* histname = 0;

  // Find histogram
  RooHist* histo = (RooHist*) frame->findObject(histname, RooHist::Class());
  if(!histo) return;

  for(Int_t i=0; i < histo->GetN(); i++){
    Double_t x,y;
    histo->GetPoint(i,x,y);

    if(fabs(y) < 0.00001 && histo->GetErrorYhigh(i) > 0.){
      histo->RemovePoint(i);
      if(i != histo->GetN()) --i;
    }
  }

  return;

}

//********************************************************************
double protoana::ProtoDUNEFitUtils::GetDataMCChi2(RooWorkspace *work, TString channelname, RooAbsData* data){
  //********************************************************************

  if(!work){
    std::cerr << "ERROR::NULL workspace. No chi2 is computed!" << std::endl;
    return -1.0;
  }

  // Silence output
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);

  // Get pdf from workspace
  RooSimultaneous* pdf = (RooSimultaneous*)work->pdf("simPdf");
  if(!pdf){
    std::cerr << "ERROR::No pdf found in workspace. No chi2 is computed!" << std::endl;
    return -1.0;
  }

  // If not provided, get data from workspace
  if(!data) data = work->data("obsData");

  // Get category components
  RooCategory* categories = work->cat("channelCat");

  std::vector<TString> categoriesName;  
  TIterator* iter = categories->typeIterator();
  RooCatType* catType;
  while( (catType = (RooCatType*) iter->Next())) {
    TString catname = catType->GetName();
    categoriesName.push_back(catname);
  }

  double chi2 = -999.0;
  for(UInt_t i = 0; i < categoriesName.size(); i++){
    TString catname = categoriesName[i];

    if(categories->setLabel(catname)){
      std::cout << "WARNING::Category " << catname.Data() << " is not a member of channelCat. Will skip." << std::endl;
      continue;
    }

    if(catname == channelname){

      RooAbsPdf* subpdf = (RooAbsPdf*)pdf->getPdf(catname.Data());
      if(!subpdf){
	std::cout << "WARNING::Can't find sub-pdf for region " << catname.Data() << ". Will skip." << std::endl;
	continue;
      }
      
      TString subdataset_str = Form("channelCat==channelCat::%s",catname.Data());
      RooAbsData* subdataset = (RooAbsData*) data->reduce(subdataset_str.Data());
      if(!subdataset){
	std::cout << "WARNING::Can't find sub-dataset for region " << catname.Data() << ". Will skip." << std::endl;
	continue;
      }
      
      RooRealVar* var =(RooRealVar*) ((RooArgSet*) subpdf->getObservables(*subdataset))->find(Form("obs_x_%s", catname.Data()));
      RooPlot* frame = var->frame();
      frame->SetName(Form("chi2frame_%s", catname.Data()));
      
      subdataset->plotOn(frame, RooFit::DataError(RooAbsData::Poisson), RooFit::MarkerColor(1), RooFit::LineColor(1));
      subpdf->plotOn(frame,RooFit::Normalization(1,RooAbsReal::RelativeExpected),RooFit::Precision(1e-5));
      chi2 = frame->chiSquare();
      break;
    }

  }

  return chi2;

}

std::vector<TH1 *> protoana::ProtoDUNEFitUtils::PlotXSecs(
        RooWorkspace * work, std::string name, /*std::string error,*/
        std::vector<TString> binnames, std::vector<double> recobins,
        std::vector<TString> incidentBinNames,
        RooAbsData * data,
        RooFitResult * result) {

  std::vector<TH1 *> xsecs;

  if (!work) {
    std::cerr << "ERROR:NULL dir. Will return empty canvas" << std::endl;
    return xsecs;
  }

  // Silence output
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);

  // Get pdf from workspace
  RooSimultaneous* pdf = (RooSimultaneous*)work->pdf("simPdf");
  if(!pdf){
    std::cerr << "ERROR::No pdf found in workspace. Will return empty vector!"
              << std::endl;
    std::cerr << "ERROR::No pdf found in workspace. Will return empty vector!"
              << std::endl;
    std::cerr << "ERROR::No pdf found in workspace. Will return empty vector!"
              << std::endl;
    return xsecs;
  }



  // Get category components
  // i.e. Incident, Abs, Cex
  RooCategory* categories = work->cat("channelCat");

  std::vector<TString> categoriesName;  
  TIterator* iter = categories->typeIterator();
  RooCatType* catType;
  while( (catType = (RooCatType*) iter->Next())) {
    TString catname = catType->GetName();
    categoriesName.push_back(catname);
    std::cout << catname << std::endl;
  }

  for (size_t i = 0; i < categoriesName.size(); ++i) {

    TString catname = categoriesName[i];

    RooAbsPdf* subpdf = (RooAbsPdf*)pdf->getPdf(catname.Data());
    if(!subpdf){
      std::cout << "WARNING::Can't find sub-pdf for region " << catname.Data()
                << ". Will skip." << std::endl;
      continue;
    }

    TString RRSumPdfName = Form("%s_model",catname.Data()); 
    RooRealSumPdf* RRSumPdf =
        (RooRealSumPdf*)subpdf->getComponents()->find(RRSumPdfName);
    RooArgList RRSumComponentsList =  RRSumPdf->funcList();
    RooLinkedListIter iter = RRSumComponentsList.iterator();
    RooProduct* component;

    std::vector<TString> compNameVec;
    while( (component = (RooProduct*) iter.Next()) ){
      TString componentName = component->GetName();
      std::cout << componentName << std::endl;
    }
  }

  /*
  TH1 * incident_signal = 0x0;
  std::vector<TH1 *> incident_backgrounds;

  for (const auto  && key : *keys) {
    std::string name = key->GetName();

    auto find_Incident = name.find("MC_ChannelIncident_Pions");
    if (find_Incident != std::string::npos) {
      incident_signal = (TH1*)key->Clone();
    }
    else {
      incident_backgrounds.push_back((TH1*)key->Clone());
    }
  }
  std::cout << incident_signal << " " << incident_backgrounds.size() << std::endl;

  for (const std::string & main_channel : {"ABS", "CEX"}) {
    std::cout << main_channel << std::endl;

    for (const auto  && key : *keys) {
       std::string name = key->GetName();
       std::cout << name << std::endl;
    }
  }
  */

  return xsecs;
}

//********************************************************************
std::vector<TCanvas*> protoana::ProtoDUNEFitUtils::PlotDatasetsAndPdfs(
    RooWorkspace *work, TString name, TString error, TString plottodraw,
    std::vector<TString> binnames, std::vector<double> recobins,
    std::vector<TString> incidentBinNames,
    std::vector<TString> sidebandBinNames, std::vector<double> sidebandBins,
    TString measurement,
    bool doNegativeReco, RooAbsData* data, RooFitResult* result) {
//********************************************************************

  std::vector<TCanvas*> rooplots;

  if(!work){
    std::cerr << "ERROR::NULL workspace. Will return empty vector!" << std::endl;
    return rooplots;
  }

  // Silence output
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);

  // Get pdf from workspace
  RooSimultaneous* pdf = (RooSimultaneous*)work->pdf("simPdf");
  if(!pdf){
    std::cerr << "ERROR::No pdf found in workspace. Will return empty vector!" << std::endl;
    std::cerr << "ERROR::No pdf found in workspace. Will return empty vector!" << std::endl;
    std::cerr << "ERROR::No pdf found in workspace. Will return empty vector!" << std::endl;
    return rooplots;
  }

  // If not provide, get data from workspace
  if(!data) data = work->data("obsData");

  // Get category components
  RooCategory* categories = work->cat("channelCat");

  // Print table with the number of data entries per category (channel)
  data->table(*((RooAbsCategory*)categories))->Print("v");

  std::vector<TString> categoriesName;  
  TIterator* iter = categories->typeIterator();
  RooCatType* catType;
  while( (catType = (RooCatType*) iter->Next())) {
    TString catname = catType->GetName();
    categoriesName.push_back(catname);
  }

  for(UInt_t i = 0; i < categoriesName.size(); i++){
    TString catname = categoriesName[i];

    if(categories->setLabel(catname)){
      std::cout << "WARNING::Category " << catname.Data() << " is not a member of channelCat. Will skip." << std::endl;
      continue;
    }
    
    RooAbsPdf* subpdf = (RooAbsPdf*)pdf->getPdf(catname.Data());
    if(!subpdf){
      std::cout << "WARNING::Can't find sub-pdf for region " << catname.Data() << ". Will skip." << std::endl;
      continue;
    }

    TString subdataset_str = Form("channelCat==channelCat::%s",catname.Data());
    RooAbsData* subdataset = (RooAbsData*) data->reduce(subdataset_str.Data());
    if(!subdataset){
      std::cout << "WARNING::Can't find sub-dataset for region " << catname.Data() << ". Will skip." << std::endl;
      continue;
    }

    RooRealVar* var =(RooRealVar*) ((RooArgSet*) subpdf->getObservables(*subdataset))->find(Form("obs_x_%s", catname.Data()));
 
    // Define canvas first
    TString canName = Form("canvas_%s_%s_%s", catname.Data(), name.Data(), plottodraw.Data());
    TCanvas* c = new TCanvas(canName, canName);

    // Legend
    TLegend* legend = new TLegend(0.10,0.70,0.90,0.90,"");
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);

    legend->SetNColumns(3);

    TLegendEntry* entry = legend->AddEntry("","Data","pl") ;
    entry->SetMarkerColor(kBlack);
    entry->SetMarkerStyle(20);
    entry=legend->AddEntry("","MC Total","lf") ;
    entry->SetLineColor(kBlack);
    entry->SetFillColor(kBlue);
    entry->SetFillStyle(3444);

    RooPlot* frame = var->frame();
    frame->SetName(Form("frame_%s_%s_%s", catname.Data(), name.Data(), plottodraw.Data()));

    // Draw data and pdf on top of error band - change to RooAbsData::SumW2 if data is weighted
    if(error.Contains("SumW2") || error.Contains("sumw2") || error.Contains("sumW2"))
      subdataset->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::MarkerColor(1), RooFit::LineColor(1));
    else
      subdataset->plotOn(frame, RooFit::DataError(RooAbsData::Poisson), RooFit::MarkerColor(1), RooFit::LineColor(1));

    TString RRSumPdfName = Form("%s_model",catname.Data()); 
    RooRealSumPdf* RRSumPdf = (RooRealSumPdf*) subpdf->getComponents()->find(RRSumPdfName);

    RooArgList RRSumComponentsList =  RRSumPdf->funcList();

    TString binWidthName =  Form("binWidth_obs_x_%s_0",catname.Data());
    RooRealVar* regionBinWidth = ((RooRealVar*) RRSumPdf->getVariables()->find(Form("binWidth_obs_x_%s_0",catname.Data()))) ;
    if(regionBinWidth == NULL)
      std::cout << "WARNING::bindWidth variable not found for region(" << catname << "). PLOTTING COMPONENTS WILL BE WRONG!" << std::endl;

    // Normalize data to expected number of events 
    double normCount = subpdf->expectedEvents(*var);
    if(result)
      std::cout << "INFO::MC events after fit = " << normCount << " for " << subpdf->GetName() << std::endl;

    RooLinkedListIter iter = RRSumComponentsList.iterator() ;
    RooProduct* component;

    std::vector<TString> compNameVec;
    std::vector <double> compFracVec;
    std::vector<TString> compStackNameVec;
    std::vector <double> compStackFracVec;
    compNameVec.clear();
    compStackNameVec.clear();
    compFracVec.clear();
    compStackFracVec.clear();

    while( (component = (RooProduct*) iter.Next()) ){
      TString  componentName = component->GetName();
      TString stackComponentName = componentName; 
      
      if(!compStackNameVec.empty())
	stackComponentName  = Form("%s,%s",compStackNameVec.back().Data() ,componentName.Data());
      compNameVec.push_back(componentName);
      compStackNameVec.push_back(stackComponentName);

      RooAbsReal*  i_RRSumPdf = ((RooAbsPdf*)work->pdf(RRSumPdfName))->createIntegral(RooArgSet(*var));
      Double_t Int_RRSumPdf = i_RRSumPdf->getVal();
      RooAbsReal*  i_component =   ((RooProduct*)work->obj(componentName))->createIntegral(RooArgSet(*var));
      Double_t Int_component = i_component->getVal();

      Double_t componentFrac = 0.;
      if(Int_RRSumPdf != 0.) 
	componentFrac =  Int_component * regionBinWidth->getVal() / Int_RRSumPdf;
      
      double stackComponentFrac = componentFrac; 
      if(!compStackFracVec.empty())
	stackComponentFrac  = compStackFracVec.back() + componentFrac; 
      
      compFracVec.push_back(componentFrac);
      compStackFracVec.push_back(stackComponentFrac);
    }

    Int_t counter = 0;
    Int_t sigcolor[13] = {2,3,4,5,6,7,8,9,kMagenta, 1, kGreen+2, kTeal, kOrange+10};
    for(int i = (compFracVec.size()-1); i > -1; i--){
      Int_t compPlotColor = i;
      if(compNameVec[i].Contains("ChannelABS_CEX") || compNameVec[i].Contains("ChannelCEX_ABS")){
	compPlotColor = 0;
      }
      else if(compNameVec[i].Contains("RecoBin")){
	compPlotColor = 40;
      }
      else{
	compPlotColor = sigcolor[counter];
	counter++;
      }
      std::cout << "INFO::Drawing " << compNameVec[i] << " " << counter << " " << compPlotColor << std::endl;
     
      subpdf->plotOn(frame,RooFit::Components(compStackNameVec[i].Data()),RooFit::FillColor(compPlotColor),RooFit::FillStyle(3001),RooFit::DrawOption("F"),RooFit::Normalization(compStackFracVec[i]*normCount,RooAbsReal::NumEvent),RooFit::Precision(1e-5));
    }

    if(result)
      subpdf->plotOn(frame, RooFit::Normalization(1,RooAbsReal::RelativeExpected), RooFit::Precision(1e-5), RooFit::VisualizeError(*result), RooFit::FillColor(kBlue), RooFit::FillStyle(3444));

    // Plot again so that it is on top of errors
    subpdf->plotOn(frame, RooFit::Normalization(1,RooAbsReal::RelativeExpected), RooFit::Precision(1e-5), RooFit::LineColor(1));
    if(error.Contains("SumW2") || error.Contains("sumw2") || error.Contains("sumW2"))
      subdataset->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::MarkerColor(1), RooFit::LineColor(1));
    else
      subdataset->plotOn(frame, RooFit::DataError(RooAbsData::Poisson), RooFit::MarkerColor(1), RooFit::LineColor(1));

    // Remove empty data bins
    RemoveEmptyBins(frame);
    
    Int_t counter2 = 0;
    bool found = false; bool found2 = false;
    //for(int i = (compNameVec.size()-1) ; i > -1; i--)
    for(unsigned int i = 0; i < compFracVec.size(); i++){
      std::cout << "NAME: " << compNameVec[i] << std::endl;
      Int_t compPlotColor = i;
      if(compNameVec[i].Contains("ChannelABS_CEX") && !found){
	compPlotColor = 0;
	entry=legend->AddEntry("","CEX Bkg","f");
	entry->SetLineColor(46);
	entry->SetFillColor(compPlotColor);
	entry->SetFillStyle(3001);
	found = true;
	continue;
      }
      else if(compNameVec[i].Contains("ChannelCEX_ABS") && !found){
	compPlotColor = 0;
	entry=legend->AddEntry("","ABS Bkg","f");
	entry->SetLineColor(46);
	entry->SetFillColor(compPlotColor);
	entry->SetFillStyle(3001);
	found = true;
	continue;
      }
      else if(compNameVec[i].Contains("RecoBin") && !found2){
	compPlotColor = 40;
	entry=legend->AddEntry("","Signal","f");
	entry->SetLineColor(40);
	entry->SetFillColor(compPlotColor);
	entry->SetFillStyle(3001);
	found2 = true;
	continue;
      }
      else{
	counter--;
	compPlotColor = sigcolor[counter];	
      }

      if(counter2 < (int)binnames.size()){
	TString legName = binnames[counter2];
	if (compNameVec[i].Contains("Incident"))
          legName = incidentBinNames[counter2];
        if (compNameVec[i].Contains("Sideband"))
          legName = sidebandBinNames[counter2];
	
	entry=legend->AddEntry("",legName.Data(),"f");
	entry->SetLineColor(compPlotColor);
	entry->SetFillColor(compPlotColor);
	entry->SetFillStyle(3001);
	counter2++;
      }
    }

    // two pads, one for standard plot, one for data/MC ratio
    float yMinP1=0.305;
    float bottomMarginP1=0.015;
    TPad *pad1 = new TPad(Form("%s_pad1",canName.Data()),Form("%s_pad1",canName.Data()),0.,yMinP1,.99,1);
    pad1->SetBottomMargin(bottomMarginP1);
    pad1->SetFillColor(kWhite);
    pad1->SetTickx();
    pad1->SetTicky();
    TPad *pad2 = new TPad(Form("%s_pad2",canName.Data()),Form("%s_pad2",canName.Data()),0.,0.01,.99,0.295);
    pad2->SetTopMargin(0.021);
    pad2->SetBottomMargin(0.3);
    pad2->SetFillColor(kWhite);
    
    pad1->Draw();
    pad2->Draw(); 
    frame->GetXaxis()->SetLabelSize(0.);
    
    pad1->cd();
    frame->SetTitle(measurement.Data()); 
    frame->Draw(); //Draw("same");
    frame->GetXaxis()->SetTitle("Reco Bin");
    frame->GetYaxis()->SetTitle("Events");


    std::string chi2_text = "#chi^{2}: " +
                            std::to_string(GetDataMCChi2(work, catname, data));
    legend->AddEntry((TObject*)0x0, chi2_text.c_str(), "");

    legend->Draw();

    pad2->cd();
    RooPlot* frame_dummy = var->frame();
 
    subdataset->plotOn(frame_dummy,RooFit::Cut(subdataset_str),RooFit::DataError(RooAbsData::Poisson));

    // Get RooHist of the data - will be used later in the ratio plot
    const char* curvename = 0;
    RooHist* nominal_hist = (RooHist*) frame_dummy->findObject(curvename, RooHist::Class());
    if(!nominal_hist){
      std::cerr << "ERROR::Drawing the data/MC histogram. Can't find nominal histogram." << std::endl;
      return rooplots;
    }

    // Normalize pdf to number of expected events, not to number of events in dataset
    subpdf->plotOn(frame_dummy,RooFit::Normalization(1,RooAbsReal::RelativeExpected),RooFit::Precision(1e-5));

    // Get RooCurve of the pdf - will be used later in the ratio plot
    RooCurve* nominal_curve = (RooCurve*) frame_dummy->findObject(curvename, RooCurve::Class());
    if(!nominal_curve){
      std::cerr << "ERROR::Drawing the data/MC histogram. Can't find nominal curve." << std::endl;
      return rooplots;
    }

    // frame_dummy->Print("v");

    RooHist* hratio = NULL;
    RooCurve* ratio_curve = new RooCurve;

    if(plottodraw == "pull"){
      hratio = (RooHist*) frame_dummy->pullHist();
      hratio->SetTitle("Pull Distribution");
    }
    else if(plottodraw == "residuals"){
      hratio = (RooHist*) frame_dummy->residHist();
      hratio->SetTitle("Residual Distribution");
    }
    else if(plottodraw == "ratio"){
      if(result)
	subpdf->plotOn(frame_dummy, RooFit::Normalization(1,RooAbsReal::RelativeExpected), RooFit::Precision(1e-5), RooFit::VisualizeError(*result, 1), RooFit::FillColor(kBlue-5), RooFit::FillStyle(3004));

      // Get error plot
      RooCurve* error_curve = (RooCurve*)frame_dummy->findObject(curvename, RooCurve::Class());
      if(!error_curve){
	std::cerr << "ERROR::Drawing the data/MC histogram. Can't find error curve." << std::endl;
	return rooplots;
      }

      ratio_curve->SetName(Form("%s_ratiocurve",nominal_curve->GetName()));
      ratio_curve->SetLineColor(kBlue-5);
      ratio_curve->SetFillColor(kBlue-5);
      ratio_curve->SetFillStyle(3004);
      ratio_curve->SetLineWidth(1);
      
      // Fill error curve
      Int_t j = 0;
      bool bottomCurve = false;
      for(Int_t i=1; i < error_curve->GetN()-1; i++){
        Double_t x = 0.;
        Double_t y = 0.;
        error_curve->GetPoint(i,x,y) ;

	if( i >= (nominal_curve->GetN()-1) ) bottomCurve = true;

        Double_t xNom = x;
        Double_t yNom = y;

        if( i == (nominal_curve->GetN() - 1) ||  i == nominal_curve->GetN() ){
	  ratio_curve->addPoint(x, 0.);   
	  continue;
        }

        if(bottomCurve){
	  nominal_curve->GetPoint(j,xNom,yNom);
	  j--;
        } 
	else{
	  j++;
	  nominal_curve->GetPoint(j,xNom,yNom);
        }

        if(fabs(yNom) > 0.00001){ 
	  ratio_curve->addPoint(x, (y / yNom));
        } 
	else{ 
	  ratio_curve->addPoint(x, 0.); 
        }
      }

      // Define ratio plot
      hratio = new RooHist(nominal_hist->getNominalBinWidth());
      
      // Determine range of curve
      Double_t xstart, xstop, y;
      nominal_curve->GetPoint(2,xstart,y);
      nominal_curve->GetPoint(nominal_curve->GetN()-1,xstop,y);
     
      for(Int_t i=0; i < nominal_hist->GetN(); i++){
	Double_t x,ypoint;
	nominal_hist->GetPoint(i,x,ypoint);
	
	if(x < xstart || x > xstop) continue;
	if(fabs(ypoint) < 0.000001 ) continue;

	Double_t yy;
	Double_t yerrorl = nominal_hist->GetErrorYlow(i);
	Double_t yerrorh = nominal_hist->GetErrorYhigh(i);

	//Double_t xerrorl = nominal_curve->GetErrorXlow(i);
	//Double_t xerrorh = nominal_curve->GetErrorXhigh(i);
	//if(xerrorl<=0 ) xerrorl = nominal_curve->GetErrorX(i);
	//if(xerrorh<=0 ) xerrorh = nominal_curve->GetErrorX(i);
	//if(xerrorl<=0 ) xerrorl = 0.5*nominal_hist->getNominalBinWidth();
	//if(yerrorh<=0 ) yerrorh = 0.5*nominal_hist->getNominalBinWidth();

	//yy = ypoint / nominal_curve->average(x-exl,x+exh);
	//yerrorl /= nominal_curve->average(x-xerrorl,x+xerrorh);
	//yerrorh /= nominal_curve->average(x-xerrorl,x+xerrorh);
	
	yy = ypoint / nominal_curve->interpolate(x);
	yerrorl /= nominal_curve->interpolate(x);
	yerrorh /= nominal_curve->interpolate(x);

	hratio->addBinWithError(x,yy,yerrorl,yerrorh);
      }
    }

    if(!hratio){
      std::cerr << "ERROR::Drawing data and MC histogram failed. RooHist is not found." << std::endl;
      return rooplots;
    } 

    hratio->SetMarkerColor(kRed);
    hratio->SetLineColor(kRed);

    // Create a new frame to draw the residual/pull/ratio distribution and add the distribution to the frame
    RooPlot* frame2 = var->frame();
    if(plottodraw == "ratio" && result) frame2->addPlotable(ratio_curve,"F");
    frame2->addPlotable(hratio,"P");
    
    if (catname.Contains("Sideband")) {
      //frame2->GetXaxis()->SetBinLabel(1, "Generic Label");
      frame2->GetXaxis()->SetTitle("Generic Label");
      for(size_t i = 1; i < sidebandBins.size(); i++){
        TString ibinstr = Form("%.1f-%.1f",sidebandBins[i-1],sidebandBins[i]);
        frame2->GetXaxis()->SetBinLabel(i, ibinstr.Data());
      }
    }
    else {
      if (doNegativeReco) {
        frame2->GetXaxis()->SetBinLabel(1, "< 0.");
        for(size_t i = 1; i < recobins.size(); i++){
          TString ibinstr = Form("%.1f-%.1f",recobins[i-1],recobins[i]);
          frame2->GetXaxis()->SetBinLabel(i+1, ibinstr.Data());
        }
      }
      else {
        for(size_t i = 1; i < recobins.size(); i++){
          TString ibinstr = Form("%.1f-%.1f",recobins[i-1],recobins[i]);
          frame2->GetXaxis()->SetBinLabel(i, ibinstr.Data());
        }
      }
      frame2->GetXaxis()->SetTitle("E_{reco} [MeV]");
    }

    // Cosmetics
    int firstbin = frame_dummy->GetXaxis()->GetFirst();
    int lastbin = frame_dummy->GetXaxis()->GetLast();
    double xmax = frame_dummy->GetXaxis()->GetBinUpEdge(lastbin);
    double xmin = frame_dummy->GetXaxis()->GetBinLowEdge(firstbin);
    
    if(plottodraw=="pull"){
      TLine* lp1 = new TLine(xmin,1.,xmax,1.);	
      TLine* lp2 = new TLine(xmin,2.,xmax,2.);	
      TLine* lp3 = new TLine(xmin,3.,xmax,3.);
      TLine* lp4 = new TLine(xmin,4.,xmax,4.);
      
      TLine* lp6 = new TLine(xmin,-1.,xmax,-1.);	
      TLine* lp7 = new TLine(xmin,-2.,xmax,-2.);	
      TLine* lp8 = new TLine(xmin,-3.,xmax,-3.);
      TLine* lp9 = new TLine(xmin,-4.,xmax,-4.);

      TLine* lp10 = new TLine(xmin,0,xmax,0);
      
      lp1->SetLineStyle(3);
      lp2->SetLineStyle(3);
      lp3->SetLineStyle(3);
      lp4->SetLineStyle(3);      
      lp6->SetLineStyle(3);
      lp7->SetLineStyle(3);
      lp8->SetLineStyle(3);
      lp9->SetLineStyle(3);
      lp10->SetLineStyle(3);
      
      frame2->addObject(lp1);
      frame2->addObject(lp2);
      frame2->addObject(lp3);
      frame2->addObject(lp4);
      frame2->addObject(lp6);
      frame2->addObject(lp7);
      frame2->addObject(lp8);
      frame2->addObject(lp9);
      frame2->addObject(lp10);

      frame2->SetMinimum(-4.0);
      frame2->SetMaximum(4.0);

      frame2->GetYaxis()->SetTitle("Pull");
    }
    else if(plottodraw == "residuals"){
      TLine* l = new TLine(xmin,0.,xmax,0.);
      l->SetLineWidth(1);
      l->SetLineStyle(2);
      frame2->addObject(l);
      frame2->GetYaxis()->SetTitle("Residual");
    }
    else if(plottodraw == "ratio"){
      TLine* lp1 = new TLine(xmin,1.,xmax,1.);
      TLine* lp2 = new TLine(xmin,0.5,xmax,0.5);
      TLine* lp3 = new TLine(xmin,1.5,xmax,1.5);
      TLine* lp4 = new TLine(xmin,2.,xmax,2.);
      TLine* lp5 = new TLine(xmin,2.5,xmax,2.5);
      lp1->SetLineWidth(1);
      lp1->SetLineStyle(3);
      lp2->SetLineStyle(3);
      lp3->SetLineStyle(3);
      lp4->SetLineStyle(3);
      lp5->SetLineStyle(3);
      frame2->addObject(lp1);
      frame2->addObject(lp2);
      frame2->addObject(lp3);
      frame2->addObject(lp4);
      frame2->addObject(lp5);

      frame2->SetMinimum(.0);
      frame2->SetMaximum(3.0);

      frame2->GetYaxis()->SetTitle("Data / MC");
    }
 
    frame2->GetYaxis()->SetLabelSize(0.10);
    frame2->GetYaxis()->SetNdivisions(504);         
    frame2->GetXaxis()->SetLabelSize(0.10);
    frame2->GetYaxis()->SetTitleSize(0.10);
    frame2->GetXaxis()->SetTitleSize(0.10);
    frame2->GetYaxis()->SetTitleOffset(0.35);
    frame2->GetXaxis()->SetTitleOffset(1.);
    frame2->GetYaxis()->SetLabelOffset(0.01);
    frame2->GetXaxis()->SetLabelOffset(0.03);
    frame2->GetXaxis()->SetTickLength(0.06);
    
    frame2->SetTitle("");
    frame2->GetYaxis()->CenterTitle(); 
    frame2->Draw();
		
    rooplots.push_back(c);
  }

  return rooplots;

}

//********************************************************************
std::vector<TCanvas*> protoana::ProtoDUNEFitUtils::PlotNLL(RooWorkspace *work, TString name, RooFitResult* result, bool plotPLL){
  //********************************************************************

  std::vector<TCanvas*> rooplots;

  if(!work){
    std::cerr << "ERROR::NULL workspace. Will return empty vector!" << std::endl;
    return rooplots;
  }

  if(!result){
    std::cerr << "ERROR::NULL fit result. Will return empty vector!" << std::endl;
    return rooplots;
  }

  RooSimultaneous* pdf = (RooSimultaneous*) work->pdf("simPdf");
  if(!pdf){
    std::cout << "ERROR::No pdf found in workspace. Will return empty vector!" << std::endl;
    return rooplots;
  }

  // Get data from workspace
  RooAbsData* data = work->data("obsData");

  // Get category components
  RooCategory* categories = work->cat("channelCat");
  // Print table with the number of data entries per category (channel)
  data->table(*((RooAbsCategory*)categories))->Print("v");

  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)work->obj("ModelConfig");
  if(!combined_config){
    std::cout << "ERROR::No model config " << " ModelConfig " << " in workspace.  Will return empty vector!" << std::endl;
    return rooplots;
  }

  const RooArgSet* obs = combined_config->GetGlobalObservables();
  if(!obs){
    std::cout << "ERROR::No observables found in ModelConfig.  Will return empty vector!" << std::endl;
    return rooplots;
  }

  RooArgList floatParsFinal = result->floatParsFinal();

  // Create NLL    
  RooAbsReal* nll = pdf->createNLL(*data, RooFit::NumCPU(2), RooFit::GlobalObservables(*obs), RooFit::Offset(true));
  
  for(Int_t i = 0; i < floatParsFinal.getSize(); i++){
    RooAbsArg* arg = floatParsFinal.at(i);
    if(!arg->InheritsFrom("RooRealVar")) continue;

    RooRealVar* par = (RooRealVar*) arg;
    TString parName = par->GetName();

    // set parameter range to readable range
    double minRange = par->getMin();
    double maxRange = par->getMax();
    if(minRange < 0.){
      par->setMin(-3.);
      par->setMax(3.);
    } 
    else{
      par->setMin(minRange);
      par->setMax(2.);
    }

    RooPlot* frame = par->frame();
    nll->plotOn(frame, RooFit::ShiftToZero());
    frame->SetMinimum(0.);
    // To be able to see the 1/2 sigma
    frame->SetMaximum(2.5);

    RooAbsReal* pll = NULL;
    if(plotPLL){ 
      pll = nll->createProfile(*par) ;
      pll->plotOn(frame, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), RooFit::NumCPU(4)); 
    }
    
    TString canName=Form("Canvas_NLL_%s_%s", name.Data(), parName.Data());
    TCanvas* c = new TCanvas(canName,canName,600,600); 
    c->cd();
    frame->Draw();

    TLegend* legend = new TLegend(0.55,0.65,0.85,0.95,"");
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);
    TLegendEntry* entry=legend->AddEntry("","NLL","l") ;
    entry->SetLineColor(kBlue);
    if(plotPLL){	
      entry=legend->AddEntry("","PLL","l") ;
      entry->SetLineColor(kRed);
      entry->SetLineStyle(kDashed);
    }
    legend->Draw();
    
    // reset parameter range to previous values
    par->setMin(minRange);
    par->setMax(maxRange);
    
    if(pll) delete pll;
    
    rooplots.push_back(c);
  }
  
  return rooplots;

}

//********************************************************************
std::vector<TCanvas*> protoana::ProtoDUNEFitUtils::PlotNuisanceParametersImpact(RooWorkspace *work, RooFitResult* fitresult, TString snapshotload, std::vector<std::string> SystsToConsider){
  //********************************************************************

  std::vector<TCanvas*> plots;

  if(!work){
    std::cerr << "ERROR::NULL workspace. Will return empty vector!" << std::endl;
    return plots;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)work->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty vector!" << std::endl;
    return plots;
  }

  // Dataset
  RooAbsData *obsdata = work->data("obsData");

  if(!fitresult){
    std::cerr << "ERROR::NULL fit result. Will return empty vector!" << std::endl;
    return plots;
  }

  protoana::ProtoDUNEFitUtils::LoadSnapshot(work, snapshotload.Data());
  protoana::ProtoDUNEFitUtils::MakeNuisanceParamsConstant(work,"");

  RooArgList floatParList = protoana::ProtoDUNEFitUtils::GetPostfitPOIList(fitresult->floatParsFinal());
  //floatParList.Print("v");

  const int nSystToConsider = SystsToConsider.size();
  const int nPOI = floatParList.getSize();

  std::vector<TH1D*> systhistovec; std::vector<TH1D*> systhistovec_prefitup; std::vector<TH1D*> systhistovec_prefitdown; std::vector<TH1D*> systhistovec_postfitup; std::vector<TH1D*> systhistovec_postfitdown; 

  for(int i=0; i < nPOI; i++){
    TH1D* systhisto             = new TH1D(Form("systhisto%i",i),            Form("systhisto%i",i),             nSystToConsider, 0, nSystToConsider);
    TH1D* systhisto_prefitup    = new TH1D(Form("systhisto%i_prefitup",i),   Form("systhisto%i_prefitup",i),    nSystToConsider, 0, nSystToConsider);
    TH1D* systhisto_prefitdown  = new TH1D(Form("systhisto%i_prefitdown",i), Form("systhisto%i_prefitdown",i),  nSystToConsider, 0, nSystToConsider);
    TH1D* systhisto_postfitup   = new TH1D(Form("systhisto%i_postfitup",i),  Form("systhisto%i_postfitup",i),   nSystToConsider, 0, nSystToConsider);
    TH1D* systhisto_postfitdown = new TH1D(Form("systhisto%i_postfitdown",i),Form("systhisto%i_postfitdown",i), nSystToConsider, 0, nSystToConsider);
    
    systhistovec.push_back(systhisto);
    systhistovec_prefitup.push_back(systhisto_prefitup);
    systhistovec_prefitdown.push_back(systhisto_prefitdown);
    systhistovec_postfitup.push_back(systhisto_postfitup);
    systhistovec_postfitdown.push_back(systhisto_postfitdown);
  }

  for(int i=0; i < nSystToConsider; i++){
    
    TString systname = TString("alpha_") + TString(SystsToConsider[i].c_str());
    RooRealVar* var = work->var(systname.Data());
    if(!var){
      std::cerr << "ERROR::Variable " << systname.Data() << " not found in workspace. Skip!" << std::endl;
      continue;
    }
    
    double originalval = var->getVal();
    double originalerror = var->getError();

    for(int j=0; j < nPOI; j++){
      systhistovec[j]->SetBinContent(i+1, originalval);
      systhistovec[j]->SetBinError(i+1, originalerror);
      systhistovec[j]->GetXaxis()->SetBinLabel(i+1, SystsToConsider[i].c_str());
      systhistovec[j]->SetTitle( (floatParList.at(j))->GetName() );
    }

    // Using the pre-fit error
    var->setVal(originalval + 1.0);
    RooFitResult *tempfitresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer("Minuit2"), RooFit::Save(), RooFit::Strategy(1), RooFit::Minos(false), RooFit::PrintLevel(1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true) );
    RooArgList paramsfit = protoana::ProtoDUNEFitUtils::GetPostfitPOIList(tempfitresult->floatParsFinal());
    //paramsfit.Print("v");

    var->setVal(originalval - 1.0);
    RooFitResult *tempfitresult2 = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer("Minuit2"), RooFit::Save(), RooFit::Strategy(1), RooFit::Minos(false), RooFit::PrintLevel(1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true) );
    RooArgList paramsfit2 = protoana::ProtoDUNEFitUtils::GetPostfitPOIList(tempfitresult2->floatParsFinal());
    
    var->setVal(originalval + originalerror);
    RooFitResult *tempfitresult3 = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer("Minuit2"), RooFit::Save(), RooFit::Strategy(1), RooFit::Minos(false), RooFit::PrintLevel(1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true) );
    RooArgList paramsfit3 = protoana::ProtoDUNEFitUtils::GetPostfitPOIList(tempfitresult3->floatParsFinal());
    
    var->setVal(originalval - originalerror);
    RooFitResult *tempfitresult4 = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer("Minuit2"), RooFit::Save(), RooFit::Strategy(1), RooFit::Minos(false), RooFit::PrintLevel(1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true) );
    RooArgList paramsfit4 = protoana::ProtoDUNEFitUtils::GetPostfitPOIList(tempfitresult4->floatParsFinal());

    TIterator* Inititer = floatParList.createIterator();
    RooRealVar* Initvar(0);
    for(int j=0; (Initvar = (RooRealVar*)Inititer->Next()); j++){
      double DMu = 0;
      TString Initvarname(Initvar->GetName());
      if(Initvarname.Contains("alpha_") || Initvarname.Contains("gamma_")) continue;
      //Double_t Initvarerror = Initvar->getError();
      // Protection
      //if(Initvarerror == 0) Initvarerror = 1.0;

      TIterator* iter = paramsfit.createIterator();
      RooRealVar* Postvar(0);
      for(int k=0; (Postvar = (RooRealVar*)iter->Next()); k++){
	TString varname(Postvar->GetName());
	if(varname == Initvarname){
	  DMu = Postvar->getVal() - Initvar->getVal();
	  break;
	}
      }
      systhistovec_postfitup[j]->SetBinContent(i+1, DMu);

      DMu = 0;
      TIterator* iter2 = paramsfit2.createIterator();
      RooRealVar* Postvar2(0);
      for(int k=0; (Postvar2 = (RooRealVar*)iter2->Next()); k++){
	TString varname(Postvar2->GetName());
	if(varname == Initvarname){
	  //DTheta = (Postvar->getVal() - Initvar->getVal())/Initvarerror;
	  DMu = Postvar2->getVal() - Initvar->getVal();
	  break;
	}
      }
      systhistovec_postfitdown[j]->SetBinContent(i+1, DMu);
      
      DMu = 0;
      TIterator* iter3 = paramsfit3.createIterator();
      RooRealVar* Postvar3(0);
      for(int k=0; (Postvar3 = (RooRealVar*)iter3->Next()); k++){
	TString varname(Postvar3->GetName());
	if(varname == Initvarname){
	  DMu = Postvar3->getVal() - Initvar->getVal();
	  break;
	}
      }
      systhistovec_prefitup[j]->SetBinContent(i+1, DMu);
      
      DMu = 0;
      TIterator* iter4 = paramsfit4.createIterator();
      RooRealVar* Postvar4(0);
      for(int k=0; (Postvar4 = (RooRealVar*)iter4->Next()); k++){
	TString varname(Postvar4->GetName());
	if(varname == Initvarname){
	  DMu = Postvar4->getVal() - Initvar->getVal();
	  break;
	}
      }
      systhistovec_prefitdown[j]->SetBinContent(i+1, DMu);
      
      delete iter;
      delete iter2;
      delete iter3;
      delete iter4;
    }

    delete Inititer;
    
    // Back to original values
    var->setVal(originalval);
    var->setError(originalerror);   
    
  }

  TLine* lp1 = new TLine(0,-1.,nSystToConsider,-1.);
  lp1->SetLineStyle(kDashed);
  TLine* lp2 = new TLine(0,1.,nSystToConsider,1.);
  lp2->SetLineStyle(kDashed);
  TLine* lp3 = new TLine(0,0,nSystToConsider,0);
  lp3->SetLineColor(kBlack);

  TLegend* legend = new TLegend(0.20,0.80,0.85,0.90,"");
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.036);
  
  legend->SetNColumns(3);
  TLegendEntry* entry = legend->AddEntry("","Pull","pl") ;
  entry->SetMarkerColor(kBlue);
  entry->SetMarkerStyle(20);
  entry=legend->AddEntry("","Pre-fit impact on #mu","f") ;
  entry->SetFillColor(kYellow);
  entry->SetLineColor(kYellow);
  entry->SetFillStyle(3001);
  entry=legend->AddEntry("","Post-fit impact on #mu","f") ;
  entry->SetFillColor(kBlue);
  entry->SetFillStyle(3005);

  TGaxis *Maxis = new TGaxis(nSystToConsider,-1.4, nSystToConsider, 1.4,-1.4,1.4,510,"+L");
  Maxis->SetLabelSize(0.03);
  Maxis->SetTextFont(72);
  Maxis->SetLabelOffset(0.025);
  Maxis->SetTitleOffset(1.1);
  Maxis->SetTitle("#Delta#mu");

  for(int i=0; i < nPOI; i++){
    systhistovec_prefitup[i]->SetLineColor(kBlue);
    systhistovec_prefitup[i]->SetFillColor(kBlue);
    systhistovec_prefitup[i]->SetFillStyle(3005);
    systhistovec_prefitdown[i]->SetLineColor(kBlue);
    systhistovec_prefitdown[i]->SetFillColor(kBlue);
    systhistovec_prefitdown[i]->SetFillStyle(3005);
    systhistovec_postfitup[i]->SetLineColor(kYellow);
    systhistovec_postfitup[i]->SetFillColor(kYellow);
    systhistovec_postfitdown[i]->SetLineColor(kYellow);
    systhistovec_postfitdown[i]->SetFillColor(kYellow);
    systhistovec[i]->SetMarkerColor(4);
    systhistovec[i]->SetMarkerStyle(20);
    systhistovec[i]->SetLineWidth(2);
    systhistovec[i]->SetStats(false);
    systhistovec[i]->GetYaxis()->SetRangeUser(-1.4,1.4);
    systhistovec[i]->GetYaxis()->SetTitle("(#theta - #theta_{0}) / #Delta#theta");
    systhistovec[i]->GetYaxis()->SetTitleOffset(1.15);

    TCanvas* can = new TCanvas(Form("can_impact%i",i),Form("can_impact%i",i));
    systhistovec[i]->Draw();
    systhistovec_postfitup[i]->Draw("same");
    systhistovec_postfitdown[i]->Draw("same");
    systhistovec_prefitup[i]->Draw("same");
    systhistovec_prefitdown[i]->Draw("same");
    systhistovec[i]->Draw("e1same");
    lp1->Draw("same");
    lp2->Draw("same");
    lp3->Draw("same");
    can->Update();
    gPad->RedrawAxis();
    Maxis->Draw("same");
    legend->Draw();

    plots.push_back(can);
  }

  return plots;

}

//********************************************************************
TCanvas* protoana::ProtoDUNEFitUtils::PlotParametersPull(TTree* tree, RooWorkspace* ws){
  //********************************************************************

  if(!tree || !ws){
    std::cerr << "ERROR::No tree or workspace found. Pull plot failed!" << std::endl;
    return NULL;
  }

  if(tree->GetEntries() < 2){
    std::cerr << "ERROR::Tree has less than two entries. No pull plots! Ignore if you run on data!" << std::endl;
    return NULL;
  }

  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. No pull plots!" << std::endl;
    return NULL;
  }
  
  std::vector<Float_t> varValsPull;
  const RooArgSet* floatParsList = combined_config->GetParametersOfInterest();
  const Int_t n = floatParsList->getSize();
  if(n == 0){
    std::cerr << "ERROR::No floating parameters found. Pull plot failed!" << std::endl;
    return NULL;
  }
  varValsPull.resize( floatParsList->getSize(), -999. );

  TCanvas* cpull = new TCanvas("PullMeanSigma","PullMeanSigma");

  TH1F* pullmeanhisto  = new TH1F("pullmeanplot",  "", n+1, 0, n+1);
  TH1F* pullsigmahisto = new TH1F("pullsigmaplot", "", n+1, 0, n+1);

  Int_t status;
  tree->SetBranchAddress("status", &status);
  
  RooRealVar* var(0);
  TIterator* Itr = floatParsList->createIterator();
  Int_t counter = 0;
  for (Int_t i=0; (var = (RooRealVar*)Itr->Next()); ++i) {
    if(var->isConstant()) continue;
    TString varName = var->GetName() + TString("pull"); 

    if(varName.Contains("Lumi") || varName.Contains("binWidth") || varName.Contains("corr")) continue;
    tree->SetBranchAddress(varName.Data(), &varValsPull[i]);

    TH1F* pullhisto = new TH1F(Form("pullhisto%i",i),Form("pullhisto%i",i),100,-5,5);
    for(Int_t j=0; j < tree->GetEntries(); j++){
      tree->GetEntry(j);
      if(status%1000 == 0)
	pullhisto->Fill(varValsPull[i]);
    }

    if(pullhisto->GetEntries() == 0){
      //delete pullhisto;
      continue;
    }

    pullhisto->Fit("gaus","Q");
    TF1* fit = pullhisto->GetFunction("gaus");
    if(!fit){
      //delete pullhisto;
      continue;
    }

    pullmeanhisto->SetBinContent(counter+1, fit->GetParameter(1));
    pullmeanhisto->SetBinError(counter+1, fit->GetParError(1));
    pullmeanhisto->GetXaxis()->SetBinLabel(counter+1, varName.Data());

    pullsigmahisto->SetBinContent(counter+1, fit->GetParameter(2));
    pullsigmahisto->SetBinError(counter+1, fit->GetParError(2));
    pullsigmahisto->GetXaxis()->SetBinLabel(counter+1, varName.Data());

    counter++;
    
    //delete pullhisto;
  }

  delete Itr;

  // Decorate histograms
  pullmeanhisto->SetLineColor(2);
  pullmeanhisto->SetMarkerStyle(20);
  pullmeanhisto->SetMarkerColor(2);

  pullsigmahisto->SetLineColor(1);
  pullsigmahisto->SetMarkerStyle(21);
  pullsigmahisto->SetMarkerColor(1);

  pullmeanhisto->GetYaxis()->SetRangeUser(-3, 3);
  pullsigmahisto->GetYaxis()->SetRangeUser(-3, 3);

  pullmeanhisto->GetXaxis()->SetLabelSize(0.015);
  pullsigmahisto->GetXaxis()->SetLabelSize(0.015);

  pullmeanhisto->SetStats(0);
  pullsigmahisto->SetStats(0);

  TLegend* legend = new TLegend(0.55,0.65,0.85,0.95,"");
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.036);
  legend->AddEntry(pullmeanhisto, "Pull Mean",  "l");
  legend->AddEntry(pullsigmahisto,"Pull Sigma", "l");

  TLine *mline = new TLine(0,0.0,n,0.0);
  mline->SetLineColor(kBlue);
  TLine *sline = new TLine(0,1.0,n,1.0);
  sline->SetLineColor(kBlue);

  cpull->cd();
  pullsigmahisto->Draw("e");
  pullmeanhisto->Draw("esame");
  mline->Draw();
  sline->Draw();
  legend->Draw();
  
  return cpull;

}

//********************************************************************
TCanvas* protoana::ProtoDUNEFitUtils::PlotNuisanceParametersPull(TTree* tree, RooWorkspace* ws){
  //********************************************************************

  if(!tree || !ws){
    std::cerr << "ERROR::No tree or workspace found. Pull plot failed!" << std::endl;
    return NULL;
  }

  if(tree->GetEntries() < 2){
    std::cerr << "ERROR::Tree has less than two entries. No pull plots! Ignore if you run on data!" << std::endl;
    return NULL;
  }

  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. No pull plots!" << std::endl;
    return NULL;
  }
  
  std::vector<Float_t> varValsPull;
  const RooArgSet* nuisanceParsList = combined_config->GetNuisanceParameters();
  if(!nuisanceParsList){
    std::cerr << "ERROR::No nuisance parameters found. Nuisance parameters plot failed!" << std::endl;
    return NULL;
  }
  const Int_t n = nuisanceParsList->getSize();
  varValsPull.resize( nuisanceParsList->getSize(), -999. );

  TCanvas* cpullnui = new TCanvas("NuisancePullMeanSigma","NuisancePullMeanSigma");

  TH1F* nuispullmeanhisto  = new TH1F("nuispullmeanplot",  "", n+1, 0, n+1);
  TH1F* nuispullsigmahisto = new TH1F("nuispullsigmaplot", "", n+1, 0, n+1);

  Int_t status;
  tree->SetBranchAddress("status", &status);
  
  RooRealVar* var(0);
  TIterator* Itr = nuisanceParsList->createIterator();
  Int_t counter = 0;
  int ngammas = 0;
  for (Int_t i=0; (var = (RooRealVar*)Itr->Next()); ++i) {
    if(var->isConstant()) continue;
    TString varName = var->GetName() + TString("pull"); 

    if(!varName.Contains("alpha") && !varName.Contains("gamma_stat")) continue;
    if(varName.Contains("Lumi") || varName.Contains("binWidth") || varName.Contains("corr")) continue;
    tree->SetBranchAddress(varName.Data(), &varValsPull[i]);

    if (varName.Contains("gamma_stat"))
      ++ngammas;

    TH1F* nuispullhisto = new TH1F(Form("nuispullhisto%i",i),Form("nuispullhisto%i",i),100,-5,5);
    for(Int_t j=0; j < tree->GetEntries(); j++){
      tree->GetEntry(j);
      if(status%1000 == 0)
	nuispullhisto->Fill(varValsPull[i]);
    }

    if(nuispullhisto->GetEntries() == 0){
      //delete pullhisto;
      continue;
    }

    nuispullhisto->Fit("gaus","Q");
    TF1* fit = nuispullhisto->GetFunction("gaus");
    if(!fit){
      //delete pullhisto;
      continue;
    }

    nuispullmeanhisto->SetBinContent(counter+1, fit->GetParameter(1));
    nuispullmeanhisto->SetBinError(counter+1, fit->GetParError(1));
    //nuispullmeanhisto->GetXaxis()->SetBinLabel(counter+1, varName.Data());

    nuispullsigmahisto->SetBinContent(counter+1, fit->GetParameter(2));
    nuispullsigmahisto->SetBinError(counter+1, fit->GetParError(2));
    //nuispullsigmahisto->GetXaxis()->SetBinLabel(counter+1, varName.Data());

    counter++;
    
    //delete pullhisto;
  }

  delete Itr;

  // Decorate histograms
  nuispullmeanhisto->SetLineColor(2);
  nuispullmeanhisto->SetMarkerStyle(20);
  nuispullmeanhisto->SetMarkerColor(2);

  nuispullsigmahisto->SetLineColor(1);
  nuispullsigmahisto->SetMarkerStyle(21);
  nuispullsigmahisto->SetMarkerColor(1);

  nuispullmeanhisto->GetYaxis()->SetRangeUser(-3, 3);
  nuispullsigmahisto->GetYaxis()->SetRangeUser(-3, 3);

  nuispullmeanhisto->SetStats(0);
  nuispullsigmahisto->SetStats(0);

  nuispullmeanhisto->GetXaxis()->SetTitle("Systematic ID");
  nuispullsigmahisto->GetXaxis()->SetTitle("Systematic ID");

  TLegend* legend = new TLegend(0.55,0.65,0.85,0.95,"");
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.036);
  legend->AddEntry(nuispullmeanhisto, "Pull Mean",  "l");
  legend->AddEntry(nuispullsigmahisto,"Pull Sigma", "l");

  TLine *mline = new TLine(0,0.0,n,0.0);
  mline->SetLineColor(kBlue);
  TLine *sline = new TLine(0,1.0,n,1.0);
  sline->SetLineColor(kBlue);
  TLine *vline = new TLine(n-ngammas, -3., n - ngammas, 3.);
  vline->SetLineColor(kGreen);

  cpullnui->cd();
  nuispullsigmahisto->Draw("e");
  nuispullmeanhisto->Draw("esame");
  mline->Draw();
  sline->Draw();
  legend->Draw();
  vline->Draw();
  
  return cpullnui;

}

//********************************************************************
TCanvas* protoana::ProtoDUNEFitUtils::PlotNuisanceParameters(TTree* tree, RooWorkspace* ws){
  //********************************************************************

  if(!tree || !ws){
    std::cerr << "ERROR::No tree or workspace found. Pull plot failed!" << std::endl;
    return NULL;
  }

  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. No nuisance parameters plots!" << std::endl;
    return NULL;
  }
 
  Int_t status;
  tree->SetBranchAddress("status", &status);

  RooRealVar* var(0);
  std::vector<Float_t> varValsNuisance;
  const RooArgSet* nuisanceParsList = combined_config->GetNuisanceParameters();
  if(!nuisanceParsList){
    std::cerr << "ERROR::No nuisance parameters found. Nuisance parameters plot failed!" << std::endl;
    return NULL;
  }
  const Int_t n = nuisanceParsList->getSize();
  varValsNuisance.resize( nuisanceParsList->getSize(), -999. );

  TH1F* nuisancehisto  = new TH1F("nuisanceplot",  "", n+1, 0, n+1);

  TIterator* Itr = nuisanceParsList->createIterator();
  for (Int_t i=0; (var = (RooRealVar*)Itr->Next()); ++i) {
    TString varName = var->GetName();
    if(!varName.Contains("alpha") && !varName.Contains("gamma_stat")) continue;
    if(varName.Contains("Lumi") || varName.Contains("nom") || varName.Contains("binWidth") || varName.Contains("err") || varName.Contains("pull")) continue;
      
    tree->SetBranchAddress(varName.Data(), &varValsNuisance[i]);

    //std::cout << "INFO::Plotting nuisance parameter " << varName << std::endl;
    TH1F* temphisto  = new TH1F("temphisto",  "temphisto", 100, -5, 5);

    for(Int_t j=0; j < tree->GetEntries(); j++){
      tree->GetEntry(j);
      if(status%1000 == 0)
	temphisto->Fill(varValsNuisance[i]);
      //nuisancehisto->SetBinContent(i+1, varValsNuisance[i]);
    }
    nuisancehisto->SetBinContent(i+1, temphisto->GetMean());
    delete temphisto;
  }

  delete Itr;

  Int_t ngammas = 0;
  // Fill the error
  TIterator* Itr2 = nuisanceParsList->createIterator();
  for (Int_t i=0; (var = (RooRealVar*)Itr2->Next()); ++i) {
    TString varName = var->GetName() + TString("err");

    if(!varName.Contains("alpha") && !varName.Contains("gamma_stat")) continue;
    if(!varName.Contains("err")) continue;
    if(varName.Contains("Lumi") || varName.Contains("nom") || varName.Contains("binWidth") || varName.Contains("pull")) continue;

    tree->SetBranchAddress(varName.Data(), &varValsNuisance[i]);

    if(varName.Contains("gamma_stat"))
      ngammas++;

    TH1F* temphisto  = new TH1F("temphisto",  "temphisto", 100, -5, 5);

    for(Int_t j=0; j < tree->GetEntries(); j++){
      tree->GetEntry(j);
      if(status%1000 == 0)
	temphisto->Fill(varValsNuisance[i]);
      //nuisancehisto->SetBinError(i+1, varValsNuisance[i]);
    }

    nuisancehisto->SetBinError(i+1, temphisto->GetMean());
    nuisancehisto->GetXaxis()->SetBinLabel(i+1, varName);
    delete temphisto;
  }

  delete Itr2;

  nuisancehisto->SetLineColor(1);
  nuisancehisto->SetMarkerStyle(21);
  nuisancehisto->SetMarkerColor(1);
  nuisancehisto->GetYaxis()->SetRangeUser(-2.0,2.0);
  nuisancehisto->GetXaxis()->SetTitle("Systematic ID");
  nuisancehisto->GetYaxis()->SetTitle("Fit Result");
  if(tree->GetEntries() > 1)
    nuisancehisto->GetYaxis()->SetTitle("Mean Fit Result From Toys");
  nuisancehisto->SetTitle("Nuisance Parameters");
  nuisancehisto->SetStats(false);

  TLine *mline = new TLine(0,0.0,n,0.0);
  mline->SetLineColor(kRed);
  TLine *sline = new TLine(0,1.0,n,1.0);
  sline->SetLineColor(kRed);
  TLine *s1line = new TLine(0,-1.0,n,-1.0);
  s1line->SetLineColor(kRed);
  TLine *vline = new TLine(n-ngammas,-2,n-ngammas,2.0);
  vline->SetLineColor(kGreen);

  TCanvas* cNuisanceParameters = new TCanvas("cNuisanceParameters","cNuisanceParamters");
  nuisancehisto->Draw("e");
  mline->Draw("same");
  sline->Draw("same");
  s1line->Draw("same");
  vline->Draw("same");

  return cNuisanceParameters;

}

//********************************************************************
TCanvas* protoana::ProtoDUNEFitUtils::PlotAverageResultsFromToys(TTree* tree, RooWorkspace* ws, TString channelname, TString catname, RooArgList * PreFit_POI){
//********************************************************************

  if(!tree || !ws){
    std::cerr << "ERROR::No tree or workspace found. Plot failed!" << std::endl;
    return NULL;
  }

  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. No plots!" << std::endl;
    return NULL;
  }

  std::vector<Float_t> varVals;
  std::vector<Float_t> varErrVals;
  const RooArgSet* floatParsList = combined_config->GetParametersOfInterest();
  if(floatParsList->getSize() == 0){
    std::cerr << "ERROR::No floating parameters found. Plot failed!" << std::endl;
    return NULL;
  }
  
  varVals.resize( floatParsList->getSize(), -999. );
  varErrVals.resize( floatParsList->getSize(), -999. );
  
  Int_t status;
  tree->SetBranchAddress("status", &status);

  Int_t NGoodToys = 0;
  for(Int_t j=0; j < tree->GetEntries(); j++){
    tree->GetEntry(j);
    if(status%1000 == 0)
      NGoodToys++;
  }

  if(NGoodToys == 0){
    std::cerr << "ERROR::Number of good fit quality is " << NGoodToys << ". Plot failed!" << std::endl;
    return NULL;
  }

  RooRealVar* var1(0);
  TIterator* Itr1 = floatParsList->createIterator();
  std::vector<TString> namevec;
  
  for(Int_t i=0; (var1 = (RooRealVar*)Itr1->Next()); ++i) {
    TString varName = TString(var1->GetName());
    if(!varName.Contains(channelname.Data())) continue;
    if(!varName.Contains(catname.Data())) continue;
    //if(varName.Contains("Lumi") || varName.Contains("binWidth") || varName.Contains("corr") || varName.Contains("Gamma")) continue;
    if(var1->isConstant()) varName += TString("_constant");
    namevec.push_back(varName);
  }
  delete Itr1;

  Int_t n = namevec.size();
  TString histoname = channelname + "_MomHisto_" + catname;
  TH1D* histo = new TH1D(histoname.Data(), histoname.Data(), n/*+1*/, 0, n/*+1*/);

  TH1D * prefit_histo = 0x0;
  if (PreFit_POI) {
    TString prefit_name = "PreFit_" + histoname;
    prefit_histo = (TH1D*)histo->Clone(prefit_name.Data());
  }

  RooRealVar* var(0);
  TIterator* Itr = floatParsList->createIterator();
  Int_t counter = n;

  for(Int_t i=0; (var = (RooRealVar*)Itr->Next()); ++i) {
    TString varName = TString(var->GetName());
    std::cout << "Ave name: " << varName << std::endl;
    if(!varName.Contains(channelname.Data())) continue;
    if(!varName.Contains(catname.Data())) continue;
    //if(varName.Contains("Lumi") || varName.Contains("binWidth") || varName.Contains("corr") || varName.Contains("Gamma")) continue;

    tree->SetBranchAddress(varName.Data(), &varVals[i]);
    TString varName2 = varName + TString("err");
    tree->SetBranchAddress(varName2.Data(), &varErrVals[i]);

    Float_t av = 0.0;
    Float_t averr = 0.0;
    for(Int_t j=0; j < tree->GetEntries(); j++){
      tree->GetEntry(j);
      if(status%1000 == 0){
	av += varVals[i]/NGoodToys;
	averr += varErrVals[i]/NGoodToys;
      }
    }

    histo->SetBinContent(counter, av);
    histo->SetBinError(counter, averr);
    
    if (PreFit_POI) {
      RooRealVar * var2(0);
      TIterator * prefit_itr = PreFit_POI->createIterator();
      while ((var2 = (RooRealVar*)prefit_itr->Next())) {
        std::cout << "Prefit:" << var2->GetName() << std::endl;
        TString var2Name(var2->GetName());
        if (var2Name == varName) {
          std::cout << "Found match " << var2->GetName() << std::endl;
          prefit_histo->SetBinContent(counter, var2->getVal());
        }
      }
    }

    counter--;
  }
  
  delete Itr;

  TLine *line = new TLine(0,1.0,n,1.0);
  line->SetLineColor(kRed);
  
  //TLine* line1 = new TLine(1,-2,1,5);
  //line1->SetLineColor(kGreen);
  //line1->SetLineStyle(kDashed);

  //TLine* line2 = new TLine(5,-2,5,5);
  //line2->SetLineColor(kGreen);
  //line2->SetLineStyle(kDashed);

  //TLine* line3 = new TLine(4,-2,4,5);
  //line3->SetLineColor(kGreen);
  //line3->SetLineStyle(kDashed);

  histoname = TString("canvas_") + histoname;
  histo->SetTitle(channelname.Data());
  //histo->GetXaxis()->SetTitle("True Momentum bin");
  histo->GetYaxis()->SetTitle("Fit result");
  if(tree->GetEntries() > 1)
    histo->GetYaxis()->SetTitle("Average Fit result");
  histo->GetYaxis()->SetRangeUser(-1.5,4.0);
  for(Int_t i=n; i > 0; i--)
    histo->GetXaxis()->SetBinLabel(i, namevec[n-i]);
  histo->SetMarkerStyle(20);
  histo->SetStats(false);
  histo->GetXaxis()->SetLabelSize(0.02);

  TCanvas* c = new TCanvas(histoname.Data(), histoname.Data());
  histo->SetFillStyle(3144);
  histo->SetFillColor(kRed);
  histo->Draw("e2");
  line->Draw("same");

  //line1->Draw("same");
  //line2->Draw("same");
  //line3->Draw("same");

  TLegend * leg = new TLegend(.65, .65, .85, .85);
  leg->AddEntry(histo, "Postfit", "lpf");

  if (PreFit_POI) {
    prefit_histo->SetMarkerColor(kBlue);
    prefit_histo->SetMarkerStyle(20);
    prefit_histo->Draw("p same");
    leg->AddEntry(prefit_histo, "Prefit", "lpf");
  }

  leg->Draw("same");

  return c;

}

//********************************************************************
RooArgList protoana::ProtoDUNEFitUtils::GetPostfitPOIList(RooArgList paramsfit, bool print){
  //********************************************************************

  if(print){
    std::cout << std::endl;
    std::cout << "INFO::Printing fractional error for each fir parameter:-" << std::endl;
  }

  RooArgList poilist;
  TIterator* iter = paramsfit.createIterator();
  RooRealVar* var(0);
  for(Int_t i=0; (var = (RooRealVar*)iter->Next()); ++i){
    TString varname(var->GetName());
    if(varname.Contains("alpha_") || varname.Contains("gamma_")) continue;
    if(varname.Contains("POI_")){
      poilist.add(*var);
      if(print)
	std::cout << varname.Data() << ": Fit result = " << var->getVal() << " +/- " << var->getError() << ". Fractional error =  " <<  var->getError()/var->getVal() << std::endl;
    }
  }

  delete iter;

  return poilist;

}

//********************************************************************
void protoana::ProtoDUNEFitUtils::MakeNuisanceParamsConstant(RooWorkspace* ws, TString exceptPar){
  //********************************************************************

  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. Exit!" << std::endl;
    return;
  }

  RooAbsPdf* pdf = combined_config->GetPdf();
  if(!pdf){
    std::cerr << "ERROR::No pdf found in ModelConfig. Exit!" << std::endl;
    return;
  }

  const RooArgSet* obs = combined_config->GetNuisanceParameters();
  if(!obs){
    std::cerr << "ERROR::No nuisance found in ModelConfig. Exit!" << std::endl;
    return;
  }

  RooArgList floatParList;
  if(obs)
    floatParList.add(*obs);

  TIterator* iter = floatParList.createIterator();
  RooAbsArg* arg;
  while((arg=(RooAbsArg*)iter->Next())){
    TString varname;
    if(arg->InheritsFrom("RooRealVar") && !arg->isConstant()){
      varname = TString(arg->GetName());
    } 
    else{
      continue; 
    }

    if(varname == exceptPar) continue;

    RooRealVar* var = ws->var(varname.Data());
    if(!var){
      std::cout << "WARNING::Can't find parameter " << varname.Data() << " in workspace." << std::endl;
      continue;
    }

    if(varname == ""){
      std::cout << "WARNING::Empty variable name. Skipping reset value." << std::endl;
      continue;
    }
    else if(varname.Contains("alpha")){
      var->setConstant(kTRUE);
    }
    else if(varname.Contains("gamma_stat")){
      var->setConstant(kTRUE);
    }
    else{
      continue;
    }
  }

  delete iter;

}

//********************************************************************
void protoana::ProtoDUNEFitUtils::ResetValues(RooWorkspace* ws, const RooArgList& parList){
  //********************************************************************

  TIterator* iter = parList.createIterator();
  RooAbsArg* arg;
  while((arg=(RooAbsArg*)iter->Next())){
    TString varname;
    if(arg->InheritsFrom("RooRealVar") && !arg->isConstant()){
      varname = TString(arg->GetName());
    } 
    else{
      continue; 
    }

    RooRealVar* var = ws->var(varname.Data());
    if(!var){
      std::cout << "WARNING::Can't find parameter " << varname.Data() << " in workspace." << std::endl;
      continue;
    }

    if(varname == ""){
      std::cout << "WARNING::Empty variable name. Skipping reset value." << std::endl;
      continue;
    }
    else if(varname.Contains("alpha")){
      var->setVal(0.0);
    }
    else if(varname.Contains("gamma_stat")){
      var->setVal(1.0);
    }
    else{
      var->setVal(1.0);
    }
  }

  delete iter;

}

//********************************************************************
void protoana::ProtoDUNEFitUtils::ResetValuesToNominal(RooWorkspace* ws, const RooArgSet& parSet){
  //********************************************************************

  TIterator* iter = parSet.createIterator();
  RooAbsArg* arg;
  while((arg=(RooAbsArg*)iter->Next())){
    TString varname;
    if(arg->InheritsFrom("RooRealVar") && !arg->isConstant()){
      varname = TString(arg->GetName());
    } 
    else{ 
      continue; 
    }

    RooRealVar* var = ws->var(varname.Data());
    if(!var){
      std::cout << "WARNING::Can't find parameter " << varname.Data() << " in workspace." << std::endl;
      continue;
    }

    if(varname == ""){
      std::cout << "WARNING::Empty variable name. Skipping reset value to nominal." << std::endl;
      continue;
    }
    else if(varname.Contains("nom_gamma_stat")){
      var->setVal(1.0);
    }
    else if(varname.Contains("nom")){
      var->setVal(0.0);
    }
    else{
      var->setVal(0.0);
    }
  }

  delete iter;

}

//********************************************************************
void protoana::ProtoDUNEFitUtils::ResetError(RooWorkspace* ws, const RooArgList& parList){
  //********************************************************************

  TIterator* iter = parList.createIterator();
  RooAbsArg* arg;
  while((arg=(RooAbsArg*)iter->Next())){
    TString varname;
    if(arg->InheritsFrom("RooRealVar") && !arg->isConstant()){
      varname = TString(arg->GetName());
    } 
    else{
      continue; 
    }

    RooRealVar* var = ws->var(varname.Data());
    if(!var){
      std::cout << "WARNING::Can't find parameter " << varname.Data() << " in workspace." << std::endl;
      continue;
    }

    if(varname == ""){
      std::cout << "WARNING::Empty variable name. Skipping reset error." << std::endl;
      continue;
    }
    else if(varname.Contains("alpha")){
      var->setError(1.0);
      if(var->getMin() < var->getVal() - 6.) var->setMin(var->getVal() - 6.);
      if(var->getMax() > var->getVal() + 6.) var->setMax(var->getVal() + 6.);
    }
    else if(varname.Contains("gamma_stat")){
      // Constraint could be either Gaus or Poisson
      RooAbsReal* constraint = (RooAbsReal*) ws->obj(Form("%s_constraint",varname.Data()));
      if(!constraint){
	std::cout << "WARNING::Constraint for variable " << varname.Data() << " not found. Skip reset error." << std::endl;
	continue;
      }

      TString constraintString = TString(constraint->IsA()->GetName());
      if(constraintString == "") continue;
      else if(constraintString.Contains("RooGaussian")){
	RooAbsReal* ErrorVar = (RooAbsReal*)ws->obj(Form("%s_sigma",varname.Data()));
	if(!ErrorVar){
	  std::cout << "WARNING::Constraint type " << constraintString.Data() << " for variable " << varname.Data() << " not found. Skip reset error." << std::endl;
	  continue;
	}
	Double_t err = ErrorVar->getVal();
	var->setError(err);
	if(var->getMin() < var->getVal() - 6.) var->setMin(var->getVal() - 6.);
	if(var->getMax() > var->getVal() + 6.) var->setMax(var->getVal() + 6.);
      }
      else if(constraintString.Contains("RooPoisson")){
	RooAbsReal* ErrorVar = (RooAbsReal*)ws->obj(Form("nom_%s",varname.Data()));
	if(!ErrorVar){
	  std::cout << "WARNING::Constraint type " << constraintString.Data() << " for variable " << varname.Data() << " not found. Skip reset error." << std::endl;
	  continue;
	}
	Double_t err = 1/sqrt(ErrorVar->getVal());
	var->setError(err);
	if(var->getMin() < var->getVal() - 6.) var->setMin(var->getVal() - 6.);
	if(var->getMax() > var->getVal() + 6.) var->setMax(var->getVal() + 6.);
      }
      else{
	std::cout << "WARNING::Unknown constraint type " << constraintString.Data() << ". Set prefit uncertainty to 0.00001." << std::endl;
      }
    }

  }

  delete iter;

}

//********************************************************************
void protoana::ProtoDUNEFitUtils::ResetAllValuesAndErrors(RooWorkspace* ws){
  //********************************************************************

  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. Will not reset values and errors." << std::endl;
    return;
  }

  RooAbsPdf* pdf = combined_config->GetPdf();
  if(!pdf){
    std::cerr << "ERROR::No pdf found in ModelConfig. Will not reset values and errors." << std::endl;
    return;
  }

  const RooArgSet* obs = combined_config->GetObservables();
  if(!obs){
    std::cerr << "ERROR::No observables found in ModelConfig. Will not reset values and errors." << std::endl;
    return;
  }

  const RooArgSet* pdfpars = pdf->getParameters(obs);
  if(!pdfpars){
    std::cerr << "ERROR::No pdf parameters found. Will not reset values and errors." << std::endl;
    return;
  }

  const RooArgSet* globalobs = combined_config->GetGlobalObservables();
  if(!globalobs){
    std::cerr << "ERROR::No global observables found in ModelConfig. Will not reset values and errors." << std::endl;
    return;
  }

  RooArgList floatParamsList;

  TIterator* iter = pdfpars->createIterator() ;
  RooAbsArg* absarg;
  while( (absarg=(RooAbsArg*)iter->Next()) ) {
    if(absarg->InheritsFrom("RooRealVar") && !absarg->isConstant()){
      floatParamsList.add(*absarg);
    }
  }
  delete iter;

  // Now reset all values and error
  protoana::ProtoDUNEFitUtils::ResetValues(ws, floatParamsList);
  protoana::ProtoDUNEFitUtils::ResetError(ws, floatParamsList);
  protoana::ProtoDUNEFitUtils::ResetValuesToNominal(ws, *globalobs);
  
}
