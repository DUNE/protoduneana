#include "LanGausFit.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFitResult.h"
#include <iostream>

using namespace std;

Double_t langaufun(Double_t *x, Double_t *par) {

  Double_t invsq2pi = 0.398942280401;// Control constants
  //Double_t mpshift = -0.22278298;
  Double_t np = 500.0;
  Double_t sc = 5.0;// convolution extends to +-sc Gaussian sigmas 
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  //mpc = par[1]- mpshift * par[0];
  mpc=par[1];
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow)/np;
 
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

 TF1 *langaufit(TH1 *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t *Status){

  Int_t i;
  Char_t FunName[100];

  sprintf(FunName,"Fitfcn_%s",his->GetName());

  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MPV","Area","GSigma");
   
  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

//  his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
  //his->SetStats(0);
  TFitResultPtr fitres = his->Fit(FunName,"RBOSQ"); // fit within specified range, use ParLimits, do not plot /////////////////// Initial code use the mode "RBO" (commented by VARUNA) ///////////

  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
 
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  Status[0] = fitres->CovMatrixStatus();
 
  return (ffit);              // return fit function
}

/////////////////////////////// Function definition /////////////////////////////////////

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

  Double_t p,x,fy,fxr,fxl;
  Double_t step;
  Double_t l,lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;

  // Search for maximum

  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l = -1.0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = langaufun(&x,params);
    if (l < lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-1);

  maxx = x;
  fy = l/2;

  // Search for right x location of fy

  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-2);

  fxr = x;

  // Search for left x location of fy

  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }

  if (i == MAXCALLS)
    return (-3);

  fxl = x;
  FWHM = fxr - fxl;
  return (0);
}

TF1 *runlangaufit(TH1 *his, int plane, double xbin_max){

  double fr[2];
  double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  fr[0]=1.;//this was originally 0.
  fr[1]=xbin_max;
  sv[1] = his->GetBinCenter(his->GetMaximumBin());
  sv[2] = 0.05*his->GetEntries();
  if (his->GetMean()>5.5){
    fr[0] = 2;
    if (plane==2){
      sv[0] = 0.7;
      sv[3] = 0.7;
    }
    else{
      sv[0] = 0.4;
      sv[3] = 1.7;
    }
  }
  else if (his->GetMean()>4){
    fr[0] = 2;
    if (plane==2){
      sv[0] = 0.3;
      sv[3] = 0.5;
    }
    else{
      sv[0] = 0.2;
      sv[3] = 0.6;
    }
  }
  else{
    sv[0] = 0.08;
    if (plane==2){//collcection plane
      sv[3] = 0.07;
    }
    else{
      sv[3] = 0.15;
    }
  }

  for(int k=0; k<4; ++k){
    pllo[k] = 0.01*sv[k];
    plhi[k] = 100*sv[k];
  if (k==0) pllo[k]=0.1;
  }
  double chisqr;
  int    ndf;
  int    status;

  TF1 *fit = langaufit(his,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status);
  //cout<<"chisqr = "<<chisqr<<" ndf = "<<ndf<<" status = "<<status<<endl;
  if (chisqr/ndf > 10 || status != 3){
    cout<<"Fit failed."<<endl; 
    cout<<"chisq = "<<chisqr<<endl;
    cout<<"ndf = "<<ndf<<endl;
    cout<<"status = "<<status<<endl;
    cout<<"hist mean = "<<his->GetMean()<<endl;
    cout<<"hist rms = "<<his->GetRMS()<<endl;
    cout<<"hist entries = "<<his->GetEntries()<<endl;
    for (int i = 0; i<4; ++i){
      cout<<"Par["<<i<<"] = "<<fit->GetParameter(i)<<"+-"<<fit->GetParError(i)<<endl;
    }

//    if(his->GetMean()>4){
//      if (plane != 2){
//        sv[0] = 0.2;
//        sv[3] = 0.6;
//      }
//      else{
//        sv[0] = 0.3;
//        sv[3] = 0.5;
//      }
//    }
//
//
    cout<<"Refitting"<<endl;
    for(int k=0; k<4; ++k){
      if (k==0) continue;
      sv[k]*=1.01;
      pllo[k] = 0.01*sv[k];
      plhi[k] = 100*sv[k];
    }
    fit = langaufit(his,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status);
    cout<<"chisq = "<<chisqr<<endl;
    cout<<"ndf = "<<ndf<<endl;
    cout<<"status = "<<status<<endl;
    for (int i = 0; i<4; ++i){
      cout<<"Par["<<i<<"] = "<<fit->GetParameter(i)<<"+-"<<fit->GetParError(i)<<endl;
    }
    if (chisqr/ndf < 10 && status == 3) cout<<"Refitting was successful"<<endl;
  }

  return fit;
}
