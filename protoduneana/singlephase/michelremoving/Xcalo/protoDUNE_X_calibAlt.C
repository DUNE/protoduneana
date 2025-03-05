#define protoDUNE_X_calibAlt_cxx
#include "protoDUNE_X_calibAlt.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <math.h> 
#include <string>
#include <TImage.h>
#include <iomanip>
#include <algorithm>
#include <TImage.h>
#include <iomanip>
#include <TSpline.h>
#include <TText.h>
#include <TFrame.h>
#include <TMinuit.h>
#include <TVectorD.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "protoduneana/singlephase/michelremoving/scripts/LanGausFit.h"
using namespace std;

////defining recombination function
double highY=375;
double lowY=225;
double highZ=460;
double lowZ=230;
std::string location="mid";
std::string tpcUsed="second";
float LAr_density=1.39;
float alp=0.93;
float bet=0.212;
//float dedx=2.08;
float dedx=1.9;
bool userecom=true;
bool sceon = true;
bool measureLifetime=false;
bool fiducialize=false;
bool lessBins=false;
bool useMPV=false;
bool useDeltaCut=false;
int numBins=144;
size_t nBins=144;

float recom_factor(float totEf){
  if (!userecom) return 1;
  float xsi=bet*dedx/(LAr_density*totEf);
  float xsi0=bet*dedx/(LAr_density*0.4867);
  float rec0=log(alp+xsi0)/xsi0;
  return (rec0*xsi)/log(alp+xsi);
}
TFile *ef = TFile::Open("$DUNE_PARDATA_DIR/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
TH3F *xneg=(TH3F*)ef->Get("Reco_ElecField_X_Neg");
TH3F *yneg=(TH3F*)ef->Get("Reco_ElecField_Y_Neg");
TH3F *zneg=(TH3F*)ef->Get("Reco_ElecField_Z_Neg");
TH3F *xpos=(TH3F*)ef->Get("Reco_ElecField_X_Pos");
TH3F *ypos=(TH3F*)ef->Get("Reco_ElecField_Y_Pos");
TH3F *zpos=(TH3F*)ef->Get("Reco_ElecField_Z_Pos");
TFile *efACA=new TFile("/exp/dune/app/users/apaudel/MC_sample_jan26/data_newfiles/Efield_correctedz.root");
TH1D *eXACA=(TH1D*)efACA->Get("med_Ef");
float tot_Ef(float xval,float yval,float zval){
  float E0value=0.4867;
  if (!sceon){
    return E0value;
  }
  else{
    if(xval>=0){
            float ex=eXACA->Interpolate(xval);
      float ey=0.0+E0value*ypos->GetBinContent(ypos->FindBin(xval,yval,zval));
      float ez=0.0+E0value*zpos->GetBinContent(zpos->FindBin(xval,yval,zval));
      return sqrt(ex*ex+ey*ey+ez*ez);
      // return ex;
    }
    else{
         float ex=eXACA->Interpolate(xval);
      float ey=0.0+E0value*yneg->GetBinContent(yneg->FindBin(xval,yval,zval));
      float ez=0.0+E0value*zneg->GetBinContent(zneg->FindBin(xval,yval,zval));
      return sqrt(ex*ex+ey*ey+ez*ez);
      // return ex;
    }
  }
}

TH3F *zpos_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Pos");
TH3F *zneg_fd=(TH3F*)ef->Get("RecoFwd_Displacement_Z_Neg");
TH3F *zpos_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Pos");
TH3F *zneg_bd=(TH3F*)ef->Get("RecoBkwd_Displacement_Z_Neg");

float zoffset(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_fd->GetBinContent(zpos_fd->FindBin(xval,yval,zval));
  }
  else{
  return zneg_fd->GetBinContent(zneg_fd->FindBin(xval,yval,zval));
  }
}
float zoffsetbd(float xval,float yval,float zval){
  if(xval>=0){
    return zpos_bd->GetBinContent(zpos_bd->FindBin(xval,yval,zval));
  }
  else{
  return zneg_bd->GetBinContent(zneg_bd->FindBin(xval,yval,zval));
  }
}

void protoDUNE_X_calibAlt::Loop(TString mn)
{

  if (fChain == 0) return;

  //int x_bin_size=5;
  //int y_bin_size = 5; // nbiny bins in y direction
  //int z_bin_size = 5; // nbinz bins in z direction
  std::cout<<"efield at the anode neg"<<tot_Ef(-352,300,300)<<std::endl;
  std::cout<<"efield at the anode pos"<<tot_Ef(352,300,300)<<std::endl;
   if(measureLifetime){fiducialize=true; useMPV=true; useDeltaCut=true; numBins=46; nBins=46;}
  ///plane_2 details
  TH1F *dqdx_X_hist_2 = new TH1F("dqdx_X_hist_2","plane_2;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
  TH1F *dedx_X_hist_2 = new TH1F("dedx_X_hist_2","plane_2;X Coordinate(cm);dE/dx(ADC/cm)",numBins,-360,360);
  TH1F *dqdx_X_correction_hist_2 = new TH1F("dqdx_X_correction_hist_2","plane_2;X Coordinate(cm);X Correction factors",numBins,-360,360);
  TH1F *dqdx_mpv_X_hist_2 = new TH1F("dqdx_mpv_X_hist_2","plane_2;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);  
  TH1F *dqdx_mpv_X_correction_hist_2 = new TH1F("dqdx_mpv_X_correction_hist_2","plane_2;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
  TH1F *corrected_dqdx_X_hist_2 = new TH1F("corrected_dqdx_X_hist_2","plane_2;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
   TH1F *corrected_mpv_dqdx_X_hist_2 = new TH1F("corrected_mpv_dqdx_X_hist_2","plane_2;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);

  TH1F *dqdx_X_hist_1 = new TH1F("dqdx_X_hist_1","plane_1;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
  TH1F *dqdx_mpv_X_hist_1 = new TH1F("dqdx_mpv_X_hist_1","plane_1;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
  TH1F *dqdx_mpv_X_correction_hist_1 = new TH1F("dqdx_mpv_X_correction_hist_1","plane_1;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
  TH1F *dedx_X_hist_1 = new TH1F("dedx_X_hist_1","plane_1;X Coordinate(cm);dE/dx(ADC/cm)",numBins,-360,360);
  TH1F *dqdx_X_correction_hist_1 = new TH1F("dqdx_X_correction_hist_1","plane_1;X Coordinate(cm);X Correction factors",numBins,-360,360);
  TH1F *corrected_dqdx_X_hist_1 = new TH1F("corrected_dqdx_X_hist_1","plane_1;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
   TH1F *corrected_mpv_dqdx_X_hist_1 = new TH1F("corrected_mpv_dqdx_X_hist_1","plane_2;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);  

  TH1F *dqdx_X_hist_0 = new TH1F("dqdx_X_hist_0","plane_0;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
  TH1F *dedx_X_hist_0 = new TH1F("dedx_X_hist_0","plane_0;X Coordinate(cm);dE/dx(ADC/cm)",numBins,-360,360);
  TH1F *dqdx_mpv_X_hist_0 = new TH1F("dqdx_mpv_X_hist_0","plane_0;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
  TH1F *dqdx_mpv_X_correction_hist_0 = new TH1F("dqdx_mpv_X_correction_hist_0","plane_0;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
  TH1F *dqdx_X_correction_hist_0 = new TH1F("dqdx_X_correction_hist_0","plane_0;X Coordinate(cm);X Correction factors",numBins,-360,360);
  TH1F *corrected_dqdx_X_hist_0 = new TH1F("corrected_dqdx_X_hist_0","plane_0;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360);
   TH1F *corrected_mpv_dqdx_X_hist_0 = new TH1F("corrected_mpv_dqdx_X_hist_0","plane_2;X Coordinate(cm);dQ/dx(ADC/cm)",numBins,-360,360); 
  //TH2F *dqdx_vs_X=new TH2F("dqdx_vs_X","dQ/dx vs X for all T0 tagged throughgoing muons;X coordinate(cm);dQ/dx(ADC/cm)",760,-380,380,1000,0,1000);
  //TH1F *max_min=new TH1F("max_min","maximum and minimum values",200,-400,400);
  //TH1F *no_hits=new TH1F("no_hits","no of hits for each bin",200,-400,400);
  TH1F *hdqdx[3];
  TH1F *hdedx[3];
  for (unsigned int i = 0; i<3; ++i){
    hdqdx[i] = new TH1F(Form("hdqdx_%d",i), Form("plane_%d;dQ/dx (ke/cm);Entries",i), 100,0,200);
    hdedx[i] = new TH1F(Form("hdedx_%d",i), Form("plane_%d;dE/dx (MeV/cm);Entries",i), 100,0,10);
  }    

  vector<vector<float>> dqdx_value_2;
  vector<vector<float>> dedx_value_2;
  vector<float> all_dqdx_value_2;
  vector<vector<float>> dqdx_frac_correction_2;
  dqdx_value_2.resize(numBins);
  dedx_value_2.resize(numBins);
  dqdx_frac_correction_2.resize(numBins);
  vector<vector<float>> dqdx_value_1;
  vector<vector<float>> dedx_value_1;
  vector<float> all_dqdx_value_1;
  vector<vector<float>> dqdx_frac_correction_1;
  dqdx_value_1.resize(numBins);
  dedx_value_1.resize(numBins);
  dqdx_frac_correction_1.resize(numBins);
  vector<vector<float>> dqdx_value_0;
  vector<vector<float>> dedx_value_0;
  vector<float> all_dqdx_value_0;
  vector<vector<float>> dqdx_frac_correction_0;
  dqdx_value_0.resize(numBins);
  dedx_value_0.resize(numBins);
  dqdx_frac_correction_0.resize(numBins);
  std::map<int, std::vector<TH1D*>> dqdx_hist;
  std::map<int, TH1D*> all_dqdx_hist;
  std::map<int, std::map<int, std::vector<TH1D*>>> dqdx_hist_fiducialize;
  std::vector<TH1D*> dqdx_mpv_X_hist_1_vec, dqdx_mpv_X_hist_0_vec, dqdx_mpv_X_hist_2_vec;
  for (size_t i=0; i<3; ++i){
  for (size_t j=0; j<3; ++j){
  dqdx_mpv_X_hist_0_vec.push_back(new TH1D(Form("dqdx_mpv_X_hist_0_%zu_%zu", i, j), Form("Plane0:%zu:%zu", i, j),numBins,-360,360)); 
  dqdx_mpv_X_hist_1_vec.push_back(new TH1D(Form("dqdx_mpv_X_hist_1_%zu_%zu", i, j), Form("Plane1:%zu:%zu", i, j),numBins,-360,360));    
  dqdx_mpv_X_hist_2_vec.push_back(new TH1D(Form("dqdx_mpv_X_hist_2_%zu_%zu", i, j), Form("Plane2:%zu:%zu", i, j),numBins,-360,360));    


  }
  }
  for (size_t i = 0; i < 3; ++i) {
    dqdx_hist[i] = std::vector<TH1D*>();
    all_dqdx_hist[i]=(new TH1D(Form("dall_qdx_%zu", i),
                                        Form("all_Plane:%zu", i),
                                        100, 0.0, 200));

  for (size_t j = 0; j < nBins; ++j) {
          dqdx_hist[i].push_back(new TH1D(Form("dqdx_%zu_%zu", i, j), 
                                        Form("Plane:%zu bin:%zu", i, j),
                                        100, 0.0, 200)); }     
    for (size_t k=0; k<9; ++k){ dqdx_hist_fiducialize[i][k]=std::vector<TH1D*>();
       for (size_t j = 0; j < nBins; ++j) {
          dqdx_hist_fiducialize[i][k].push_back(new TH1D(Form("dqdx_%zu_%zu_%zu", i,k, j),
                                        Form("Plane:%zu bin:%zu:%zu", i, k,j),
                                        100, 0.0, 200));



       }
  }}
  ////////////////////// Importing Y-Z plane fractional corrections /////////////
  fChain->GetEntry(0);
  if (run>10000) run = 0;
  TFile my_file(Form("YZcalo_mich%s_r%d.root",mn.Data(), run));
  TH2F *YZ_negativeX_hist_2= (TH2F*)my_file.Get("correction_dqdx_ZvsY_negativeX_hist_2");
  TH2F *YZ_positiveX_hist_2= (TH2F*)my_file.Get("correction_dqdx_ZvsY_positiveX_hist_2");
  ////////////////////// Importing Y-Z plane fractional corrections /////////////
  TH2F *YZ_negativeX_hist_1= (TH2F*)my_file.Get("correction_dqdx_ZvsY_negativeX_hist_1");
  TH2F *YZ_positiveX_hist_1= (TH2F*)my_file.Get("correction_dqdx_ZvsY_positiveX_hist_1");
  ////////////////////// Importing Y-Z plane fractional corrections /////////////
  TH2F *YZ_negativeX_hist_0= (TH2F*)my_file.Get("correction_dqdx_ZvsY_negativeX_hist_0");
  TH2F *YZ_positiveX_hist_0= (TH2F*)my_file.Get("correction_dqdx_ZvsY_positiveX_hist_0");
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  std::string fileName=Form("Xcalo_mich%s_r%d_alt.root",mn.Data(), run);
  if (measureLifetime==true) fileName=Form("Xcalo_mich%s_r%d_alt.root",mn.Data(), run); 
  TFile *file = new TFile(Form("%s",fileName.c_str()),"recreate");
  TTree t1("t1","a simple Tree with simple variables");//creating a tree example
  Int_t run_number;
  Double_t event_time1;
  Int_t year_month_date1;
  Int_t hour_minute_second1;
  Float_t global_med_0,global_med_1,global_med_2;
  t1.Branch("run_number",&run_number,"run_number/I");
  t1.Branch("event_time1",&event_time1,"event_time1/D");
  t1.Branch("year_month_date1",&year_month_date1,"year_month_date1/I");
  t1.Branch("hour_minute_second1",&hour_minute_second1,"hour_minute_second1/I");
  t1.Branch("global_med_0",&global_med_0,"global_med_0/F");
  t1.Branch("global_med_1",&global_med_1,"global_med_1/F");
  t1.Branch("global_med_2",&global_med_2,"global_med_2/F");
  Int_t runvalue=0;
  Double_t time1=0;
  //Filling the TTree

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t real_nentries = fChain->GetEntries();
  if (real_nentries >200000){
    cout<<"Total entries = "<<real_nentries<<endl;
    cout<<"Only use 200000 events."<<endl;
    nentries = 200000;
    real_nentries = nentries;
  }
  //Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
  // for (Long64_t jentry=0; jentry<10000;jentry++)
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%10000==0) cout<<jentry<<"/"<<real_nentries<<endl;
    if(jentry==0){
      time1=evttime;
      runvalue=run;
      year_month_date1=year_month_date;
      hour_minute_second1=hour_min_sec;
    }
    int rest_counter=0;
    int x_bin;
    int yLoc; int zLoc;
    for(int i=0; i<cross_trks; ++i){
      bool testneg=0;
      bool testpos=0;
      if(trkstartx[i]*trkendx[i]>0) continue;
      if(trkstartx[i]<-350 or trkendx[i]<-350) testneg=1;      
      if(trkstartx[i]>350 or trkendx[i]>350) testpos=1;      
      //plane 2
      if(!((TMath::Abs(trkstartx[i])>350||trkstarty[i]<50||trkstarty[i]>550||trkstartz[i]<50||trkstartz[i]>645)&&(TMath::Abs(trkendx[i])>350||trkendy[i]<50||trkendy[i]>550||trkendz[i]<50||trkendz[i]>645))) continue;
      //  if(!(((TMath::Abs(trackthetaxz[i])>1.13) && (TMath::Abs(trackthetaxz[i])<2.0))||(TMath::Abs(trackthetayz[i])>1.22 && TMath::Abs(trackthetayz[i])<1.92)))
      if(!((abs(180/TMath::Pi()*trackthetaxz[i])>60 && abs(180/TMath::Pi()*trackthetaxz[i])<120)||(abs(180/TMath::Pi()*trackthetayz[i])>80 && abs(180/TMath::Pi()*trackthetayz[i])<100))){
      for(int j=1; j<TMath::Min(ntrkhits[i][2]-1,3000); ++j){
	  if((trkhity[i][2][j]<600)&&(trkhity[i][2][j]>0)){
	    if((trkhitz[i][2][j]<695)&&(trkhitz[i][2][j]>0)){
	        yLoc=-1; zLoc=-1;
                if (trkhity[i][2][j]>25 && trkhity[i][2][j]<175) yLoc=0;
                if (trkhity[i][2][j]>225 && trkhity[i][2][j]<375) yLoc=1;
                if (trkhity[i][2][j]>425 && trkhity[i][2][j]<575) yLoc=2;

                if (trkhitz[i][2][j]>50 && trkhitz[i][2][j]<220) zLoc=0;
                if (trkhitz[i][2][j]>240 && trkhitz[i][2][j]<450) zLoc=1;
                if (trkhitz[i][2][j]>470 && trkhitz[i][2][j]<645) zLoc=2;





               if(trkhitx[i][2][j]<0 && trkhitx[i][2][j]>-360 && testneg){//negative drift
		if(trkhitx[i][2][j]<0 && trkhitx[i][2][j+1]>0) continue;
		if(trkhitx[i][2][j]<0 && trkhitx[i][2][j-1]>0) continue;
		x_bin=dqdx_X_hist_2->FindBin(trkhitx[i][2][j]);
		float YZ_correction_factor_negativeX_2=YZ_negativeX_hist_2->GetBinContent(YZ_negativeX_hist_2->FindBin(trkhitz[i][2][j],trkhity[i][2][j]));
		float recom_correction=recom_factor(tot_Ef(trkhitx[i][2][j],trkhity[i][2][j],trkhitz[i][2][j]));
		float corrected_dqdx_2=trkdqdx[i][2][j]*YZ_correction_factor_negativeX_2*recom_correction;
		if(useDeltaCut){
                if(rest_counter>0){ rest_counter=rest_counter-1; continue;}
                if (j<ntrkhits[i][2]-3){

                double deltaQ1=recom_correction*YZ_correction_factor_negativeX_2*trkdqdx[i][2][j+1];
                double deltaQ2=recom_correction*YZ_correction_factor_negativeX_2*trkdqdx[i][2][j+2];

                if (corrected_dqdx_2>80 && deltaQ1>80 && deltaQ2>80 ) {rest_counter=2; continue;}
                }
                }
                if (mn!="3") dqdx_value_2[x_bin-1].push_back(corrected_dqdx_2);
                else dqdx_value_2[x_bin-1].push_back(trkdqdx[i][2][j]);
		dedx_value_2[x_bin-1].push_back(trkdedx[i][2][j]);
                hdqdx[2]->Fill(trkdqdx[i][2][j]);
                hdedx[2]->Fill(trkdedx[i][2][j]);
                dqdx_hist[2][x_bin-1]->Fill(corrected_dqdx_2);
                if (yLoc==0 && zLoc>-1) dqdx_hist_fiducialize[2][zLoc][x_bin-1]->Fill(corrected_dqdx_2);
                else if (yLoc==1 && zLoc>-1) dqdx_hist_fiducialize[2][zLoc+3][x_bin-1]->Fill(corrected_dqdx_2);
                else if (yLoc==2 && zLoc>-1) dqdx_hist_fiducialize[2][zLoc+6][x_bin-1]->Fill(corrected_dqdx_2);

              }//X containment
	      if(trkhitx[i][2][j]>0 && trkhitx[i][2][j]<360 && testpos){//positive drift
		if(trkhitx[i][2][j]>0 && trkhitx[i][2][j+1]<0) continue;
		if(trkhitx[i][2][j]>0 && trkhitx[i][2][j-1]<0) continue;
		x_bin=dqdx_X_hist_2->FindBin(trkhitx[i][2][j]);
                  
		float YZ_correction_factor_positiveX_2=YZ_positiveX_hist_2->GetBinContent(YZ_positiveX_hist_2->FindBin(trkhitz[i][2][j],trkhity[i][2][j]));
		float recom_correction=recom_factor(tot_Ef(trkhitx[i][2][j],trkhity[i][2][j],trkhitz[i][2][j]));
		float corrected_dqdx_2=trkdqdx[i][2][j]*YZ_correction_factor_positiveX_2*recom_correction;
                if(useDeltaCut){
                if(rest_counter>0){ rest_counter=rest_counter-1; continue;}
                if (j<ntrkhits[i][2]-3){ 
                 
                double deltaQ1=recom_correction*YZ_correction_factor_positiveX_2*trkdqdx[i][2][j+1];
                double deltaQ2=recom_correction*YZ_correction_factor_positiveX_2*trkdqdx[i][2][j+2];
     
                if (corrected_dqdx_2>80 && deltaQ1>80 && deltaQ2>80 ) {rest_counter=2; continue;}
                }
                }
                //std::cout<<TrkID[i]<<" "<<trkhitx[i][2][j]<<" "<<trkhity[i][2][j]<<" "<<trkhitz[i][2][j]<<" "<<YZ_correction_factor_positiveX_2<<" "<<trkdqdx[i][2][j]<<" "<<corrected_dqdx_2<<std::endl;
		if (mn!="3") dqdx_value_2[x_bin-1].push_back(corrected_dqdx_2);
                else dqdx_value_2[x_bin-1].push_back(trkdqdx[i][2][j]);
		dedx_value_2[x_bin-1].push_back(trkdedx[i][2][j]);
                hdqdx[2]->Fill(trkdqdx[i][2][j]);
                hdedx[2]->Fill(trkdedx[i][2][j]);
                dqdx_hist[2][x_bin-1]->Fill(corrected_dqdx_2);
                if (yLoc==0 && zLoc>-1) dqdx_hist_fiducialize[2][zLoc][x_bin-1]->Fill(corrected_dqdx_2);
                else if (yLoc==1 && zLoc>-1) dqdx_hist_fiducialize[2][zLoc+3][x_bin-1]->Fill(corrected_dqdx_2);
	        else if (yLoc==2 && zLoc>-1) dqdx_hist_fiducialize[2][zLoc+6][x_bin-1]->Fill(corrected_dqdx_2);

              }//X containment
	    } // Z containment
	  } // Y containment
	} // loop over hits of the track in the given plane
      }//angular cut
    


      ////plane_1
	rest_counter=0;
      for(int j=1; j<TMath::Min(ntrkhits[i][1]-1,3000); ++j){
	if((trkhity[i][1][j]<600)&&(trkhity[i][1][j]>0)){
	  if((trkhitz[i][1][j]<695)&&(trkhitz[i][1][j]>0)){
	  
                yLoc=-1; zLoc=-1;
                if (trkhity[i][1][j]>25 && trkhity[i][1][j]<175) yLoc=0;
                if (trkhity[i][1][j]>225 && trkhity[i][1][j]<375) yLoc=1;
                if (trkhity[i][1][j]>425 && trkhity[i][1][j]<575) yLoc=2;

                if (trkhitz[i][1][j]>50 && trkhitz[i][1][j]<220) zLoc=0;
                if (trkhitz[i][1][j]>240 && trkhitz[i][1][j]<450) zLoc=1;
                if (trkhitz[i][1][j]>470 && trkhitz[i][1][j]<645) zLoc=2;







            if(trkhitx[i][1][j]<0 && trkhitx[i][1][j]>-360 && testneg){
	      if(abs(180/TMath::Pi()*trackthetaxz[i])>140){
	    	if(trkhitx[i][1][j]<0 && trkhitx[i][1][j+1]>0) continue;
		if(trkhitx[i][1][j]<0 && trkhitx[i][1][j-1]>0) continue;
		x_bin=dqdx_X_hist_1->FindBin(trkhitx[i][1][j]);
		float recom_correction=recom_factor(tot_Ef(trkhitx[i][1][j],trkhity[i][1][j],trkhitz[i][1][j]));
	     	float YZ_correction_factor_negativeX_1=YZ_negativeX_hist_1->GetBinContent(YZ_negativeX_hist_1->FindBin(trkhitz[i][1][j],trkhity[i][1][j]));
		float corrected_dqdx_1=trkdqdx[i][1][j]*YZ_correction_factor_negativeX_1*recom_correction;
                if(useDeltaCut){
                if(rest_counter>0){ rest_counter=rest_counter-1; continue;}

                if (j<ntrkhits[i][1]-3){

                double deltaQ1=recom_correction*YZ_correction_factor_negativeX_1*trkdqdx[i][1][j+1];
                double deltaQ2=recom_correction*YZ_correction_factor_negativeX_1*trkdqdx[i][1][j+2];

                if (corrected_dqdx_1>80 && deltaQ1>80 && deltaQ2>80 ) {rest_counter=2; continue;}
                }
                }
		if (mn!="3") dqdx_value_1[x_bin-1].push_back(corrected_dqdx_1);

                else dqdx_value_1[x_bin-1].push_back(trkdqdx[i][1][j]);
		dedx_value_1[x_bin-1].push_back(trkdedx[i][1][j]);
                hdqdx[1]->Fill(trkdqdx[i][1][j]);
                hdedx[1]->Fill(trkdedx[i][1][j]);
                dqdx_hist[1][x_bin-1]->Fill(corrected_dqdx_1);
                if (yLoc==0 && zLoc>-1) dqdx_hist_fiducialize[1][zLoc][x_bin-1]->Fill(corrected_dqdx_1);
                else if (yLoc==1 && zLoc>-1) dqdx_hist_fiducialize[1][zLoc+3][x_bin-1]->Fill(corrected_dqdx_1);
                else if (yLoc==2 && zLoc>-1) dqdx_hist_fiducialize[1][zLoc+6][x_bin-1]->Fill(corrected_dqdx_1);






              }
	    }
	    if(trkhitx[i][1][j]>0 && trkhitx[i][1][j]<360 && testpos){
	      if(trkhitx[i][1][j]>0 && trkhitx[i][1][j+1]<0) continue;
	      if(trkhitx[i][1][j]>0 && trkhitx[i][1][j-1]<0) continue;
	      if(abs(180/TMath::Pi()*trackthetaxz[i])<40){
		x_bin=dqdx_X_hist_1->FindBin(trkhitx[i][1][j]);
		float recom_correction=recom_factor(tot_Ef(trkhitx[i][1][j],trkhity[i][1][j],trkhitz[i][1][j]));
		float YZ_correction_factor_positiveX_1=YZ_positiveX_hist_1->GetBinContent(YZ_positiveX_hist_1->FindBin(trkhitz[i][1][j],trkhity[i][1][j]));
		float corrected_dqdx_1=trkdqdx[i][1][j]*YZ_correction_factor_positiveX_1*recom_correction;
                if(useDeltaCut){
                if(rest_counter>0){ rest_counter=rest_counter-1; continue;}

                if (j<ntrkhits[i][1]-3){

                double deltaQ1=recom_correction*YZ_correction_factor_positiveX_1*trkdqdx[i][1][j+1];
                double deltaQ2=recom_correction*YZ_correction_factor_positiveX_1*trkdqdx[i][1][j+2];

                if (corrected_dqdx_1>80 && deltaQ1>80 && deltaQ2>80 ) {rest_counter=2; continue;}
                }
                }
		if (mn!="3") dqdx_value_1[x_bin-1].push_back(corrected_dqdx_1);
                else dqdx_value_1[x_bin-1].push_back(trkdqdx[i][1][j]);
		dedx_value_1[x_bin-1].push_back(trkdedx[i][1][j]);
                hdqdx[1]->Fill(trkdqdx[i][1][j]);
                hdedx[1]->Fill(trkdedx[i][1][j]);
                dqdx_hist[1][x_bin-1]->Fill(corrected_dqdx_1);
                if (yLoc==0 && zLoc>-1) dqdx_hist_fiducialize[1][zLoc][x_bin-1]->Fill(corrected_dqdx_1);
                else if (yLoc==1 && zLoc>-1) dqdx_hist_fiducialize[1][zLoc+3][x_bin-1]->Fill(corrected_dqdx_1);
                else if (yLoc==2 && zLoc>-1) dqdx_hist_fiducialize[1][zLoc+6][x_bin-1]->Fill(corrected_dqdx_1);





              }
	    }
	  } // Z containment
	} // Y containment
	
      } // loop over hits of the track in the given plane
   	rest_counter=0;
      /////plane_0
      for(int j=1; j<TMath::Min(ntrkhits[i][0]-1,3000); ++j){
	if((trkhity[i][0][j]<600)&&(trkhity[i][0][j]>0)){
	  if((trkhitz[i][0][j]<695)&&(trkhitz[i][0][j]>0)){

                yLoc=-1; zLoc=-1;
                if (trkhity[i][0][j]>25 && trkhity[i][0][j]<175) yLoc=0;
                if (trkhity[i][0][j]>225 && trkhity[i][0][j]<375) yLoc=1;
                if (trkhity[i][0][j]>425 && trkhity[i][0][j]<575) yLoc=2;

                if (trkhitz[i][0][j]>50 && trkhitz[i][0][j]<220) zLoc=0;
                if (trkhitz[i][0][j]>240 && trkhitz[i][0][j]<450) zLoc=1;
                if (trkhitz[i][0][j]>470 && trkhitz[i][0][j]<645) zLoc=2;






	    if(trkhitx[i][0][j]<0 && trkhitx[i][0][j]>-360 && testneg){
	      if(abs(180/TMath::Pi()*trackthetaxz[i])<40){
		if(trkhitx[i][0][j]<0 && trkhitx[i][0][j+1]>0) continue;
		if(trkhitx[i][0][j]<0 && trkhitx[i][0][j-1]>0) continue;
		x_bin=dqdx_X_hist_0->FindBin(trkhitx[i][0][j]);
		float recom_correction=recom_factor(tot_Ef(trkhitx[i][0][j],trkhity[i][0][j],trkhitz[i][0][j]));
		float YZ_correction_factor_negativeX_0=YZ_negativeX_hist_0->GetBinContent(YZ_negativeX_hist_0->FindBin(trkhitz[i][0][j],trkhity[i][0][j]));
		float corrected_dqdx_0=trkdqdx[i][0][j]*YZ_correction_factor_negativeX_0*recom_correction;
                if(useDeltaCut){
                if(rest_counter>0){ rest_counter=rest_counter-1; continue;}
                if (j<ntrkhits[i][0]-3){

                double deltaQ1=recom_correction*YZ_correction_factor_negativeX_0*trkdqdx[i][0][j+1];
                double deltaQ2=recom_correction*YZ_correction_factor_negativeX_0*trkdqdx[i][0][j+2];

                if (corrected_dqdx_0>80 && deltaQ1>80 && deltaQ2>80 ) {rest_counter=2; continue;}
                }
                }

		if (mn!="3") dqdx_value_0[x_bin-1].push_back(corrected_dqdx_0);
                else  dqdx_value_0[x_bin-1].push_back(trkdqdx[i][0][j]);
		dedx_value_0[x_bin-1].push_back(trkdedx[i][0][j]);
                hdqdx[0]->Fill(trkdqdx[i][0][j]);
                hdedx[0]->Fill(trkdedx[i][0][j]);
                dqdx_hist[0][x_bin-1]->Fill(corrected_dqdx_0);
                //if (x_bin == 1) cout<<event<<" neg "<<x_bin<<" "<<trkhitx[i][0][j]<<" "<<trkhity[i][0][j]<<" "<<trkhitz[i][0][j]<<endl;
               if (yLoc==0 && zLoc>-1) dqdx_hist_fiducialize[0][zLoc][x_bin-1]->Fill(corrected_dqdx_0);
                else if (yLoc==1 && zLoc>-1) dqdx_hist_fiducialize[0][zLoc+3][x_bin-1]->Fill(corrected_dqdx_0);
                else if (yLoc==2 && zLoc>-1) dqdx_hist_fiducialize[0][zLoc+6][x_bin-1]->Fill(corrected_dqdx_0);



	     


 }
	    }
	    if(trkhitx[i][0][j]>0 && trkhitx[i][0][j]<360 && testpos){
	      if(abs(180/TMath::Pi()*trackthetaxz[i])>140){
		if(trkhitx[i][0][j]>0 && trkhitx[i][0][j+1]<0) continue;
		if(trkhitx[i][0][j]>0 && trkhitx[i][0][j-1]<0) continue;
		x_bin=dqdx_X_hist_0->FindBin(trkhitx[i][0][j]);
		float recom_correction=recom_factor(tot_Ef(trkhitx[i][0][j],trkhity[i][0][j],trkhitz[i][0][j]));
		float YZ_correction_factor_positiveX_0=YZ_positiveX_hist_0->GetBinContent(YZ_positiveX_hist_0->FindBin(trkhitz[i][0][j],trkhity[i][0][j]));
		float corrected_dqdx_0=trkdqdx[i][0][j]*YZ_correction_factor_positiveX_0*recom_correction;
                if(useDeltaCut){
                if(rest_counter>0){ rest_counter=rest_counter-1; continue;}
                if (j<ntrkhits[i][0]-3){

                double deltaQ1=recom_correction*YZ_correction_factor_positiveX_0*trkdqdx[i][0][j+1];
                double deltaQ2=recom_correction*YZ_correction_factor_positiveX_0*trkdqdx[i][0][j+2];

                if (corrected_dqdx_0>80 && deltaQ1>80 && deltaQ2>80 ) {rest_counter=2; continue;}
                }
                }

		if (mn!="3") dqdx_value_0[x_bin-1].push_back(corrected_dqdx_0);
                else  dqdx_value_0[x_bin-1].push_back(trkdqdx[i][0][j]);
		dedx_value_0[x_bin-1].push_back(trkdedx[i][0][j]);
                hdqdx[0]->Fill(trkdqdx[i][0][j]);
                hdedx[0]->Fill(trkdedx[i][0][j]);
                dqdx_hist[0][x_bin-1]->Fill(corrected_dqdx_0);
                if (yLoc==0 && zLoc>-1) dqdx_hist_fiducialize[0][zLoc][x_bin-1]->Fill(corrected_dqdx_0);
                else if (yLoc==1 && zLoc>-1) dqdx_hist_fiducialize[0][zLoc+3][x_bin-1]->Fill(corrected_dqdx_0);
                else if (yLoc==2 && zLoc>-1) dqdx_hist_fiducialize[0][zLoc+6][x_bin-1]->Fill(corrected_dqdx_0);



                //if (x_bin == 1) cout<<event<<" pos "<<x_bin<<" "<<trkhitx[i][0][j]<<" "<<trkhity[i][0][j]<<" "<<trkhitz[i][0][j]<<endl;
	      }
	    }
	  } // Z containment
	} // Y containment
      } // loop over hits of the track in the given plane

    } // loop over crossing tracks in the event
  } // loop over jentries

  std::cout << "*************** Calculating the local median dQ/dx values for each Y-Z cell ******************" << std::endl;

  ////////////////////// Getting inforamtion from dqdx_value vector ////////////////////
  for(size_t i=0; i<dqdx_value_2.size(); i++){
    //std::cout<<"no of entries in each bin plane 2 "<<i<<"  "<<dqdx_value_2[i].size()<<std::endl;
    if(dqdx_value_2[i].size()>5){
      for(size_t k=0; k<dqdx_value_2[i].size(); k++){
	all_dqdx_value_2.push_back(dqdx_value_2[i][k]);
        all_dqdx_hist[2]->Fill(dqdx_value_2[i][k]);
      }
      float local_median_dqdx_2=TMath::Median(dqdx_value_2[i].size(),&dqdx_value_2[i][0]);
      dqdx_X_hist_2->SetBinContent(i+1,local_median_dqdx_2);
      float local_median_dedx_2=TMath::Median(dedx_value_2[i].size(),&dedx_value_2[i][0]);
      dedx_X_hist_2->SetBinContent(i+1,local_median_dedx_2);
    }
  }
  if (run>10000) run = 0;
  float global_median_dqdx_2=TMath::Median(all_dqdx_value_2.size(),&all_dqdx_value_2[0]); 
  global_med_2=global_median_dqdx_2;//Filling the Tree variable
  ofstream outfile0,outfile1,outfile2;
  outfile2.open(Form("global_median_2_r%d.txt", run));
  outfile2<<run<<"\t"<<global_median_dqdx_2<<std::endl;
  //////////////////////////////////////////////////////////////////////////////////////
 
  std::cout << "**************** Calculating fractional correction for each x cell *********************" << std::endl;
 
  ///////////////////////// Calculating fractional corrections /////////////////////////////
 
  for(size_t i=0; i<dqdx_value_2.size(); i++){
    if(dqdx_value_2[i].size()>5){
      float local_median_dqdx_2=TMath::Median(dqdx_value_2[i].size(),&dqdx_value_2[i][0]);
      float fractional_dqdx_2=float(global_median_dqdx_2)/local_median_dqdx_2;
      dqdx_X_correction_hist_2->SetBinContent(i+1,fractional_dqdx_2);
      dqdx_frac_correction_2[i].push_back(fractional_dqdx_2);
    }
  }
 
  ////////////////////////////////////////////////////////////////////////////////////
  std::cout << "**************** Calculating XYZ corrected dQ/dx values ********************" << std::endl;
  //////////////// How corrected X dqdx distribution looks like /////////////////////
  for(size_t i=0; i<dqdx_value_2.size(); i++){
    if(dqdx_value_2[i].size()>5){
      float local_median_dqdx_2=TMath::Median(dqdx_value_2[i].size(),&dqdx_value_2[i][0]);
      float Xcorrected_dqdx_2=local_median_dqdx_2*dqdx_frac_correction_2[i][0];
      corrected_dqdx_X_hist_2->SetBinContent(i+1,Xcorrected_dqdx_2);
    }
  }
 
  //////////////////////////////////////////////////////////////////////////////////
 
  dqdx_X_hist_2->Write();
  dedx_X_hist_2->Write();
  dqdx_X_correction_hist_2->Write();
  corrected_dqdx_X_hist_2->Write();
 ////////////////////// Getting inforamtion from dqdx_value vector ////////////////////
  for(size_t i=0; i<dqdx_value_1.size(); i++){
    //std::cout<<"no of entries in each bin plane 1 "<<i<<"  "<<dqdx_value_1[i].size()<<std::endl;
    if(dqdx_value_1[i].size()>5){
      for(size_t k=0; k<dqdx_value_1[i].size(); k++){
	all_dqdx_value_1.push_back(dqdx_value_1[i][k]);
        all_dqdx_hist[1]->Fill(dqdx_value_1[i][k]);
      }
      float local_median_dqdx_1=TMath::Median(dqdx_value_1[i].size(),&dqdx_value_1[i][0]);
      dqdx_X_hist_1->SetBinContent(i+1,local_median_dqdx_1);
      float local_median_dedx_1=TMath::Median(dedx_value_1[i].size(),&dedx_value_1[i][0]);
      dedx_X_hist_1->SetBinContent(i+1,local_median_dedx_1);
    }
  }
 
  float global_median_dqdx_1=TMath::Median(all_dqdx_value_1.size(),&all_dqdx_value_1[0]);
 global_med_1=global_median_dqdx_1;//Filling the Tree variable
  outfile1.open(Form("global_median_1_r%d.txt", run));
  outfile1<<run<<"\t"<<global_median_dqdx_1<<std::endl; 
 
  //////////////////////////////////////////////////////////////////////////////////////
 
  std::cout << "**************** Calculating fractional correction for each x cell *********************" << std::endl;
 
  ///////////////////////// Calculating fractional corrections /////////////////////////////
 
  for(size_t i=0; i<dqdx_value_1.size(); i++){
    if(dqdx_value_1[i].size()>5){
      float local_median_dqdx_1=TMath::Median(dqdx_value_1[i].size(),&dqdx_value_1[i][0]);
      float fractional_dqdx_1=float(global_median_dqdx_1)/local_median_dqdx_1;
      dqdx_X_correction_hist_1->SetBinContent(i+1,fractional_dqdx_1);
      dqdx_frac_correction_1[i].push_back(fractional_dqdx_1);
    }
  }
 
  ////////////////////////////////////////////////////////////////////////////////////
  std::cout << "**************** Calculating XYZ corrected dQ/dx values ********************" << std::endl;
  //////////////// How corrected X dqdx distribution looks like /////////////////////
  for(size_t i=0; i<dqdx_value_1.size(); i++){
    if(dqdx_value_1[i].size()>5){
      float local_median_dqdx_1=TMath::Median(dqdx_value_1[i].size(),&dqdx_value_1[i][0]);
      float Xcorrected_dqdx_1=local_median_dqdx_1*dqdx_frac_correction_1[i][0];
      corrected_dqdx_X_hist_1->SetBinContent(i+1,Xcorrected_dqdx_1);
    }
  }
  //////////////////////////////////////////////////////////////////////////////////
  dqdx_X_hist_1->Write();
  dedx_X_hist_1->Write();
  dqdx_X_correction_hist_1->Write();
  corrected_dqdx_X_hist_1->Write();


 ////////////////////// Getting inforamtion from dqdx_value vector ////////////////////
  for(size_t i=0; i<dqdx_value_0.size(); i++){
    //std::cout<<"no of entries in each bin plane 0 "<<i<<"  "<<dqdx_value_0[i].size()<<std::endl;
    if(dqdx_value_0[i].size()>5){
      for(size_t k=0; k<dqdx_value_0[i].size(); k++){
	all_dqdx_value_0.push_back(dqdx_value_0[i][k]);
        all_dqdx_hist[0]->Fill(dqdx_value_0[i][k]);
      }
      float local_median_dqdx_0=TMath::Median(dqdx_value_0[i].size(),&dqdx_value_0[i][0]);
      dqdx_X_hist_0->SetBinContent(i+1,local_median_dqdx_0);
      float local_median_dedx_0=TMath::Median(dedx_value_0[i].size(),&dedx_value_0[i][0]);
      dedx_X_hist_0->SetBinContent(i+1,local_median_dedx_0);
    }
  }
  float global_median_dqdx_0=TMath::Median(all_dqdx_value_0.size(),&all_dqdx_value_0[0]); 
 global_med_0=global_median_dqdx_0;//Filling the Tree variable
  outfile0.open(Form("global_median_0_r%d.txt",run));
  outfile0<<run<<"\t"<<global_median_dqdx_0<<std::endl; 
  run_number=runvalue;
  event_time1=time1;
  t1.Fill();//Filling the Tree
  //////////////////////////////////////////////////////////////////////////////////////
 
  std::cout << "**************** Calculating fractional correction for each x cell *********************" << std::endl;
 
  ///////////////////////// Calculating fractional corrections /////////////////////////////
 
  for(size_t i=0; i<dqdx_value_0.size(); i++){
    if(dqdx_value_0[i].size()>5){
      float local_median_dqdx_0=TMath::Median(dqdx_value_0[i].size(),&dqdx_value_0[i][0]);
      float fractional_dqdx_0=float(global_median_dqdx_0)/local_median_dqdx_0;
      dqdx_X_correction_hist_0->SetBinContent(i+1,fractional_dqdx_0);
      dqdx_frac_correction_0[i].push_back(fractional_dqdx_0);
    }
  }
 
  ////////////////////////////////////////////////////////////////////////////////////
  std::cout << "**************** Calculating XYZ corrected dQ/dx values ********************" << std::endl;
  //////////////// How corrected X dqdx distribution looks like /////////////////////
  for(size_t i=0; i<dqdx_value_0.size(); i++){
    if(dqdx_value_0[i].size()>5){
      float local_median_dqdx_0=TMath::Median(dqdx_value_0[i].size(),&dqdx_value_0[i][0]);
      float Xcorrected_dqdx_0=local_median_dqdx_0*dqdx_frac_correction_0[i][0];
      corrected_dqdx_X_hist_0->SetBinContent(i+1,Xcorrected_dqdx_0);
    }
  }
 
  //////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////
  std::cout << "**************** Calculating XYZ corrected dQ/dx values ********************" << std::endl;
  //////////////// How corrected X dqdx distribution looks like /////////////////////
  

  for(size_t i=0; i<dqdx_value_0.size(); i++){
    if(dqdx_value_0[i].size()>5){
      float local_median_dqdx_0=TMath::Median(dqdx_value_0[i].size(),&dqdx_value_0[i][0]);
      float Xcorrected_dqdx_0=local_median_dqdx_0*dqdx_frac_correction_0[i][0];
      corrected_dqdx_X_hist_0->SetBinContent(i+1,Xcorrected_dqdx_0);
    }
  }
 
  if(useMPV){
  std::vector<double> global_mpv;




  for (int i = 0; i < 3; ++i) {
    std::vector<double> chi2_vals;
      vector<double> chi_denominator;
      vector<double> chi_numerator;
      TF1 *fitsnr_global = runlangaufit(all_dqdx_hist[i],i,200);
      global_mpv.push_back(fitsnr_global->GetParameter(1));

      for (size_t j = 0; j < nBins; j++){       
        if (dqdx_hist[i][j]->GetEntries()<301) continue;
    
   // std::cout<<dqdx_hist[i][j]->GetEntries()<<std::endl;
    TF1 *fitsnr = runlangaufit(dqdx_hist[i][j], i, 200);  
     if(i==0){ dqdx_mpv_X_hist_0->SetBinContent(j+1,fitsnr->GetParameter(1));
dqdx_mpv_X_hist_0->SetBinError(j+1,fitsnr->GetParError(1));


      dqdx_mpv_X_correction_hist_0->SetBinContent(j+1,global_mpv.at(0)/fitsnr->GetParameter(1));
      corrected_mpv_dqdx_X_hist_0->SetBinContent(j+1,dqdx_mpv_X_correction_hist_0->GetBinContent(j+1)*dqdx_mpv_X_hist_0->GetBinContent(i+1));
      }

     if(i==1){ dqdx_mpv_X_hist_1->SetBinContent(j+1,fitsnr->GetParameter(1));
dqdx_mpv_X_hist_1->SetBinError(j+1,fitsnr->GetParError(1));

      dqdx_mpv_X_correction_hist_1->SetBinContent(j+1,global_mpv.at(1)/fitsnr->GetParameter(1));
      corrected_mpv_dqdx_X_hist_1->SetBinContent(j+1,dqdx_mpv_X_correction_hist_1->GetBinContent(j+1)*dqdx_mpv_X_hist_1->GetBinContent(j+1));
      }

     if(i==2){ dqdx_mpv_X_hist_2->SetBinContent(j+1,fitsnr->GetParameter(1));
dqdx_mpv_X_hist_2->SetBinError(j+1,fitsnr->GetParError(1));

      dqdx_mpv_X_correction_hist_2->SetBinContent(j+1,global_mpv.at(2)/fitsnr->GetParameter(1));
      corrected_mpv_dqdx_X_hist_2->SetBinContent(j+1,dqdx_mpv_X_correction_hist_2->GetBinContent(j+1)*dqdx_mpv_X_hist_2->GetBinContent(j+1));
      }

        
}
      for (size_t k = 0; k < 9; k++){

      for (size_t j = 0; j < nBins; j++){
    if (dqdx_hist_fiducialize[i][k][j]->GetEntries()<301) continue;
    TF1 *fitsnr2 = runlangaufit(dqdx_hist_fiducialize[i][k][j], i, 200);
      if(i==0){ dqdx_mpv_X_hist_0_vec.at(k)->SetBinContent(j+1,fitsnr2->GetParameter(1));
dqdx_mpv_X_hist_0_vec.at(k)->SetBinError(j+1,fitsnr2->GetParError(1));}   
      
      if(i==1){ dqdx_mpv_X_hist_1_vec.at(k)->SetBinContent(j+1,fitsnr2->GetParameter(1));
dqdx_mpv_X_hist_1_vec.at(k)->SetBinError(j+1,fitsnr2->GetParError(1));}

      if(i==2){ dqdx_mpv_X_hist_2_vec.at(k)->SetBinContent(j+1,fitsnr2->GetParameter(1));
dqdx_mpv_X_hist_2_vec.at(k)->SetBinError(j+1,fitsnr2->GetParError(1));}

}
}
}
}

  //////////////////////////////////////////////////////////////////////////////////
  dqdx_mpv_X_hist_0->Write();
  dqdx_mpv_X_hist_1->Write();
  dqdx_mpv_X_hist_2->Write();  
  corrected_mpv_dqdx_X_hist_0->Write();
  corrected_mpv_dqdx_X_hist_1->Write();
  corrected_mpv_dqdx_X_hist_2->Write(); 
  dqdx_mpv_X_correction_hist_0->Write();
  dqdx_mpv_X_correction_hist_1->Write();
  dqdx_mpv_X_correction_hist_2->Write();
 
  dqdx_X_hist_0->Write();
  dedx_X_hist_0->Write();
  dqdx_X_correction_hist_0->Write();
  corrected_dqdx_X_hist_0->Write();

  for (int i = 0; i<3; ++i){
    hdqdx[i]->Write();
    hdedx[i]->Write();}
  for (int j =0; j<9; ++j){
  dqdx_mpv_X_hist_0_vec.at(j)->Write();
  dqdx_mpv_X_hist_1_vec.at(j)->Write();
  dqdx_mpv_X_hist_2_vec.at(j)->Write();
  }
  t1.Write();
  file->Close();  
  dqdx_X_hist_2->Draw();
  TFile treefile(Form("globalmedians_cathanode_r%d.root",run),"RECREATE");
  t1.Write();


  std::cout << "*************** X_Correction_make_class.C macro has ended ******************" << std::endl; 
}

int main(int argc, char *argv[]) {
  
  if (!argv[2]) {
    cout << "Error: No input file or michelremoving tree number was provided!" << endl;
    cout << "Usage: " << endl;
    cout << "make_x_correction root_file_or_list [michelremoving_tree_number]" << endl;

    return 0;
  }
  
  string infile = argv[1];
  string michelnumber = argv[2];
  string sce = argv[3];
  if(argv[4]){
  string lifetimevar=argv[4];
  if(lifetimevar=="1") measureLifetime=true;
  }

  if (!(michelnumber == "0"||michelnumber == "1"||michelnumber == "2"||michelnumber == "3")){
    cout << "Error: Michel tree number must be 0,1, or 2" << endl;
    return 0;
    }

 if (michelnumber=="0") michelnumber = "";
  cout << Form("michelremoving%s/Event", michelnumber.c_str()) << endl;

  TChain* shtree = new TChain("Event");

  if (infile.substr(infile.find_last_of(".") + 1) == "root"){
    shtree->Add(Form("%s/michelremoving%s/Event", infile.c_str(), michelnumber.c_str()));
  }

  else /*if(infile.substr(infile.find_last_of(".") + 1) == "list" || infile.substr(infile.find_last_of(".") + 1) == "txt")*/{
    std::ifstream in;
    in.open(infile.c_str());
    char line[1024];

    while(1){
      in.getline(line,1024);
      if (!in.good()) break;
      shtree->Add(Form("%s/michelremoving%s/Event", line, michelnumber.c_str()));
    }
    in.close();
    in.clear();
  }

  if (sce=="0"){
    cout<<"SCE off"<<endl;
    sceon = false;
    userecom = false;
  }
  else{
    sceon = true;
    userecom = true;
    cout<<"SCE on"<<endl;
  }
  protoDUNE_X_calibAlt t(shtree);
  
  t.Loop(michelnumber.c_str());
} // main
