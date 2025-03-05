#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TTimeStamp.h"
#include <fstream>
#include "TMinuit.h"
#include "TString.h"
#include <vector>
#include <string.h>
#include "TLatex.h"
#include "TPaveStats.h"
#include "TDatime.h"
#include "TColor.h"
#include "TProfile.h"
#include "TProfile2D.h"
void plotXCaloLifetimeAlt(std::string runNumber, std::string txtLabel, std::string location, std::string tpc, int plane=2)
{
        int yLoc=1; int zLoc=1;
        if (location=="low") yLoc=0;
        if (location=="mid") yLoc=1;
        if (location=="high") yLoc=2;
        if (tpc=="first") zLoc=0;
        if (tpc=="second") zLoc=1;
        if (tpc=="third") zLoc=2;  
	std::cout<<yLoc<<","<<zLoc<<std::endl;
        ofstream outfile;
	outfile.open(Form("mich1_lifetimeresults_%s_%s_%s_plane%d_alt.txt",txtLabel.c_str(),location.c_str(),tpc.c_str(),plane),std::ios_base::app);
	gROOT->LoadMacro("protoDUNEStyle.C");
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.35);
	gStyle->SetOptFit(111);
        gStyle->SetStatX(0.85);
	gStyle->SetStatY(0.88); 
	TCanvas c1=TCanvas();
        std::string fileName="Xcalo_mich1_r13857466.root";
	std::string set=runNumber;
	std::string s=runNumber;
	if (runNumber.find(Form("MC"))!=std::string::npos){
	std::string set="MC";
	std::string s="MC";
	fileName=Form("Xcalo_mich1_mc_sceOn_alt.root");
	}
        else if (runNumber.find(Form("NoDiff"))!=std::string::npos){
        std::string set="MC No Diff";
        std::string s="MC No Diff";
        fileName=Form("Xcalo_mich1_mc_sceOn_noDiff_alt.root");
        }





        else{  
	fileName=Form("Xcalo_mich1_r%s_alt.root",runNumber.c_str());
	}
	TFile my_file1(Form("%s",fileName.c_str()));
	TTree* t1=(TTree*)my_file1.Get("t1");
	Int_t hour_min_sec; Int_t year_month_day;
	t1->SetBranchAddress("hour_minute_second1",&hour_min_sec);
	t1->SetBranchAddress("year_month_date1",&year_month_day);
	t1->GetEntry(0);
	std::cout<<hour_min_sec<<','<<year_month_day<<std::endl;

	TH1F *x_correction_factor=(TH1F*)my_file1.Get(Form("dqdx_mpv_X_hist_%d_%d_%d",plane,yLoc,zLoc));
        int nBins=x_correction_factor->GetNbinsX();
	TProfile* rebinHist=new TProfile("rebinHist","rebinHist",nBins/2-3,0.0,2.0,"S");
	TH1F* rebinHistL=new TH1F("rebinHistL","rebinHistL",nBins/2-3,0.0,2.0);
	TH1F* rebinHistR=new TH1F("rebinHistR","rebinHistR",nBins/2-3,0.0,2.0);
	double velocity=156;
    std::vector<double> time_BL;
    std::vector<double> time_BR;
    std::vector<double> beamLeft;
    std::vector<double> beamRight;
    
for (int i=1; i<=x_correction_factor->GetNbinsX();i++){
    
    double x_value=(abs(x_correction_factor->GetBinCenter(1))-abs(x_correction_factor->GetBinCenter(i)))/velocity;
    double position=rebinHist->FindBin(x_value);
     if (position>20 || position<2) continue; 

     if (x_correction_factor->GetBinContent(i)<10 || x_correction_factor->GetBinContent(i)>200) continue;
    
    
    rebinHist->Fill(x_value,x_correction_factor->GetBinContent(i));
    
    if (x_correction_factor->GetBinCenter(i)>0) {rebinHistL->SetBinContent(position,x_correction_factor->GetBinContent(i));
    rebinHistL->SetBinError(position,x_correction_factor->GetBinError(i));
    time_BL.push_back(rebinHistL->GetBinCenter(position));
    beamLeft.push_back(rebinHistL->GetBinContent(position));    
    
    }
    if (x_correction_factor->GetBinCenter(i)<0){ rebinHistR->SetBinContent(position,x_correction_factor->GetBinContent(i));
 rebinHistR->SetBinError(position,x_correction_factor->GetBinError(i));
   time_BR.push_back(rebinHistR->GetBinCenter(position));
   beamRight.push_back(rebinHistR->GetBinContent(position));
    }
    }
    if (time_BR.size()<10 || time_BL.size()<10) return;
	rebinHistR->GetXaxis()->SetTitle("Drift Time [ms]");
	rebinHistR->GetYaxis()->SetTitle("dQ/dx [(ADC count)#timestick/cm]");
	rebinHistR->GetYaxis()->SetRangeUser(0,100);
	rebinHistR->GetXaxis()->SetRangeUser(0,2.0);
	rebinHistR->GetXaxis()->CenterTitle();
	rebinHistR->GetYaxis()->CenterTitle(); 
	rebinHistR->SetTitle(Form("Run %s BR",s.c_str()));
	if (runNumber.find(Form("MC"))!=std::string::npos)  rebinHist->SetTitle(Form("%s BR",s.c_str()));
	double lifetimeEst_BR=(time_BR.at(time_BR.size()-1)-time_BR.at(0))/(TMath::Log(beamRight.at(0)/beamRight.at(beamRight.size()-1)));
        
        double initial_BR=beamRight.at(0);
        //std::reverse(beamLeft.begin(),beamLeft.end());
        //std::reverse(time_BL.begin(),time_BL.end());
        double initial_BL=beamLeft.at(0);
        double lifetimeEst_BL=(time_BL.at(time_BL.size()-1)-time_BL.at(0))/(TMath::Log(beamLeft.at(0)/beamLeft.at(beamLeft.size()-1)));
        if(lifetimeEst_BL<0 || isnan(lifetimeEst_BL)) lifetimeEst_BL=100;
        if(lifetimeEst_BR<0 || isnan(lifetimeEst_BR)) lifetimeEst_BR=100;
        std::cout<<initial_BL<<","<<lifetimeEst_BL<<std::endl;
      
        TH1D* rebinHistRCopy=(TH1D*)rebinHistR->Clone("rebinHistRCopy");
        double chi2_dof=1000000;
        TF1 *fit= new TF1("exp","[0]*exp(x/(-[1]))",0.0,2.0);
	fit->SetParName(0,"Constant");
	fit->SetParName(1,"e^{-} Lifetime [ms]");
	int iter=0;
        while ((chi2_dof>1.05 || chi2_dof<0.95) && iter<100000){
        fit->SetParameter(0,initial_BR);
	fit->SetParameter(1,lifetimeEst_BR);
        
        rebinHistRCopy->Fit("exp","Q");
        chi2_dof=fit->GetChisquare()/fit->GetNDF();  
        //std::cout<<chi2_dof<<std::endl; 
       if (chi2_dof>1.05){
        for (int bin=0; bin<rebinHistR->GetXaxis()->GetNbins(); bin++){  if (rebinHistRCopy->GetBinContent(bin+1)>0) rebinHistRCopy->SetBinError(bin+1,0.001+rebinHistRCopy->GetBinError(bin+1));}
       
        }
	if (chi2_dof<0.95){
        for (int bin=0; bin<rebinHistR->GetXaxis()->GetNbins(); bin++){ if (rebinHistRCopy->GetBinContent(bin+1)>0)  rebinHistRCopy->SetBinError(bin+1,abs(rebinHistRCopy->GetBinError(bin+1)-0.001));}        }
        iter++;
        }
        TLatex tL;

	
	double fullLife_R=fit->GetParameter(1); double fullConst_R=fit->GetParameter(0);
	double fullSystTau_R=fit->GetParError(1); double fullConstErr_R=fit->GetParError(0);
	rebinHistR->Draw("P E");
        rebinHistRCopy->Draw("P SAMES");
	tL.SetNDC();  
	tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
	tL.Draw("SAME");
        c1.Print(Form("fitXCalo%s_BR_%s_%s_plane%d_alt.png",s.c_str(),location.c_str(),tpc.c_str(),plane));
	c1.Print(Form("fitXCalo%s_BR_%s_%s_plane%d_alt.pdf",s.c_str(),location.c_str(),tpc.c_str(),plane));
	rebinHistL->GetXaxis()->SetTitle("Drift Time [ms]");
	rebinHistL->GetYaxis()->SetTitle("dQ/dx [(ADC count)#timestick/cm]");
	rebinHistL->GetYaxis()->SetRangeUser(0,100);
	rebinHistL->GetXaxis()->SetRangeUser(0,2.0);
        rebinHistL->GetXaxis()->CenterTitle();
	rebinHistL->GetYaxis()->CenterTitle(); 
	rebinHistL->SetTitle(Form("Run %s BL",s.c_str()));
	if (runNumber.find(Form("MC"))!=std::string::npos)  rebinHist->SetTitle(Form("%s BL",s.c_str()));

	TH1D* rebinHistLCopy=(TH1D*)rebinHistL->Clone("rebinHistLCopy");
        fit->SetParameter(0,initial_BL);
	fit->SetParameter(1,lifetimeEst_BL); 

	chi2_dof=1000000;
        iter=0;
        while ((chi2_dof>1.05 || chi2_dof<0.95) && iter<100000){
        fit->SetParameter(0,initial_BL);
        fit->SetParameter(1,lifetimeEst_BL);

        rebinHistLCopy->Fit("exp","Q");
        chi2_dof=fit->GetChisquare()/fit->GetNDF();
        if (chi2_dof>1.05){
        for (int bin=0; bin<rebinHistR->GetXaxis()->GetNbins(); bin++){ if (rebinHistLCopy->GetBinContent(bin+1)>0) rebinHistLCopy->SetBinError(bin+1,rebinHistLCopy->GetBinError(bin+1)+0.001);}
        }
        if (chi2_dof<0.95){
        for (int bin=0; bin<rebinHistL->GetXaxis()->GetNbins(); bin++){ if (rebinHistLCopy->GetBinContent(bin+1)>0) rebinHistLCopy->SetBinError(bin+1,abs(rebinHistLCopy->GetBinError(bin+1)-0.001));}
        }
        iter++;
        }
        double fullLife_L=fit->GetParameter(1); double fullConst_L=fit->GetParameter(0);
	double fullSystTau_L=fit->GetParError(1); double fullConstErr_L=fit->GetParError(0);
	rebinHistL->Draw("P E");
        rebinHistLCopy->Draw("P  SAMES");
	tL.SetNDC();  
	tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
	tL.Draw("SAME");
        c1.Print(Form("fitXCalo%s_BL_%s_%s_plane%d_alt.png",s.c_str(),location.c_str(),tpc.c_str(),plane));
	c1.Print(Form("fitXCalo%s_BL_%s_%s_plane%d_alt.pdf",s.c_str(),location.c_str(),tpc.c_str(),plane));
	rebinHist->GetXaxis()->SetTitle("Drift Time [ms]");
	rebinHist->GetYaxis()->SetTitle("dQ/dx [(ADC count)#timestick/cm]");
	rebinHist->GetYaxis()->SetRangeUser(0,100);
	rebinHist->GetXaxis()->SetRangeUser(0,2.0);
        rebinHist->GetXaxis()->CenterTitle();
	rebinHist->GetYaxis()->CenterTitle(); 
	rebinHist->SetTitle(Form("Run %s",s.c_str()));
	if (runNumber.find(Form("MC"))!=std::string::npos)  rebinHist->SetTitle(Form("%s",s.c_str()));
	rebinHist->Draw("p e0");
	tL.SetNDC();  
	tL.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
	tL.Draw("SAME");
	fit->SetParameter(0,initial_BL);
	fit->SetParameter(1,lifetimeEst_BL);
        chi2_dof=1000000;
        
        fit->SetParameter(0,initial_BL);
        fit->SetParameter(1,lifetimeEst_BL);

        rebinHist->Fit("exp","WQ");
        chi2_dof=fit->GetChisquare()/fit->GetNDF();
        
	double fullLife=fit->GetParameter(1); double fullConst=fit->GetParameter(0);
	double fullSystTau=fit->GetParError(1); double fullConstErr=fit->GetParError(0);
	outfile<<runNumber<<','<<year_month_day<<','<<hour_min_sec<<","<<fullLife<<","<<fullSystTau<<","<<fullConst<<","<<fullConstErr<<','<<fullLife_R<<","<<fullSystTau_R<<","<<fullConst_R<<","<<fullConstErr_R<<','<<fullLife_L<<","<<fullSystTau_L<<","<<fullConst_L<<","<<fullConstErr_L<<std::endl;

	c1.Print(Form("fitXCalo%s_%s_%s_plane%d_alt.png",s.c_str(),location.c_str(),tpc.c_str(),plane));
	c1.Print(Form("fitXCalo%s_%s_%s_plane%d_alt.pdf",s.c_str(),location.c_str(),tpc.c_str(),plane));
}

