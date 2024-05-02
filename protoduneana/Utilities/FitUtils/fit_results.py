import ROOT as RT
import numpy as np
import math
from array import array
import sys
from argparse import ArgumentParser as ap

def get_new_xsec(f, name):
  
  if name == 'abs': title = 'Abs'
  elif name == 'cex': title = 'Cex'
  elif name == 'other': title = 'OtherInel'

  hist = f.Get(f"PostFitXSec/PostFit{title}XSec")
  xs = []
  ys = []
  eyhs = []
  eyls = []
  for i in range(1, hist.GetNbinsX()+1):
    xs.append(hist.GetBinCenter(i))
    ys.append(hist.GetBinContent(i))
    eyhs.append(0.)
    eyls.append(0.)
  
  n = len(xs)
  print(n)
  new_gr = RT.TGraphAsymmErrors(
    n,
    np.array(xs),
    np.array(ys),
    np.array([0.]*n),
    np.array([0.]*n),
    np.array(eyls),
    np.array(eyhs)
  )
  return new_gr

def combine(hists, title):
  combined = RT.TH1D(title, "", 3, 0, 3)
  data1 = hists[-3]
  data2 = hists[-2]
  data3 = hists[-1]
  combined.SetBinContent(1, data1.GetBinContent(1))
  combined.SetBinContent(2, data2.GetBinContent(1))
  combined.SetBinContent(3, data3.GetBinContent(1))
  combined.GetXaxis().SetBinLabel(1, "Failed Beam Cuts")
  combined.GetXaxis().SetBinLabel(2, "No Beam Track")
  combined.GetXaxis().SetBinLabel(3, "Michel Vertex Cut")
  combined.GetXaxis().SetLabelSize(.05)
  combined.SetMinimum(0.)
  return combined
  
def rebin(hists, title):
  results = []
  for i in range(0, 4):
    h = hists[i]
    bins = []
    for j in range(1, h.GetNbinsX()+1):
      label = h.GetXaxis().GetBinLabel(j)
      if j == 1:
        bins.append(float(label.split(' ')[0]))
      bins.append(float(label.split(' ')[-1]))
    rebinned = RT.TH1D(hists[i].GetName() + title + "_rebinned", hists[i].GetTitle(),
                       len(bins) -1, array('d', bins)) 
    for j in range(1, h.GetNbinsX()+1):
      rebinned.SetBinContent(j, h.GetBinContent(j))
    results.append(rebinned)
  return results 




def fix_params_plot(args):
  #RT.gROOT.LoadMacro("~/protoDUNEStyle.C")
  #RT.gROOT.SetBatch();
  #RT.gStyle.SetOptStat(00000)
  #RT.gStyle.SetTitleAlign(33)
  #RT.gStyle.SetTitleX(.915)
  #RT.gStyle.SetTitleY(.95)
  #RT.gROOT.ForceStyle()
  RT.gROOT.LoadMacro("~/style.C")
  RT.gROOT.ProcessLine('style()')
  RT.gStyle.SetErrorX(0.5)

  #if not args.noplotstyle:
  #  import dunestyle.root as dunestyle

  f = RT.TFile.Open(args.i)

  h = f.Get('postFitPars')
  hpre = f.Get('preFitPars')
  hnew = RT.TH1D('hnew', '', args.npars, 0, args.npars)
  hnew_pre = RT.TH1D('hnew_pre', '', args.npars, 0, args.npars)
  for i in range(1, args.npars+1):
    hnew.SetBinContent(i, h.GetBinContent(i))
    hnew.SetBinError(i, h.GetBinError(i))
    print('error:', h.GetBinError(i))

    hnew_pre.SetBinContent(i, hpre.GetBinContent(i))
    hnew_pre.SetBinError(i, hpre.GetBinError(i))

  minval = 0.
  maxval = 2.
  for i in range(1, args.npars+1):
    a = hnew.GetBinContent(i) + 2.*hnew.GetBinError(i)
    if a > maxval: maxval = a

    b = hnew.GetBinContent(i) - 2.*hnew.GetBinError(i)
    if b < minval: minval = b

  print('Setting', minval, maxval)
  hnew.SetMinimum(minval)
  hnew.SetMaximum(maxval)
  hnew.GetXaxis().SetTitleOffset(.8);
  hnew.GetYaxis().SetTitleOffset(.8);
  hnew.SetTitleSize(.05, "XY");
  hnew.SetMarkerColor(RT.kRed);
  hnew.SetMarkerStyle(20);
  hnew.GetXaxis().SetRangeUser(0, args.npars)
  hnew.SetFillStyle(3001)
  hnew.SetFillColor(RT.kRed)
  hnew.SetLineWidth(0)
  hnew.SetTitle(';Parameter;Parameter Value')

  c = RT.TCanvas('cParameters')
  c.SetTicks()
  hnew.Draw('pe2')
  hnew_pre.SetLineColor(RT.kBlue)
  hnew_pre.SetMarkerColor(RT.kBlue)
  hnew_pre.SetMarkerStyle(20)
  hnew_pre.Draw('p same')
  #fOut.cd()
  l = RT.TLine(0., 1., args.npars, 1.)
  #l.SetLineWidth(2)
  l.Draw()
  hnew.Draw('pe2 same')
  leg = RT.TLegend()
  leg.AddEntry(hnew_pre, 'Pre-Fit', 'lp')
  leg.AddEntry(hnew, 'Post-Fit #pm1#sigma', 'lpf')
  leg.SetFillStyle(0)
  leg.SetLineWidth(0)
  leg.Draw()
  c.RedrawAxis()

  c2D = f.Get('cCorr') 
  c2D.SetTicks()
  c2D.Draw()

  prims = c2D.GetListOfPrimitives()
  for i in range(prims.GetSize()):
    if 'hist' in prims.At(i).GetName().lower(): break
  prims.At(i).GetXaxis().SetRangeUser(0, args.npars)
  prims.At(i).GetYaxis().SetRangeUser(0, args.npars)
  c2D.Update()
  c2D.RedrawAxis()



  if args.save:
    c.SaveAs('pars.pdf')
    c2D.SaveAs('corr.pdf')
  #else:
  #  fOut = RT.TFile(args.o, 'recreate')
  #  c.Write()
  #  fOut.Close()

  f.Close()

def normal(args):
  #RT.gROOT.LoadMacro("~/protoDUNEStyle.C")
  #RT.gROOT.LoadMacro("~/style.C")
  #RT.gROOT.ProcessLine('style()')
  
  RT.gROOT.SetBatch();
  RT.gStyle.SetOptStat(00000)
  RT.gStyle.SetErrorX(1.e-4)
  RT.gStyle.SetTitleAlign(33)
  RT.gStyle.SetTitleX(.915)
  RT.gStyle.SetTitleY(.95)
  RT.gROOT.ForceStyle()
  #fIn = RT.TFile(sys.argv[1], "OPEN");
  
  #if not args.noplotstyle:
  #  import dunestyle.root as dunestyle

  prefit_names = [
    "NominalAbsTotal", "NominalCexTotal",
    "NominalRejectedIntTotal", "NominalAPA2Total",
    "NominalFailedBeamCutsTotal", "NominalNoBeamTrackTotal",
    "NominalMichelCutTotal"
  ]
  
  prefit_hists = []
  
  total_prefit = 0.
  if args.pretune is not None:
    fPretune = RT.TFile.Open(args.pretune)
    for n in prefit_names:
      prefit_hists.append(fPretune.Get(n))
      prefit_hists[-1].SetDirectory(0)
      total_prefit += prefit_hists[-1].Integral()
    fPretune.Close()
  
  fIn = RT.TFile.Open(args.i);
  
  
  nominal_canvas = fIn.Get("cNominalAbs")
  nominal_canvas.Draw()
  n_prim = RT.gPad.GetListOfPrimitives().GetSize()
  for i in range(0, n_prim):
    if RT.gPad.GetListOfPrimitives().At(i).GetName() == "TPave":
      nominal_leg = RT.gPad.GetListOfPrimitives().At(i)
      break
  
  nominal_label = nominal_leg.GetListOfPrimitives().At(nominal_leg.GetListOfPrimitives().GetSize()-1).GetLabel()
  print(nominal_label)
  
  postfit_canvas = fIn.Get("cPostFitAbs")
  postfit_canvas.Draw()
  n_prim = RT.gPad.GetListOfPrimitives().GetSize()
  for i in range(0, n_prim):
    if RT.gPad.GetListOfPrimitives().At(i).GetName() == "TPave":
      postfit_leg = RT.gPad.GetListOfPrimitives().At(i)
      break
  
  postfit_label = postfit_leg.GetListOfPrimitives().At(postfit_leg.GetListOfPrimitives().GetSize()-1).GetLabel()
  print(postfit_label)
  
  '''
  nominal_xsec_canvas = fIn.Get("Throws/cXSecThrowAbsUnderflow")
  nominal_xsec_canvas.Draw()
  n_prim = RT.gPad.GetListOfPrimitives().GetSize()
  for i in range(0, n_prim):
    if RT.gPad.GetListOfPrimitives().At(i).GetName() == "TPave":
      nominal_xsec_leg = RT.gPad.GetListOfPrimitives().At(i)
      break
  
  nominal_xsec_label = nominal_xsec_leg.GetListOfPrimitives().At(nominal_xsec_leg.GetListOfPrimitives().GetSize()-2).GetLabel()
  print(nominal_xsec_label)
  fake_xsec_label = nominal_xsec_leg.GetListOfPrimitives().At(nominal_xsec_leg.GetListOfPrimitives().GetSize()-1).GetLabel()
  print(fake_xsec_label)
  '''
  
  data_names = [
    "Data/Data_selected_Abs_hist", "Data/Data_selected_Cex_hist",
    "Data/Data_selected_RejectedInt_hist", "Data/Data_selected_APA2_hist",
    "Data/Data_selected_FailedBeamCuts_hist", "Data/Data_selected_NoBeamTrack_hist",
    "Data/Data_selected_MichelCut_hist"
  ]
  
  data_hists = []
  total_data = 0.
  for n in data_names:
    data_hists.append(fIn.Get(n))
    total_data += data_hists[-1].Integral()
  
  if args.pretune is None:
    for n in prefit_names:
      prefit_hists.append(fIn.Get(n))
      total_prefit += prefit_hists[-1].Integral()
  
  postfit_names = [
    "PostFitAbsTotal", "PostFitCexTotal",
    "PostFitRejectedIntTotal", "PostFitAPA2Total",
    "PostFitFailedBeamCutsTotal", "PostFitNoBeamTrackTotal",
    "PostFitMichelCutTotal"
  ]
  postfit_hists = []
  for n in postfit_names:
    postfit_hists.append(fIn.Get(n))
  
  new_data_hists = []
  for h, hMC in zip(data_hists, prefit_hists):
    h.Scale(total_prefit/total_data)
    new_h = hMC.Clone()
    new_h.Reset()
    for i in range(1, h.GetNbinsX()+1):
      new_h.SetBinContent(i, h.GetBinContent(i))
    new_data_hists.append(new_h)
  
  names = [
    "Absorption", "Charge Exchange",
    "Other Interactions", "Past Fiducial Volume",
    "Failed Beam Cuts", "No Beam Track",
    "Michel Cut",
    ""
    ,"Absorption", "Charge Exchange",
    "Other Interactions", "Past Fiducial Volume",
  ]
  
  short_names = [
    "Abs", "Cex",
    "RejectedInt", "APA2",
    "FailedBeamCuts", "NoBeamTrack"
  ]
  
  short_names = [
    "Abs", "Cex",
    "RejectedInt", "APA2",
    "FailedBeamCuts", "NoBeamTrack",
    "MichelCut", "BadEvents"
    ,"AbsRebin", "CexRebin",
    "RejectedIntRebin", "APA2Rebin",
  ]
  
  XTitles = [
    "Reconstructed KE [MeV]",
    "Reconstructed KE [MeV]",
    "Reconstructed KE [MeV]",
    "", #"Reconstructed End Z [cm]",
    "",
    "",
    "", ""
    ,"Reconstructed KE [MeV]",
    "Reconstructed KE [MeV]",
    "Reconstructed KE [MeV]",
    "", #"Reconstructed End Z [cm]",
  ]
  
  #fOut = RT.TFile(sys.argv[2], "RECREATE");
  fOut = RT.TFile(args.o, "RECREATE");
  
  data_combined = combine(new_data_hists, "data_comb") 
  #data_combined = RT.TH1D("data_comb", "", 2, 0, 2)
  #data1 = new_data_hists[-2]
  #data2 = new_data_hists[-1]
  #data_combined.SetBinContent(1, data1.GetBinContent(1))
  #data_combined.SetBinContent(2, data2.GetBinContent(1))
  
  prefit_combined = combine(prefit_hists, "prefit_comb") 
  #prefit_combined = RT.TH1D("prefit_comb", "", 2, 0, 2)
  #prefit1 = prefit_hists[-2]
  #prefit2 = prefit_hists[-1]
  #prefit_combined.SetBinContent(1, prefit1.GetBinContent(1))
  #prefit_combined.SetBinContent(2, prefit2.GetBinContent(1))
  
  postfit_combined = combine(postfit_hists, "postfit_comb") 
  #postfit_combined = RT.TH1D("postfit_comb", "", 2, 0, 2)
  #postfit1 = postfit_hists[-2]
  #postfit2 = postfit_hists[-1]
  #postfit_combined.SetBinContent(1, postfit1.GetBinContent(1))
  #postfit_combined.SetBinContent(2, postfit2.GetBinContent(1))
  
  new_data_hists.append(data_combined)
  prefit_hists.append(prefit_combined)
  postfit_hists.append(postfit_combined)
  
  new_data_hists += rebin(new_data_hists, 'data')
  prefit_hists += rebin(prefit_hists, 'prefit')
  postfit_hists += rebin(postfit_hists, 'postfit')
  tt = RT.TLatex();
  tt.SetNDC();
  for i in range(0, len(prefit_hists)):
    data = new_data_hists[i]
    data.SetMarkerStyle(20);
    data.SetMarkerColor(RT.kBlack);
    #data.SetMarkerSize(.5)
    data.SetLineColor(RT.kBlack);
    data.SetTitle(names[i] + ";" + XTitles[i] + ";Events per bin [x10^{3}]")
    data.GetYaxis().SetMaxDigits(3)
    #data.SetTitleSize(.05, "XY")
    data.GetXaxis().CenterTitle()
    data.GetYaxis().CenterTitle()
    data.Scale(1.e-3)
    #data.GetXaxis().SetTitleOffset(.8)
    data.GetYaxis().SetTitleOffset(.95)
    data.SetMinimum(0.)
    #if (i == len(prefit_hists)-1):
    #  data.SetTitle(names[i])
    #else:
    #  data.SetTitle(names[i] + ";" + XTitles[i])
  
    prefit = prefit_hists[i]
    prefit.SetFillColor(0);
    prefit.SetLineStyle(9);
    prefit.SetLineColor(RT.kBlue);
    prefit.SetTitle(names[i] + ";" + XTitles[i] + ";Events per bin [x10^{3}]")
    prefit.GetXaxis().CenterTitle()
    prefit.GetYaxis().CenterTitle()
    prefit.Scale(1.e-3)
    #prefit.SetTitleSize(.05, "XY")
    #prefit.GetXaxis().SetTitleOffset(.8)
    prefit.GetYaxis().SetTitleOffset(.95)
    prefit.GetYaxis().SetMaxDigits(3)
    prefit.SetMinimum(0.)
    #if (i == len(prefit_hists)-1):
    #  prefit.SetTitle(names[i])
    #else:
    #  prefit.SetTitle(names[i] + ";Reconstructed KE (MeV)")
  
    postfit = postfit_hists[i]
    postfit.SetFillColor(0);
    postfit.SetLineColor(RT.kRed)
    postfit.SetTitle(names[i] + ";" + XTitles[i] + ";Events per bin [x10^{3}]")
    postfit.GetYaxis().SetMaxDigits(3)
    postfit.GetXaxis().CenterTitle()
    postfit.GetYaxis().CenterTitle()
    #postfit.SetTitleSize(.05, "XY")
    #postfit.GetXaxis().SetTitleOffset(.8)
    postfit.GetYaxis().SetTitleOffset(.95)
    postfit.SetMinimum(0.)
    postfit.Scale(1.e-3)
    #if (i == len(prefit_hists)-1):
    #  postfit.SetTitle(names[i])
    #else:
    #  postfit.SetTitle(names[i] + ";Reconstructed KE (MeV)")
  
    if i == 3:
      data.GetXaxis().SetLabelOffset(.01)
      prefit.GetXaxis().SetLabelOffset(.01)
      postfit.GetXaxis().SetLabelOffset(.01)
      data.GetXaxis().SetLabelSize(.07)
      prefit.GetXaxis().SetLabelSize(.07)
      postfit.GetXaxis().SetLabelSize(.07)
  
  
    c = RT.TCanvas("c" + short_names[i], "c" + short_names[i])
    c.SetTicks()
  
    post_bins = [postfit.GetBinContent(i) for i in range(1, postfit.GetNbinsX()+1)]
    pre_bins = [prefit.GetBinContent(i) for i in range(1, prefit.GetNbinsX()+1)]
    data_bins = [data.GetBinContent(i) + data.GetBinError(i) for i in range(1, data.GetNbinsX()+1)]
    if 'Michel' in short_names[i] or 'APA2' in short_names[i]:
      postfit.GetXaxis().SetBinLabel(1, '')
      prefit.GetXaxis().SetBinLabel(1, '')
      data.GetXaxis().SetBinLabel(1, '')
      
    max_post = max(post_bins)
    max_pre = max(pre_bins)
    max_data = max(data_bins)
  
    if (max_post > max_pre and max_post > max_data):
      postfit.Draw("hist")
      prefit.Draw("hist same")
      data.Draw("same e")
    elif (max_pre > max_post and max_pre > max_data):
      prefit.Draw("hist")
      postfit.Draw("hist same")
      prefit.Draw("hist same")
      data.Draw("same e")
    elif (max_data > max_pre and max_data > max_post):
      data.Draw("e")
      postfit.Draw("hist same")
      prefit.Draw("hist same")
      data.Draw("same e")
  
    tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  
  
    if (i == 0):
      l = RT.TLegend()
      l.AddEntry(postfit, "Post-Fit", "l")
      l.AddEntry(prefit, "Nominal MC", "l")
      l.AddEntry(data, "Data", "le")
      h = RT.TH1D("dummy", "", 1, 0, 1)
      h.SetLineColor(0)
      l.Write("leg")
      l.SetFillStyle(0)
      l.SetLineWidth(0)
      #l.AddEntry(h, "Nominal " + nominal_label, "")
      #l.AddEntry(h, "Post-Fit " + postfit_label, "")
    if i in [0, 8]:
      l.Draw()
    RT.gPad.RedrawAxis()
    c.Write()
    #if args.save:
    #  c.SaveAs("c" + short_names[i] + ".pdf")
  
  
  if not args.nothrows:
    abs_hist = fIn.Get("PostFitXSec/PostFitAbsHist")
    abs_xs = []
    for i in range(1, abs_hist.GetNbinsX()+1):
      abs_xs.append(abs_hist.GetBinCenter(i))
    
    if not args.nothrows:
      gr_abs = fIn.Get("Throws/grXSecThrowAbsUnderflow");
      n_abs = gr_abs.GetN()
      print(n_abs)
      abs_ys = []
      abs_eyhs = []
      abs_eyls = []
      for i in range(0, n_abs):
        abs_ys.append(gr_abs.GetY()[i])
        abs_eyhs.append(gr_abs.GetEYhigh()[i])
        abs_eyls.append(gr_abs.GetEYlow()[i])
      new_gr_abs = RT.TGraphAsymmErrors(n_abs, array('d', abs_xs), array('d', abs_ys), array('d', [0.]*n_abs), array('d', [0.]*n_abs), array('d', abs_eyls), array('d', abs_eyhs))
    
    
    cex_hist = fIn.Get("PostFitXSec/PostFitCexHist")
    cex_xs = []
    for i in range(1, cex_hist.GetNbinsX()+1):
      cex_xs.append(cex_hist.GetBinCenter(i))
    
    gr_cex = fIn.Get("Throws/grXSecThrowCexUnderflow");
    n_cex = gr_cex.GetN()
    cex_ys = []
    cex_eyhs = []
    cex_eyls = []
    for i in range(0, n_cex):
      cex_ys.append(gr_cex.GetY()[i])
      cex_eyhs.append(gr_cex.GetEYhigh()[i])
      cex_eyls.append(gr_cex.GetEYlow()[i])
    new_gr_cex = RT.TGraphAsymmErrors(n_cex, array('d', cex_xs), array('d', cex_ys), array('d', [0.]*n_cex), array('d', [0.]*n_cex), array('d', cex_eyls), array('d', cex_eyhs))
    
    
    other_hist = fIn.Get("PostFitXSec/PostFitOtherInelHist")
    if other_hist:
      other_xs = []
      for i in range(1, other_hist.GetNbinsX()+1):
        other_xs.append(other_hist.GetBinCenter(i))
      
      if not args.nothrows:
        gr_other = fIn.Get("Throws/grXSecThrowOtherInelUnderflow");
        n_other = gr_other.GetN()
        other_ys = []
        other_eyhs = []
        other_eyls = []
        for i in range(0, n_other):
          other_ys.append(gr_other.GetY()[i])
          other_eyhs.append(gr_other.GetEYhigh()[i])
          other_eyls.append(gr_other.GetEYlow()[i])
        new_gr_other = RT.TGraphAsymmErrors(n_other, array('d', other_xs), array('d', other_ys), array('d', [0.]*n_other), array('d', [0.]*n_other), array('d', other_eyls), array('d', other_eyhs))
        new_gr_other.Write('other_xsec')
    
    new_gr_abs.Write('abs_xsec')
    new_gr_cex.Write('cex_xsec')
    
    corr = fIn.Get('Throws/xsec_corr')
    cov = fIn.Get('Throws/xsec_cov')
    
    if len(args.add_files) > 0:
    #if len(sys.argv) > 3:
      add_covs = (
        [args.add_covs[0]]*len(args.add_files) if len(args.add_covs) == 1
        else args.add_covs
      )
      for af, ac in zip(args.add_files, add_covs):
        #extra_cov_file = RT.TFile(sys.argv[3], 'open')
        extra_cov_file = RT.TFile(af, 'open')
        #extra_cov = extra_cov_file.Get(sys.argv[4])
        extra_cov = extra_cov_file.Get(ac)
    
        for i in range(1, extra_cov.GetNbinsX()+1):
          for j in range(1, extra_cov.GetNbinsX()+1):
            cov.SetBinContent(i, j, cov.GetBinContent(i, j) + extra_cov.GetBinContent(i, j))
    
      #for i in range(1, extra_cov.GetNbinsX()+1):
      #  for j in range(1, extra_cov.GetNbinsX()+1):
      for i in range(1, cov.GetNbinsX()+1):
        for j in range(1, cov.GetNbinsX()+1):
          corr.SetBinContent(i, j, cov.GetBinContent(i, j)/math.sqrt(cov.GetBinContent(i, i)*cov.GetBinContent(j, j)))
    
    if args.fixed:
      hs = [fIn.Get('FixedPlots/FixedXSec/Fixed%sXSec'%n) for n in ['Abs', 'Cex', 'OtherInel']]
      vals = [[h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)] for h in hs]
      bins = [[h.GetBinCenter(i) for i in range(1, h.GetNbinsX()+1)] for h in hs]
    
      grs = [RT.TGraph(len(vs), array('d', bs), array('d', vs)) for vs, bs in zip(vals, bins)]
    
    
    fOut.cd()
    corr.Write()
    cov.Write()
    if args.fixed:
      for i in range(3):
        grs[i].Write('gr_fixed_%s'%['Abs', 'Cex', 'OtherInel'][i])
  else:
    '''
    abs_hist = fIn.Get("PostFitXSec/PostFitAbsHist")
    abs_xs = []
    abs_ys = []
    abs_eyhs = []
    abs_eyls = []
    for i in range(1, abs_hist.GetNbinsX()+1):
      abs_xs.append(abs_hist.GetBinCenter(i))
      abs_ys.append(abs_hist.GetBinContent(i))
      abs_eyhs.append(0.)
      abs_eyls.append(0.)
    
    n_abs = len(abs_xs)
    print(n_abs)
    new_gr_abs = RT.TGraphAsymmErrors(
      n_abs,
      np.array(abs_xs),
      np.array(abs_ys),
      np.array([0.]*n_abs),
      np.array([0.]*n_abs),
      np.array(abs_eyls),
      np.array(abs_eyhs)
    )
    '''
    new_gr_abs = get_new_xsec(fIn, 'abs')
    new_gr_cex = get_new_xsec(fIn, 'cex')
    new_gr_other = get_new_xsec(fIn, 'other')
    
    n = new_gr_abs.GetN() + new_gr_cex.GetN() + new_gr_other.GetN()
    cov = RT.TH2D('xsec_cov', '', n, 0, n, n, 0, n)
    corr = RT.TH2D('xsec_corr', '', n, 0, n, n, 0, n)
  
    if len(args.add_files) > 0:
    #if len(sys.argv) > 3:
      add_covs = (
        [args.add_covs[0]]*len(args.add_files) if len(args.add_covs) == 1
        else args.add_covs
      )
      for af, ac in zip(args.add_files, add_covs):
        #extra_cov_file = RT.TFile(sys.argv[3], 'open')
        extra_cov_file = RT.TFile(af, 'open')
        #extra_cov = extra_cov_file.Get(sys.argv[4])
        extra_cov = extra_cov_file.Get(ac)
    
        print(af)
        for i in range(1, extra_cov.GetNbinsX()+1):
          for j in range(1, extra_cov.GetNbinsX()+1):
            cov.SetBinContent(i, j, cov.GetBinContent(i, j) + extra_cov.GetBinContent(i, j))
    
      #for i in range(1, extra_cov.GetNbinsX()+1):
      #  for j in range(1, extra_cov.GetNbinsX()+1):
      for i in range(1, cov.GetNbinsX()+1):
        for j in range(1, cov.GetNbinsX()+1):
          corr.SetBinContent(i, j, cov.GetBinContent(i, j)/math.sqrt(cov.GetBinContent(i, i)*cov.GetBinContent(j, j)))
  
  
    
    fOut.cd()
    new_gr_abs.Write('abs_xsec')
    new_gr_cex.Write('cex_xsec')
    new_gr_other.Write('other_xsec')
    cov.Write()
    corr.Write()
  fOut.Close();


def save(args):
  RT.gROOT.LoadMacro("~/style.C")
  RT.gROOT.ProcessLine('style()')
  f = RT.TFile.Open(args.i)
  cs = [
    'cAbsRebin', 'cCexRebin', 'cRejectedIntRebin',
    'cAPA2', 'cBadEvents', 'cMichelCut',
  ]
  for n in cs:
    c = f.Get(n)
    c.Draw()
    c.SaveAs(n + '.pdf') 
  f.Close()

if __name__ == '__main__':
  parser = ap()
  parser.add_argument('-o', type=str, default='')
  parser.add_argument('-i', type=str, required=True)
  parser.add_argument('--pretune', type=str, default=None)
  parser.add_argument('--add_files', nargs='+', default=[])
  parser.add_argument('--add_covs', nargs='+', default=[])
  parser.add_argument('--fixed', action='store_true')
  parser.add_argument('--nothrows', action='store_true')
  parser.add_argument('--save', action='store_true')
  parser.add_argument('--routine', type=str, default='normal',
                      choices=[
                        'normal', 'pars', 'save'
                      ])
  parser.add_argument('--npars', default=0, type=int)
  parser.add_argument('--noplotstyle', action='store_true')
  args = parser.parse_args()
  
  routines = {
    'normal':normal,
    'pars':fix_params_plot,
    'save':save,
  }
  routines[args.routine](args)
