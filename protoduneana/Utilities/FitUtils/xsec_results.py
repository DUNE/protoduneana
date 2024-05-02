import ROOT as RT 
from array import array
from argparse import ArgumentParser as ap
import sys
from math import sqrt


parser = ap()

parser.add_argument('-o', type=str, required=True)
parser.add_argument('--total', action='store_true')
parser.add_argument('-f', type=str, required=True)
parser.add_argument('-x', type=str, default='/exp/dune/data/users/calcuttj/old_data2/PiAnalysis_G4Prediction/thresh_abscex_xsecs.root')
parser.add_argument('--genie', type=str, default=None)
parser.add_argument('--xt', type=str, default=None)
parser.add_argument('--al', type=int, default=0)
parser.add_argument('--ah', type=int, default=-1)
parser.add_argument('--cl', type=int, default=0)
parser.add_argument('--ch', type=int, default=-1)
parser.add_argument('--ol', type=int, default=0)
parser.add_argument('--oh', type=int, default=-1)
parser.add_argument('-t', action='store_true')
parser.add_argument('-p', action='store_true')
parser.add_argument('-m', action='store_true', help='Set max to ultimate max of xsec')
parser.add_argument('--LADS', type=int, help='0: Use results from Kotlinski. 1: Use results from Rowntree. 2: Both', default = -1)
parser.add_argument('--add', action='store_true', help='Set error bars to cov') 
parser.add_argument('--fixed', action='store_true')
parser.add_argument('--nochi2', action='store_true')
parser.add_argument('-v', action='store_true')

args = parser.parse_args()

RT.gROOT.LoadMacro("~/protoDUNEStyle.C")
RT.gROOT.SetBatch();
RT.gStyle.SetOptStat(00000)
RT.gStyle.SetErrorX(1.e-4)
RT.gStyle.SetTitleAlign(33)
RT.gStyle.SetTitleX(.5)
tt = RT.TLatex();
tt.SetNDC();



if args.p:
  t_prelim = RT.TLatex()
  t_prelim.SetNDC()
  t_prelim.SetTextColor(17)
  t_prelim.SetTextSize(0.1);
  t_prelim.SetTextAngle(26.15998);
  t_prelim.SetLineWidth(2);

t_xsecs = RT.TLatex()
t_xsecs.SetNDC()
t_xsecs.SetTextAngle(90.);

fResults = RT.TFile(args.f, 'open')
cov_hist = fResults.Get('xsec_cov')
gr_abs = fResults.Get('abs_xsec')
gr_cex = fResults.Get('cex_xsec')
gr_other = fResults.Get('other_xsec')
result_xsecs = [gr_abs, gr_cex, gr_other]
#for x in result_xsecs:
#  for i in range(x.GetN()):
#    x.SetPointY(i, x.GetY()[i]*0.9628*0.9628)
if args.add:
  count = 1
  for x in result_xsecs:
    for i in range(0, x.GetN()):
      print(count, count, sqrt(cov_hist.GetBinContent(count, count)))
      x.SetPointEYhigh(i, sqrt(cov_hist.GetBinContent(count, count)))
      x.SetPointEYlow(i, sqrt(cov_hist.GetBinContent(count, count)))
      count += 1

#if args.fixed:
#  fixed_xsecs = []
#  for n in ['Abs', 'Cex', 'OtherInel']:

fG4 = RT.TFile(args.x, 'open')

n_abs = gr_abs.GetN()
n_cex_abs = n_abs + gr_cex.GetN()
n_total = n_cex_abs + gr_other.GetN()

g4_xsecs = [fG4.Get('abs_KE').Clone(), fG4.Get('cex_KE').Clone()]
#other_grs = [fG4.Get('dcex_KE'), fG4.Get('inel_KE'), fG4.Get('prod_KE')]
#total = [other_grs[0].GetY()[i] + other_grs[1].GetY()[i] + other_grs[2].GetY()[i] for i in range(0, other_grs[0].GetN())]
#xs = [x for x in other_grs[0].GetX()]
#g4_xsecs.append(RT.TGraph(len(xs), array('d', xs), array('d', total)))
g4_xsecs.append(fG4.Get('other_KE').Clone())


g4_maxes = [max([y for y in g.GetY()]) for g in g4_xsecs]
print(g4_maxes)


genie_maxes = []
if args.genie:
  fGenie = RT.TFile(args.genie, 'open')
  genie_xsecs = [fGenie.Get('abs').Clone(), fGenie.Get('cex').Clone(), fGenie.Get('other').Clone()]
  genie_maxes = [max([y for y in g.GetY()]) for g in genie_xsecs]
  fGenie.Close()


if args.xt:
  fG4Thresh = RT.TFile(args.xt, 'open')
  g4_xsecs_thresh = [fG4Thresh.Get('abs_KE').Clone(), fG4Thresh.Get('cex_KE').Clone()]
  other_grs_thresh = [fG4Thresh.Get('dcex_KE'), fG4Thresh.Get('inel_KE'), fG4Thresh.Get('prod_KE')]
  total_thresh = [other_grs_thresh[0].GetY()[i] + other_grs_thresh[1].GetY()[i] + other_grs_thresh[2].GetY()[i] for i in range(0, other_grs_thresh[0].GetN())]
  xs_thresh = [x for x in other_grs_thresh[0].GetX()]
  g4_xsecs_thresh.append(RT.TGraph(len(xs), array('d', xs_thresh), array('d', total_thresh)))


result_maxes = []
for i in range(0, len(result_xsecs)):
  g = result_xsecs[i]
  eys = [g.GetErrorYhigh(j) for j in range(0, g.GetN())]
  ys = [g.GetY()[j] for j in range(0, g.GetN())]
  tots = [y + e for y,e in zip(ys, eys)]
  result_maxes.append(max(tots))

names = ['abs', 'cex', 'other']
titles = ['Absorption', 'Charge Exchange', 'Other']

#the_max = max([i for i in g4_maxes] + [i for i in result_maxes])
the_max = max(g4_maxes + result_maxes + genie_maxes)
fOut = RT.TFile(args.o, 'recreate')
for i in [0, 1, 2]:
  c = RT.TCanvas('c%s'%names[i], '')
  c.SetTicks()
  g4_xsecs[i].SetMinimum(0.)
  result_xsecs[i].SetMinimum(0.)
  result_xsecs[i].SetLineWidth(2)

  if args.m:
    g4_xsecs[i].SetMaximum(1.05*the_max)
  else:
    g4_xsecs[i].SetMaximum(1.05*max([g4_maxes[i], result_maxes[i]]))

  g4_xsecs[i].SetLineColor(RT.kRed)
  g4_xsecs[i].SetTitle('%s;Kinetic Energy [MeV];#sigma [mb]'%titles[i])
  g4_xsecs[i].GetXaxis().CenterTitle()
  g4_xsecs[i].GetYaxis().CenterTitle()
  g4_xsecs[i].GetXaxis().SetRangeUser(0., 999.)
  g4_xsecs[i].Draw('AC')
  g4_xsecs[i].SetLineWidth(2)

  if args.genie:
    genie_xsecs[i].SetLineColor(RT.kBlue)
    genie_xsecs[i].SetTitle('%s;Kinetic Energy [MeV];#sigma [mb]'%titles[i])
    genie_xsecs[i].GetXaxis().CenterTitle()
    genie_xsecs[i].GetYaxis().CenterTitle()
    genie_xsecs[i].GetXaxis().SetRangeUser(0., 999.)
    genie_xsecs[i].Draw('C same')
    genie_xsecs[i].SetLineWidth(2)


  if args.xt:
    g4_xsecs_thresh[i].SetLineColor(RT.kRed)
    g4_xsecs_thresh[i].SetLineWidth(2)
    g4_xsecs_thresh[i].SetLineStyle(9)
    g4_xsecs_thresh[i].Draw('C same')

  result_xsecs[i].Draw('pez same')
  result_xsecs[i].SetMarkerStyle(20)
  result_xsecs[i].SetMarkerColor(RT.kBlack)
  if i == 0:
    leg = RT.TLegend()
    leg.SetFillStyle(0)
    leg.SetLineWidth(0)
    if args.xt:
      leg.AddEntry(g4_xsecs[i], 'Geant4 10.6 Thresholds', 'l')
    elif args.v:
      leg.AddEntry(g4_xsecs[i], 'Geant4 10.6 Varied', 'l')
    else:
      leg.AddEntry(g4_xsecs[i], 'Geant4 10.6', 'l')

    if args.genie:
      leg.AddEntry(genie_xsecs[i], 'Genie', 'l')
    if args.xt:
      leg.AddEntry(g4_xsecs_thresh[i], 'Geant4 10.6 No Thresholds', 'l')
    leg.AddEntry(result_xsecs[i], 'ProtoDUNE-SP', 'pez')
    leg.Draw()
  tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  if args.p:
    t_prelim.DrawLatex(0.33, .5, 'Preliminary')
  RT.gPad.RedrawAxis()

  c.Write()

xsec_corr = fResults.Get('xsec_corr')
#xsec_corr.SetTitle(';Cross Section Bin;Cross Section Bin')
xsec_corr.SetTitle(';;')
cCorr = RT.TCanvas('cCorr', '')
xsec_corr.GetXaxis().CenterTitle()
xsec_corr.GetYaxis().CenterTitle()
xsec_corr.GetXaxis().SetNdivisions(13)
xsec_corr.GetYaxis().SetNdivisions(13)
xsec_corr.SetMaximum(1.)
xsec_corr.SetMinimum(-1.)
xsec_corr.Draw('colz')
RT.gStyle.SetPaintTextFormat('.2f')
xsec_corr.SetMarkerSize(1.5)
if args.t:
  xsec_corr.Draw('text same')
tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
tt.DrawLatex(0.23, 0.025, "#bf{Abs}")
tt.DrawLatex(0.50, 0.025, "#bf{Cex}")
tt.DrawLatex(0.73, 0.025, "#bf{Other}")
t_xsecs.DrawLatex(0.030, 0.25, "#bf{Abs}")
t_xsecs.DrawLatex(0.030, 0.50, "#bf{Cex}")
t_xsecs.DrawLatex(0.030, 0.73, "#bf{Other}")

#corr_labels = [450, 550, 650, 750, 850,
#  450, 550, 650, 750,
#  450, 550, 650, 750
#]
#xsec_corr.GetXaxis().SetLabelOffset(.01)
#xsec_corr.GetYaxis().SetLabelOffset(.01)
#for i in range(0, len(corr_labels)):
#  xsec_corr.GetXaxis().SetBinLabel(i+1, str(corr_labels[i]))
#  xsec_corr.GetYaxis().SetBinLabel(i+1, str(corr_labels[i]))

l1 = RT.TLine(n_abs, 0., n_abs, n_total)
l2 = RT.TLine(n_cex_abs, 0., n_cex_abs, n_total)

l3 = RT.TLine(0., n_abs, n_total, n_abs)
l4 = RT.TLine(0., n_cex_abs, n_total, n_cex_abs)

l1.Draw()
l2.Draw()
l3.Draw()
l4.Draw()

if args.p:
  t_prelim.DrawLatex(0.33, .5, 'Preliminary')
cCorr.Write()

for x,n in zip(result_xsecs, names):
  print(n, [i for i in x.GetY()])
  print(n, [x.GetErrorYhigh(i) for i in range(0, x.GetN())])

###Total
if args.total:
  total_xs = [i for i in result_xsecs[0].GetX()][args.al:]
  total_results = []
  total_results += [i for i in result_xsecs[0].GetY()][args.al:]
  total_errs = []
  for i in range(0, result_xsecs[0].GetN()):
    total_errs.append(result_xsecs[0].GetEYhigh()[i]**2)
  total_errs = total_errs[args.al:]
  #print(total_errs)
  for i in range(0, len(total_results)):
    total_results[i] += result_xsecs[1].GetY()[args.cl+i]
    total_results[i] += result_xsecs[2].GetY()[args.ol+i]
  
    total_errs[i] += result_xsecs[1].GetEYhigh()[args.cl+i]**2
    total_errs[i] += result_xsecs[2].GetEYhigh()[args.ol+i]**2
  #total_errs = [sqrt(e) for e in total_errs]
  print(total_results)
  print([sqrt(i) for i in total_errs])
  
  #cov_hist = fResults.Get('xsec_cov')
  total_cov_hist = RT.TH2D("total_cov", "", len(total_results), 0, len(total_results), len(total_results), 0, len(total_results))
  for i in range(0, 3):
    for j in range(0, 3):
      for k in range(0, len(total_results)):
        for l in range(0, len(total_results)):
          bin_i = i*len(total_results) + k + 1 + (args.al if i == 0 else 0)
          bin_j = j*len(total_results) + l + 1 + (args.al if j == 0 else 0)
          #print(bin_i, bin_j)
          total_cov_hist.SetBinContent(k+1, l+1, total_cov_hist.GetBinContent(k+1, l+1) + cov_hist.GetBinContent(bin_i, bin_j))
  total_cov_hist.Write("total_cov")
  
  added_errs = [sqrt(total_cov_hist.GetBinContent(i+1, i+1)) for i in range(0, len(total_results))]
  gr_total = RT.TGraphErrors(len(total_results), array('d', total_xs), array('d', total_results),
                             array('d', [0.]*len(total_results)), array('d', added_errs))
  gr_total.Write('total_gr')
  
  total_g4 = fG4.Get('total_inel_KE')
  total_g4.SetMinimum(0.)
  total_g4.SetLineWidth(2)
  gr_total.SetLineWidth(2)
  total_g4.SetLineColor(RT.kRed)
  total_g4.SetTitle('Total Inelastic;Kinetic Energy [MeV];#sigma [mb]')
  
  c = RT.TCanvas('cTotal', '')
  total_g4.Draw('AC')
  gr_total.Draw('pez same')
  total_g4.GetXaxis().CenterTitle()
  total_g4.GetYaxis().CenterTitle()
  total_g4.GetXaxis().SetRangeUser(0., 999.)
  tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  leg = RT.TLegend()
  leg.AddEntry(total_g4, 'Geant4 10.6', 'l')
  leg.AddEntry(gr_total, 'ProtoDUNE-SP', 'pez')
  leg.Draw()
  if args.p:
    t_prelim.DrawLatex(0.33, .5, 'Preliminary')
  c.Write()
  
  total_corr = total_cov_hist.Clone('total_corr')
  for i in range(1, total_corr.GetNbinsX()+1):
    for j in range(1, total_corr.GetNbinsX()+1):
      total_corr.SetBinContent(i, j, total_corr.GetBinContent(i, j)/(gr_total.GetEY()[i-1]*gr_total.GetEY()[j-1]))
  c = RT.TCanvas('cTotalCorr', '')
  total_corr.Draw('colz')
  total_corr.GetXaxis().SetNdivisions(4)
  total_corr.GetYaxis().SetNdivisions(4)
  total_corr.GetXaxis().CenterTitle()
  total_corr.GetYaxis().CenterTitle()
  total_corr.SetTitle(';Cross Section Bin;Cross Section Bin;')
  total_corr.GetZaxis().SetRangeUser(-1., 1.)
  total_corr.SetMarkerSize(1.5)
  if args.t:
    total_corr.Draw('text same')
  tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");
  if args.p:
    t_prelim.DrawLatex(0.33, .5, 'Preliminary')
  c.Write()

##WWith LADS
#Rowntree
LADS_abs = [393., 366., 282.]
LADS_errs = [21., 22., 21.]
LADS_xs = [118., 162., 239.]
LADS_1 = RT.TGraphErrors(len(LADS_abs), array('d', LADS_xs), array('d', LADS_abs),
                               array('d', [0.]*len(LADS_abs)),
                               array('d', LADS_errs))

#Kotlinski
LADS_abs = [180., 320., 351., 283., 225]
LADS_errs = [43., 65., 40., 28., 17.]
LADS_xs = [70., 118., 162., 239., 330.]
LADS_0 = RT.TGraphErrors(len(LADS_abs), array('d', LADS_xs), array('d', LADS_abs),
                               array('d', [0.]*len(LADS_abs)),
                               array('d', LADS_errs))

c = RT.TCanvas("cabs_lads")
c.SetTicks()
g4_xsecs[0].Draw('AC')
g4_xsecs[0].SetLineWidth(2)
if args.genie:
  genie_xsecs[0].SetLineWidth(2)
  genie_xsecs[0].Draw('C same')
if args.xt:
  g4_xsecs_thresh[0].SetLineColor(RT.kRed)
  g4_xsecs_thresh[0].SetLineWidth(2)
  g4_xsecs_thresh[0].SetLineStyle(9)
  g4_xsecs_thresh[0].Draw('C same')
result_xsecs[0].Draw('pez same')
tt.DrawLatex(0.10,0.94,"#bf{DUNE:ProtoDUNE-SP}");

LADS_0.SetLineWidth(2)
LADS_0.SetMarkerStyle(25)
LADS_1.SetLineWidth(2)
LADS_1.SetMarkerStyle(26)

LADS_0.SetLineColor(4)
LADS_0.SetMarkerColor(4)
LADS_1.SetLineColor(8)
LADS_1.SetMarkerColor(8)

if args.LADS in [0, 2]:
  LADS_0.Draw('pez same')

if args.LADS in [1, 2]:
  LADS_1.Draw('pez same')

KEs = [450, 550, 650, 750, 850]
lines = []
#for i in range(0, len(KEs)):
#  if i == 0: line = '%i & %.2f $\pm$ %.2f & - & - & -\\\\'%(KEs[i], result_xsecs[0].GetY()[i], result_xsecs[0].GetEYhigh()[i])
#  else: line = '%i & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f\\\\'%(KEs[i], result_xsecs[0].GetY()[i], result_xsecs[0].GetEYhigh()[i], result_xsecs[1].GetY()[i-1], result_xsecs[1].GetEYhigh()[i-1], result_xsecs[2].GetY()[i-1], result_xsecs[2].GetEYhigh()[i-1], total_results[i-1], sqrt(total_errs[i-1]))
#  print(line)

all_xs = [result_xsecs[j].GetX()[i] for j in range(0, 3) for i in range(0, result_xsecs[j].GetN())]
all_ys = [result_xsecs[j].GetY()[i] for j in range(0, 3) for i in range(0, result_xsecs[j].GetN())]
all_g4s = [g4_xsecs[j].Eval(result_xsecs[j].GetX()[i]) for j in range(0, 3) for i in range(0, result_xsecs[j].GetN())]
print(all_xs)
print(all_ys)
print(all_g4s)
xsec_cov_mat = RT.TMatrixD(cov_hist.GetNbinsX(), cov_hist.GetNbinsX())
for i in range(1, cov_hist.GetNbinsX()+1):
  for j in range(1, cov_hist.GetNbinsX()+1):
    xsec_cov_mat[i-1][j-1] = cov_hist.GetBinContent(i, j) 
xsec_cov_mat = xsec_cov_mat.Invert()

xsec_chi2 = 0.
for i in range(0, len(all_xs)):
  for j in range(0, len(all_xs)):
    print((all_ys[i] - all_g4s[i])*xsec_cov_mat[i][j]*(all_ys[j] - all_g4s[j]), xsec_cov_mat[i][j])
    xsec_chi2 += (all_ys[i] - all_g4s[i])*xsec_cov_mat[i][j]*(all_ys[j] - all_g4s[j])
print('cross section chi2: %.2f'%xsec_chi2)
leg = RT.TLegend()
leg.AddEntry(g4_xsecs[0], 'Geant4 10.6' if not args.xt else 'Geant4 10.6 Thresholds', 'l')
if args.xt:
  leg.AddEntry(g4_xsecs_thresh[0], 'Geant4 10.6 No Thresholds', 'l')
leg.AddEntry(result_xsecs[0], 'ProtoDUNE-SP', 'pez')



if args.LADS >= 0:
  if args.LADS in [0, 2]:
    leg.AddEntry(LADS_0, "Kotlinski et al. (2000)", 'pez')
  if args.LADS in [1, 2]:
    leg.AddEntry(LADS_1, "Rowntree et al. (1999)", 'pez')

if not args.nochi2:
  leg.AddEntry('', '#chi^{2} = %.2f'%xsec_chi2, '')
leg.Draw('same')
if args.p:
  t_prelim.DrawLatex(0.33, .5, 'Preliminary')
RT.gPad.RedrawAxis()
c.Write()

fOut.Close()
