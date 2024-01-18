import ROOT as RT
from math import sqrt
from argparse import ArgumentParser as ap
from array import array
RT.gROOT.SetBatch()
RT.gStyle.SetOptStat(0)

parser = ap()

parser.add_argument("-i", type=str, help='Input file', default = "")
parser.add_argument("-p", type=str, help='Plus file', default = "")
parser.add_argument("-m", type=str, help='Minus file', default = "")
parser.add_argument("-o", type=str, help='Output file', default='beam_shift_results.root')
parser.add_argument('--uncorrelated', action='store_true')
parser.add_argument('--post', action='store_true')
args = parser.parse_args()


names = ['grXSecThrowAbsUnderflow', 'grXSecThrowCexUnderflow',
         'grXSecThrowOtherInelUnderflow']

fN = RT.TFile(args.i, 'open')
if args.post:
  hNs = [fN.Get('PostFitXSec/PostFitAbsXSec'),
         fN.Get('PostFitXSec/PostFitCexXSec'),
         fN.Get('PostFitXSec/PostFitOtherInelXSec')] 
  xs = [array('d', [h.GetBinCenter(i) for i in range(1, h.GetNbinsX()+1)]) for h in hNs]
  yNs = [array('d', [h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)]) for h in hNs]
  gNs = [RT.TGraph(len(x), x, y) for x,y in zip(xs, yNs)]


else:
  gNs = [fN.Get('Throws/%s'%n) for n in names]

  xs = [array('d', [gN.GetX()[i] for i in range(gN.GetN())]) for gN in gNs]
fP = RT.TFile(args.p, 'open')
fM = RT.TFile(args.m, 'open')

if args.post:
  hPs = [fP.Get('PostFitXSec/PostFitAbsXSec'),
         fP.Get('PostFitXSec/PostFitCexXSec'),
         fP.Get('PostFitXSec/PostFitOtherInelXSec')]

  hMs = [fM.Get('PostFitXSec/PostFitAbsXSec'),
         fM.Get('PostFitXSec/PostFitCexXSec'),
         fM.Get('PostFitXSec/PostFitOtherInelXSec')]

  yPs = [array('d', [hP.GetBinContent(i) for i in range(1, hP.GetNbinsX()+1)]) for hP in hPs]
  yMs = [array('d', [hM.GetBinContent(i) for i in range(1, hM.GetNbinsX()+1)]) for hM in hMs]

  print(yPs)
  print(yMs)

  gPs = [RT.TGraph(len(x), x, y) for x,y in zip(xs, yPs)]
  gMs = [RT.TGraph(len(x), x, y) for x,y in zip(xs, yMs)]
else:
  gPs = [fP.Get('Throws/%s'%n) for n in names]
  gMs = [fM.Get('Throws/%s'%n) for n in names]

fOut = RT.TFile(args.o, 'recreate')

leg = RT.TLegend()
for i in range(0, len(gNs)):
  gN = gNs[i]
  gP = gPs[i]
  gM = gMs[i]

  themax = max([max([gN.GetY()[j], gP.GetY()[j], gM.GetY()[j]]) for j in range(gNs[i].GetN())])
  print('Max:', themax)

  c = RT.TCanvas('c%i'%i, '')
  c.SetTicks()

  gP.SetLineColor(RT.kRed)
  gP.SetMarkerColor(RT.kRed)
  gM.SetLineColor(RT.kBlue)
  gM.SetMarkerColor(RT.kBlue)

  gN.SetMaximum(themax)
  gN.SetMinimum(0.)
  gN.Draw('AP')
  gN.Draw('P same')
  gP.SetMarkerStyle(20)
  gP.Draw('P same')
  gM.SetMarkerStyle(20)
  gM.Draw('P same')
  gN.SetTitle(';Cross Section Bin;#sigma (mb)');

  if i == 0:
    leg.AddEntry(gN, 'Nominal Fit')
    leg.AddEntry(gP, '+1#sigma Fit')
    leg.AddEntry(gM, '-1#sigma Fit')
    leg.Write('leg')
  c.Write()
  gP.Write(f'plus_{i}')
  gM.Write(f'minus_{i}')
  gN.Write(f'nom_{i}')

nom_vals = [gN.GetY()[i] for gN in gNs for i in range(0, gN.GetN())]
plus_vals = [gP.GetY()[i] for gP in gPs for i in range(0, gP.GetN())]
minus_vals = [gM.GetY()[i] for gM in gMs for i in range(0, gM.GetN())]
plus_diffs = [n - p for n, p in zip(nom_vals, plus_vals)]
minus_diffs = [n - p for n, p in zip(nom_vals, minus_vals)]
plus_cov_vals = [[di*dj for di in plus_diffs] for dj in plus_diffs]
minus_cov_vals = [[di*dj for di in minus_diffs] for dj in minus_diffs]

plus_cov = RT.TH2D('plus_cov', '', len(nom_vals), 0, len(nom_vals), len(nom_vals), 0, len(nom_vals))
minus_cov = RT.TH2D('minus_cov', '', len(nom_vals), 0, len(nom_vals), len(nom_vals), 0, len(nom_vals))
ave_cov = RT.TH2D('ave_cov', '', len(nom_vals), 0, len(nom_vals), len(nom_vals), 0, len(nom_vals))
ave_var = RT.TH1D('ave_var', '', len(nom_vals), 0, len(nom_vals))

for i in range(0, len(nom_vals)):
  ave_var.SetBinContent(i+1, sqrt((plus_cov_vals[i][i] + minus_cov_vals[i][i])/2.))
  for j in range(0, len(nom_vals)):
    plus_cov.SetBinContent(i+1, j+1, plus_cov_vals[i][j])
    minus_cov.SetBinContent(i+1, j+1, minus_cov_vals[i][j])
    if args.uncorrelated:
      if i == j:
        ave_cov.SetBinContent(i+1, j+1, (plus_cov_vals[i][j] + minus_cov_vals[i][j])/2.)
      else:
        ave_cov.SetBinContent(i+1, j+1, 0.)
    else:
      ave_cov.SetBinContent(i+1, j+1, (plus_cov_vals[i][j] + minus_cov_vals[i][j])/2.)

plus_cov.SetTitle(';Cross Section Bin;Cross Section Bin')
plus_cov.Write()
minus_cov.SetTitle(';Cross Section Bin;Cross Section Bin')
minus_cov.Write()
ave_cov.SetTitle(';Cross Section Bin;Cross Section Bin')
ave_var.SetTitle(';Cross Section Bin;Uncertainty [mb]')
ave_var.SetTitleSize(.05, 'XY')
ave_var.SetTitleOffset(.8, 'XY')
ave_cov.Write()
ave_var.Write()
c_var = RT.TCanvas()
c_var.SetTicks()
ave_var.Draw()
c_var.Write('c_var')

ave_corr = ave_cov.Clone('ave_corr')
for i in range(0, len(nom_vals)):
  val_i = ave_cov.GetBinContent(i+1, i+1)
  for j in range(0, len(nom_vals)):
    val_j = ave_cov.GetBinContent(j+1, j+1)
    ave_corr.SetBinContent(i+1, j+1,
                           ave_cov.GetBinContent(i+1, j+1)/sqrt(val_i*val_j))

ave_corr.SetTitleSize(.05, 'XY')
ave_corr.SetTitleOffset(.8, 'XY')
ave_corr.Write()
c_corr = RT.TCanvas()
c_corr.SetTicks()
ave_corr.Draw('colz')
c_corr.Write('c_corr')

fOut.Close()
