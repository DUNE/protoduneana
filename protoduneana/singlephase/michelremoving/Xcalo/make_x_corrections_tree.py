import ROOT as RT
import yaml
import sys
from statistics import median
from argparse import ArgumentParser as ap
from array import array
import numpy as np
from math import log, sqrt

def get_recomb_factor(ef):
  E0 = 0.4867
  beta = .212
  alpha = .93
  rho = 1.38
  dedx0 = 1.9
  xsi0 = beta*dedx0/(rho*E0)
  rec0 = log(alpha + xsi0)/xsi0
  xsi = .212*1.9/(1.38*ef) 
  return (rec0*xsi)/log(alpha + xsi)

def get_efield(x, y, z, efield_pos_hists, efield_neg_hists):
  E0 = 0.4867
  if x > 0: hists = efield_pos_hists
  else: hists = efield_neg_hists

  #print(x, y, z)
  ex = E0 + E0*hists[0].Interpolate(x, y, z)
  ey = E0*hists[1].Interpolate(x, y, z)
  ez = E0*hists[2].Interpolate(x, y, z)
  return sqrt(ex**2 + ey**2 + ez**2)
 
def get_config(args):
  with open(args.c, 'r') as fin:
    config = yaml.safe_load(fin)

  #for k,v in defaults.items():
  #  if k not in config.keys(): config[k] = v

  return config


if __name__ == '__main__':
  parser = ap()
  parser.add_argument('-f', type=str, required=True)
  parser.add_argument('-c', type=str, required=True)
  parser.add_argument('-o', type=str, required=True)
  parser.add_argument('-e', type=str, required=True)
  parser.add_argument('--yz', type=str, required=True)
  parser.add_argument('--bin_max', type=int, default=10000)
  parser.add_argument('-n', type=int, default=-1)
  parser.add_argument('--nx', type=int, default=144)##replace with config
  parser.add_argument('--xmax', type=float, default=350.) ##Replace with config
  parser.add_argument('--run', type=int, default=0)
  args = parser.parse_args()
  
  
  config = get_config(args)
  efield_file = RT.TFile.Open(args.e, 'open')
  efield_pos_hists = [efield_file.Get('Reco_ElecField_%s_Pos'%i) for i in ['X', 'Y', 'Z']]
  efield_neg_hists = [efield_file.Get('Reco_ElecField_%s_Neg'%i) for i in ['X', 'Y', 'Z']]
  ymax = efield_pos_hists[0].GetYaxis().GetBinUpEdge(efield_pos_hists[0].GetNbinsY())
  zmax = efield_pos_hists[0].GetZaxis().GetBinUpEdge(efield_pos_hists[0].GetNbinsZ())
  print(ymax, zmax)

 
  yz_file = RT.TFile.Open(args.yz, 'open')
  yz_pos_correction = [yz_file.Get('correction_dqdx_ZvsY_positiveX_hist_%i'%i) for i in range(3)]
  yz_neg_correction = [yz_file.Get('correction_dqdx_ZvsY_negativeX_hist_%i'%i) for i in range(3)]
  
  
  
  
  ts = [RT.TChain('x_tree_%i'%i) for i in range(3)]
  with open(args.f, 'r') as f:
    files = [l.strip() for l in f.readlines() if '#' not in l]
  print(files)
  for f in files:
    for t in ts:
      t.AddFile(f)
      #print(t.GetEntries())
  for t in ts:
    print(t.GetEntries())
  '''
  t = RT.TChain('x_tree')
  with open(args.f, 'r') as f:
    files = [l.strip() for l in f.readlines() if '#' not in l]
  print(files)
  for f in files: t.AddFile(f)
  print(t.GetEntries())
  '''
  
  nx = config['nx']
  xmax = config['xmax']
  coverage_hists = [RT.TH1D('coverage_hist_%i'%i, '', nx, -xmax, xmax) for i in range(3)]
  correction_hists = [RT.TH1D('dqdx_X_correction_hist_%i'%i, '', nx, -xmax, xmax) for i in range(3)]
  dqdx_hists = [RT.TH1D('dqdx_hist_%i'%i, '', nx, -xmax, xmax) for i in range(3)]
  
  for h in coverage_hists: h.SetDirectory(0)
  
  global_median = []
  for i in range(3):
  
    all_dqdx = []
    n_cell = []
    dqdx = []
  
    for j in range(nx):
      dqdx.append(array('f', [0.])*args.bin_max)
      n_cell.append(0)
  
    t = ts[i]
    a = 0
    nevents = t.GetEntries()
    n_full = 0
    for e in t:
      if not a % 100000: print('%i/%i'%(a, nevents), end='\r')
      a += 1
      if args.n > 0 and a > args.n: break
    
      x = e.hit_x
      y = e.hit_y
      z = e.hit_z
    
      hit_plane = e.hit_plane
    
      x_bin = coverage_hists[hit_plane].FindBin(x) - 1
      if x_bin >= nx: continue
      if y < 0. or y > ymax: continue
      if z < 0. or z > zmax: continue
      if n_cell[x_bin] >= args.bin_max: continue
      recomb_factor = get_recomb_factor(get_efield(x, y, z, efield_pos_hists, efield_neg_hists))
      if x < 0.:
    
        yz_factor = yz_neg_correction[hit_plane].GetBinContent(
            yz_neg_correction[hit_plane].FindBin(z, y))
      else:
        yz_factor = yz_pos_correction[hit_plane].GetBinContent(
            yz_pos_correction[hit_plane].FindBin(z, y))
    
      #print(hit_plane, x_bin, n_cell[x_bin]) 
      dqdx[x_bin][n_cell[x_bin]] = e.dq_dx*yz_factor*recomb_factor
    
      n_cell[x_bin] += 1
    
      if n_cell[x_bin] == args.bin_max:
        print(x_bin, 'full')
        n_full += 1
    
      all_dqdx.append(e.dq_dx*yz_factor*recomb_factor)
      coverage_hists[hit_plane].Fill(x)
  
    this_global_median = median(all_dqdx) if len(all_dqdx) >= 5 else 1.
    print(this_global_median)
    global_median.append(this_global_median)
    for j in range(nx):
      #print(i, j, end='\r')
      nonzero = [i for i in dqdx[j] if i > 0.]
  
      nonzero_med = median(nonzero)  if (n_cell[j] >= 5) else 1.
      f = this_global_median/nonzero_med if (n_cell[j] >= 5) else 1.
      dqdx_hists[i].SetBinContent(j+1, nonzero_med)
      correction_hists[i].SetBinContent(j+1, f)
  
  
  
      #print(i, j, k, med_neg, med_pos, (global_median_neg[i]/med_neg), (global_median_pos[i]/med_pos))
  
  
  with open('globalmedians_run%i.csv'%args.run, 'w') as fcsv:
    fcsv.write('channel,tv,norm,norm_err\n')
    for i in range(3): 
      fcsv.write('%i,%i,%f,%f\n'%(i, args.run, global_median[i], global_median[i]/10.))
  
  fOut = RT.TFile(args.o, 'recreate')
  for h in correction_hists: h.Write()
  for h in coverage_hists: h.Write()
  for h in dqdx_hists: h.Write()
  fOut.Close()
