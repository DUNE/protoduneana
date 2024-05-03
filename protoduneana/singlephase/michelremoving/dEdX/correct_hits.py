import ROOT as RT
import sys
from statistics import median
from argparse import ArgumentParser as ap
from array import array
import numpy as np

parser = ap()


parser.add_argument('-f', type=str, required=True)
parser.add_argument('-o', type=str, required=True)
parser.add_argument('-x', type=str, required=True)
parser.add_argument('--yz', type=str, required=True)
args = parser.parse_args()

x_corr_file = RT.TFile(args.x, 'open')
yz_corr_file = RT.TFile(args.yz, 'open')
x_corr_hists = [x_corr_file.Get('dqdx_X_correction_hist_%i'%i) for i in range(0, 3)]
#x_corr_hists = [x_corr_file.Get('correction_dqdx_hist_%i'%i) for i in range(0, 3)]
yz_corr_hists_neg = [yz_corr_file.Get('correction_dqdx_ZvsY_negativeX_hist_%i'%i) for i in range(0, 3)]
yz_corr_hists_pos = [yz_corr_file.Get('correction_dqdx_ZvsY_positiveX_hist_%i'%i) for i in range(0, 3)]

t = RT.TChain('dedx_tree')
with open(args.f, 'r') as f:
  files = [l.strip() for l in f.readlines() if '#' not in l]
print(files)
for f in files:
  t.AddFile(f)
nentries = t.GetEntries()
print(nentries)

fOut = RT.TFile(args.o, 'recreate')
new_tree = RT.TTree('dedx_tree', '')
corrected_dq_dx = array('d', [0.])
new_tree.Branch('corrected_dq_dx', corrected_dq_dx, 'corrected_dq_dx/D')

dq_dx = array('d', [0.])
hit_plane = array('i', [0])
hit_x = array('d', [0.])
hit_y = array('d', [0.])
hit_z = array('d', [0.])
efield = array('d', [0.])
resrange = array('d', [0.])
range_bin = array('i', [0])
ke = array('d', [0.])

new_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
new_tree.Branch('hit_plane', hit_plane, 'hit_plane/I')
new_tree.Branch('hit_x', hit_x, 'hit_x/D')
new_tree.Branch('hit_y', hit_y, 'hit_y/D')
new_tree.Branch('hit_z', hit_z, 'hit_z/D')
new_tree.Branch('efield', efield, 'efield/D')
new_tree.Branch('resrange', resrange, 'resrange/D')
new_tree.Branch('range_bin', range_bin, 'range_bin/I')
new_tree.Branch('ke', ke, 'ke/D')

a = 0
for e in t:
  if not a % 1000: print('%i/%i'%(a, nentries), end='\r')  
  e.hit_plane
  
  x_bin = x_corr_hists[e.hit_plane].FindBin(e.hit_x)
  #if (x_bin != int((e.hit_x[0] + 360.)/5.)+1):
  #  print("ERROR")
  #  break 
  Cx = x_corr_hists[e.hit_plane].GetBinContent(x_bin)

  if e.hit_x < 0.:
    yz_hist = yz_corr_hists_neg[e.hit_plane]
  else:
    yz_hist = yz_corr_hists_pos[e.hit_plane]
  yz_bin = yz_hist.FindBin(e.hit_z, e.hit_y) 
  Cyz = yz_hist.GetBinContent(yz_bin)

  corrected_dq_dx[0] = e.dq_dx*Cx*Cyz
  dq_dx[0] = e.dq_dx
  hit_plane[0] = e.hit_plane
  hit_x[0] = e.hit_x
  hit_y[0] = e.hit_y
  hit_z[0] = e.hit_z
  efield[0] = e.efield
  resrange[0] = e.resrange
  range_bin[0] = e.range_bin
  ke[0] = e.ke
  new_tree.Fill()

fOut.cd()
new_tree.Write()
fOut.Close()
