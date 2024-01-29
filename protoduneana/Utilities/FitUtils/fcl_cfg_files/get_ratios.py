import ROOT as RT
import numpy as np
import sys

if __name__ == '__main__':
  filename = sys.argv[1]
  retune = len(sys.argv) > 2
  if filename[:5] == '/pnfs':
    filename = filename.replace('/pnfs', 'root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr')
    print(filename)
  fIn = RT.TFile.Open(filename)

  all_vals = []
  for n in ['Abs', 'Cex', 'OtherInel']:
    vals = []
    #f.write('1.\n')
    hPost = fIn.Get(f"PostFitXSec/PostFit{n}XSec")
    hPre = fIn.Get(f"PreFitXSec/PreFit{n}XSec")
    for i in range(1, hPost.GetNbinsX()+1):
      vals.append(hPost.GetBinContent(i)/hPre.GetBinContent(i))
    vals = [np.mean(vals[:2])] + vals + [np.mean(vals[-2:])]
    all_vals.append(vals)

  if (retune):
    old_tune = fIn.Get('TunePars/tune_pars')
    ipar = 0
    for vals in all_vals:
      for i in range(len(vals)):
        vals[i] *= old_tune[ipar]
        ipar += 1

  fIn.Close()

  with open('tune.txt', 'w') as f:
    for vals in all_vals:
      for v in vals:
        f.write(f'{v}\n')
      
  
  '''with open('tune.txt', 'w') as f:
    for n in ['Abs', 'Cex', 'OtherInel']:
      f.write('1.\n')
      hPost = fIn.Get(f"PostFitXSec/PostFit{n}XSec")
      hPre = fIn.Get(f"PreFitXSec/PreFit{n}XSec")
      for i in range(1, hPost.GetNbinsX()+1):
        f.write(f'{hPost.GetBinContent(i)/hPre.GetBinContent(i)}\n')
      f.write('1.\n')'''

  with open('file.txt', 'w') as f:
    f.write(filename)

  #fIn.Close()

