import sys
import numpy as np
from argparse import ArgumentParser as ap

if __name__ == '__main__':

  parser = ap()
  parser.add_argument('-n', nargs=3, default=[5, 5, 5])
  parser.add_argument('-w', type=float, default=.33)
  parser.add_argument('-s', type=float, default=.10)
  args = parser.parse_args()

  centrals = args.w * np.random.randn(3) + 1.
  print(centrals)

  with open('tune.txt', 'w') as f:
    for i in range(3):
      
      n = args.n[i]
      si = (np.random.rand() - .5)*args.s/.5
      for j in range(args.n[i]):
        v = centrals[i] + (si/(n/2))*(-n/2 + j*(n-1)/n)
        if v < .1: v = .1
        f.write(f'{v}\n')

