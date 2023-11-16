import ROOT as RT
from array import array
from argparse import ArgumentParser as ap
from math import sqrt
import yaml

defaults = {
  'ymax':600.,
  'zmax':695.,
  'level':'michelremoving2',
  'track_xmax':350.,
  'track_ymin':40.,
  'track_ymax':560.,
  'track_zmin':40.,
  'track_zmax':655.,
  'plane2_theta_xz_min':60.,
  'plane2_theta_xz_max':120.,
  'plane2_theta_yz_min':80.,
  'plane2_theta_yz_max':100.,
}

def get_config(args):
  with open(args.c, 'r') as fin:
    config = yaml.safe_load(fin)

  for k,v in defaults.items():
    if k not in config.keys(): config[k] = v

  return config

if __name__ == '__main__':
  parser = ap()
  
  parser.add_argument('-o', type=str, required=True)
  parser.add_argument('-i', type=str, required=True)
  parser.add_argument('-c', type=str, required=True)

  parser.add_argument('--ymax', type=float, default=600.)
  parser.add_argument('--zmax', type=float, default=695.)
  #parser.add_argument('-e', type=str, required=True)
  args = parser.parse_args()

  config = get_config(args) 
  track_xmax = config['track_xmax']
  track_ymax = config['track_ymax']
  track_ymin = config['track_ymin']
  track_zmax = config['track_zmax']
  track_zmin = config['track_zmin']
  
  with open(args.i, 'r') as f:
    input_files = [l.strip() for l in f.readlines()]
  print(input_files)
  
  tree = RT.TChain(f'{config["level"]}/Event')
  for i in input_files:
    tree.AddFile(i)
  print(tree.GetEntries())
  
  fOut = RT.TFile(args.o, 'recreate')
  out_tree = RT.TTree('tree', '')
  dq_dx = array('d', [0])
  hit_x = array('d', [0])
  hit_y = array('d', [0])
  hit_z = array('d', [0])
  hit_plane_out = array('i', [0])
  x_corr_filter = array('i', [0])
  out_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
  out_tree.Branch('hit_x', hit_x, 'hit_x/D')
  out_tree.Branch('hit_y', hit_y, 'hit_y/D')
  out_tree.Branch('hit_z', hit_z, 'hit_z/D')
  out_tree.Branch('hit_plane', hit_plane_out, 'hit_plane/I')
  out_tree.Branch('x_corr_filter', x_corr_filter, 'x_corr_filter/I')
  
  nevent = 0
  for e in tree:
    print(nevent, end='\r')
    nevent += 1
    for i in range(0, e.cross_trks):
      #print('\tTrack %i'%i)
      track_index = i*9000
  
      x_cross_filter = e.trkstartx[i]*e.trkendx[i] < 0
      testneg = (e.trkstartx[i]<-track_xmax or e.trkendx[i]<-track_xmax)
      testpos = (e.trkstartx[i]>track_xmax or e.trkendx[i]>track_xmax)

      ###Check if the start or endpoint are in the FV. Both must be outside to enter sample
      #if ((abs(e.trkstartx[i]) < 350 and e.trkstarty[i] > 40 and
      #     e.trkstarty[i] < 560 and e.trkstartz[i] > 40 and e.trkstartz[i] < 655) or
      #    (abs(e.trkendx[i]) < 350 and e.trkendy[i] > 40 and e.trkendy[i] < 560 and
      #     e.trkendz[i] > 40 and e.trkendz[i] < 655)):
      if ((abs(e.trkstartx[i]) < track_xmax and e.trkstarty[i] > track_ymin and
           e.trkstarty[i] < track_ymax and e.trkstartz[i] > track_zmin and
           e.trkstartz[i] < track_zmax) or
          (abs(e.trkendx[i]) < track_xmax and e.trkendy[i] > track_ymin and
           e.trkendy[i] < track_ymax and e.trkendz[i] > track_zmin and
           e.trkendz[i] < track_zmax)):
        continue
  
      theta_xz_deg = 180./RT.TMath.Pi()*e.trackthetaxz[i]
      theta_yz_deg = 180./RT.TMath.Pi()*e.trackthetayz[i]
      for hit_plane in range(0, 3):
        ##Angular cuts
        if (hit_plane == 2 and
            ((abs(theta_xz_deg) > 60 and abs(theta_xz_deg) < 120) or
             (abs(theta_yz_deg) > 80 and abs(theta_yz_deg) < 100))):
          continue
  
        hit_index = hit_plane*3000
        hit_plane_out[0] = hit_plane
  
        nhits = e.ntrkhits[i*3 + hit_plane]
        for j in range(0, min(nhits, 3000)):
          hit_x[0] = e.trkhitx[track_index + hit_index + j]
          hit_y[0] = e.trkhity[track_index + hit_index + j]
          hit_z[0] = e.trkhitz[track_index + hit_index + j]
  
          if ((hit_y[0] > config['ymax']) or (hit_y[0] < 0.) or
              (hit_z[0] > config['zmax']) or (hit_z[0] < 0.)):
            continue
  
          ##negative x
          if (hit_x[0] > -360 and hit_x[0] < 0.):
            hit_pos = False
            if hit_plane == 1 and abs(180/RT.TMath.Pi()*e.trackthetaxz[i]) < 140: continue
            elif hit_plane == 0 and abs(180/RT.TMath.Pi()*e.trackthetaxz[i]) > 40: continue
          elif (e.trkhitx[track_index + hit_index + j] > 0. and e.trkhitx[track_index + hit_index + j] < 360.):
            hit_pos = True
            if hit_plane == 1 and abs(180/RT.TMath.Pi()*e.trackthetaxz[i]) > 40: continue
            elif hit_plane == 0 and abs(180/RT.TMath.Pi()*e.trackthetaxz[i]) < 140: continue
          else: continue
  
          dq_dx[0] = e.trkdqdx[track_index + hit_index + j]
  
          if (j < nhits-1) and (j > 0):
            if testpos and hit_pos:
              good_x_hit = (e.trkhitx[track_index + hit_index + j+1] > 0 and e.trkhitx[track_index + hit_index + j-1] > 0)
            elif testneg and not hit_pos:
              good_x_hit = (e.trkhitx[track_index + hit_index + j+1] < 0 and e.trkhitx[track_index + hit_index + j-1] < 0)
            else: good_x_hit = False
          else: good_x_hit = False
          x_corr_filter[0] = good_x_hit and x_cross_filter
  
          out_tree.Fill()
  fOut.cd()
  out_tree.Write()
  fOut.Close()
