import ROOT as RT
from array import array
from argparse import ArgumentParser as ap
from statistics import median
from math import sqrt
import cProfile
import pstats
import numpy as np

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
  'peakT_min':100.,
  'peakT_max':5900.,
  'track_len_max':700.,
  'track_len_min':100.,
  'track_ediv_min':226.,
  'track_ediv_max':236.,
  'track_ediv2_min':456.,
  'track_ediv2_max':472.,
  'do_true_efield':False,
  #'do_combined_efield':false,
}

def get_config(args):
  with open(args.c, 'r') as fin:
    config = yaml.safe_load(fin)

  for k,v in defaults.items():
    if k not in config.keys(): config[k] = v

  return config

def get_efield(x, y, z, pos_hists, neg_hists):
  if x > 0: hists = pos_hists
  else: hists = neg_hists

  #print('hists', hists)

  ex = E0 + E0*hists[0].Interpolate(x, y, z)
  ey = E0*hists[1].Interpolate(x, y, z)
  ez = E0*hists[2].Interpolate(x, y, z)
  return sqrt(ex**2 + ey**2 + ez**2)

#def filter_event_for_dedx(e, i, config):
def filter_event_for_dedx(e, i, ediv_min, ediv_max, ediv2_min, ediv2_max,
                          peakT_min, peakT_max, track_len_min, track_len_max):
  if e.trkstartx[i]*e.trkendx[i] > 0: return False #skip non-crossers

  #ediv_min = config['track_ediv_min']
  #ediv_max = config['track_ediv_max']
  #ediv2_min = config['track_ediv2_min']
  #ediv2_max = config['track_ediv2_max']

  #if ((e.peakT_min[i] < config['peakT_min']) or
  #    (e.peakT_max[i] > config['peakT_max']) or
  #    (e.trklen[i] < config['track_len_min']) or
  #    (e.trklen[i] > config['track_len_max']) or
  if ((e.peakT_min[i] < peakT_min) or
      (e.peakT_max[i] > peakT_max) or
      (e.trklen[i] < track_len_min) or
      (e.trklen[i] > track_len_max) or
      ((e.trkendz[i] > ediv_min) and (e.trkendz[i] < ediv_max)) or
      ((e.trkstartz[i] > ediv_min) and (e.trkstartz[i] < ediv_max)) or
      ((e.trkendz[i] > ediv2_min) and (e.trkendz[i] < ediv2_max)) or
      ((e.trkstartz[i] > ediv2_min) and (e.trkstartz[i] < ediv2_max))):
    return False #filter for plane 2
  if e.adjacent_hits[i] != 0 or e.dist_min[i] > 5: return False;
  return True

def filter_plane_for_dedx(e, i, hit_plane, track_index, hit_index, res, dq):
  if len(res) == 0: return False
  if (hit_plane == 2 and (e.lastwire[i] <= 5 or e.lastwire[i] >= 475)):
    return False 
  if (hit_plane != 2 and (e.lastwire[i] <= 5 or e.lastwire[i] >= 795)):
    return False 

  ##Skipping flipped tracks because they have no effect anyways
  nhits = e.ntrkhits[i*3 + hit_plane]
  if (nhits < 0): return False

  hit_y = e.trkhity[track_index + hit_index]
  next_hit_y = e.trkhity[track_index + hit_index + 1]
  resrange = e.trkresrange[track_index + hit_index]
  next_resrange = e.trkresrange[track_index + hit_index + 1]

  if (next_hit_y < hit_y and next_resrange > resrange): return False

  prev_resrange = e.trkresrange[track_index + hit_index - 1]
  prev_hit_y = e.trkhity[track_index + hit_index - 1]
  if (hit_y < prev_hit_y and resrange > prev_resrange): return False
  #if ((e.trkhity[track_index + hit_index + 1] < e.trkhity[track_index + hit_index] and
  #     e.trkresrange[track_index + hit_index + 1] > e.trkresrange[track_index + hit_index]) or
  #    (e.trkhity[track_index + hit_index] < e.trkhity[track_index + hit_index + nhits - 1] and
  #     e.trkresrange[track_index + hit_index ] > e.trkresrange[track_index + hit_index + nhits - 1])):
  #  return False
  #if ((next_hit_y < hit_y and next_resrange > resrange) or
  #    (hit_y < prev_hit_y and resrange > prev_resrange)):
  #  return False 


  first_5_dq = dq[np.where(res < 5)]
  if len(first_5_dq) < 5: return False
  max_res = max(res)
  last_5_dq = dq[np.where(res > max_res - 5)]
  #first_5_dq = [q for r, q in zip(res, dq) if r < 5]
  #last_5_dq = [q for r, q in zip(res, dq) if r > max_res - 5]
  #if len(first_5_dq) < 5: return False


  med1 = np.median(first_5_dq)
  med2 = np.median(last_5_dq)
  #med1 = RT.TMath.Median(len(first_5_dq), array('d', first_5_dq))
  #med2 = RT.TMath.Median(len(last_5_dq), array('d', last_5_dq))
  if med1/med2 < 1.4: return False
  return True

#def filter_event_for_yz_corr(e, i, config):
def filter_event_for_yz_corr(e, i, xmax, ymax, ymin, zmax, zmin):
  ###Check if the start or endpoint are in the FV. Both must be outside to enter sample
  #xmax = config['track_xmax']
  #ymax = config['track_ymax']
  #ymin = config['track_ymin']
  #zmax = config['track_zmax']
  #zmin = config['track_zmin']

  if ((abs(e.trkstartx[i]) < xmax and e.trkstarty[i] > ymin and
       e.trkstarty[i] < ymax and e.trkstartz[i] > zmin and e.trkstartz[i] < zmax) or
      (abs(e.trkendx[i]) < xmax and e.trkendy[i] > ymin and e.trkendy[i] < ymax and
       e.trkendz[i] > zmin and e.trkendz[i] < zmax)):
    return False
  return True



if __name__ == '__main__':
  parser = ap()
  
  parser.add_argument('-o', type=str, required=True)
  parser.add_argument('-c', type=str, required=True)
  parser.add_argument('-i', type=str, required=True)
  parser.add_argument('-e', type=str, required=True)
  parser.add_argument('-n', type=int, default=100)
  parser.add_argument('--nskip', type=int, default=0)
  #parser.add_argument('--ymax', type=float, default=600.)
  #parser.add_argument('--zmax', type=float, default=695.)
  parser.add_argument('--ywidth', type=float, default=5.)
  parser.add_argument('--zwidth', type=float, default=5.)
  args = parser.parse_args()
  
  config = get_config(args) 
  #track_xmax = config['track_xmax']
  #track_ymax = config['track_ymax']
  #track_ymin = config['track_ymin']
  #track_zmax = config['track_zmax']
  #track_zmin = config['track_zmin']

  ymax = config['ymax']
  zmax = config['zmax']
  nbins_y = int(ymax/args.ywidth)
  nbins_z = int(zmax/args.zwidth)
  
  efield_file = RT.TFile.Open(args.e, 'open')
  if config['do_true_efield']:
    pos_hists = [efield_file.Get('True_ElecField_%s'%i) for i in ['X', 'Y', 'Z']]
    neg_hists = pos_hists 
    #print(pos_hists, neg_hists)
  else:
    pos_hists = [efield_file.Get('Reco_ElecField_%s_Pos'%i) for i in ['X', 'Y', 'Z']]
    neg_hists = [efield_file.Get('Reco_ElecField_%s_Neg'%i) for i in ['X', 'Y', 'Z']]

  E0 = 0.4867
  with open(args.i, 'r') as f:
    input_files = [l.strip() for l in f.readlines()]
  print(input_files)
  
  #tree = RT.TChain('michelremoving2/Event')
  level = config['level']
  tree = RT.TChain(f'{level}/Event')
  for i in input_files:
    if i[0] == '#': continue
    if i[:5] == '/pnfs':
      i = i.replace('/pnfs', 'root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr')
    tree.AddFile(i)
  print(tree.GetEntries())
  
  spline_range = [0.70437, 1.27937, 2.37894, 4.72636, 7.5788, 22.0917, 30.4441, 48.2235, 76.1461, 123.567, 170.845, 353.438, 441.476]
  spline_ke = [10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000]
  ke_spline = RT.TSpline3("Cubic Spline", array('d', spline_range), array('d', spline_ke), 13, "b2e2", 0, 0)
  
  fOut = RT.TFile.Open(args.o, 'recreate')
  out_tree = RT.TTree('tree', '')
  corrected_dq_dx = array('d', [0])
  dq_dx = array('d', [0])
  hit_plane_out = array('i', [0])
  hit_x = array('d', [0])
  hit_y = array('d', [0])
  hit_z = array('d', [0])
  efield_out = array('d', [0.])
  ke = array('d', [0.])
  resrange_out = array('d', [0.])
  range_bin = array('i', [0])
  
  dedx_filter = array('i', [0])
  yz_corr_filter = array('i', [0])
  x_corr_filter = array('i', [0])
  Cx = array('i', [0])
  Cyz = array('i', [0])
  
  out_tree.Branch('corrected_dq_dx', corrected_dq_dx, 'corrected_dq_dx/D')
  out_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
  out_tree.Branch('hit_plane', hit_plane_out, 'hit_plane/I')
  out_tree.Branch('hit_x', hit_x, 'hit_x/D')
  out_tree.Branch('hit_y', hit_y, 'hit_y/D')
  out_tree.Branch('hit_z', hit_z, 'hit_z/D')
  out_tree.Branch('efield', efield_out, 'efield/D')
  out_tree.Branch('resrange', resrange_out, 'resrange/D')
  out_tree.Branch('range_bin', range_bin, 'range_bin/I')
  out_tree.Branch('ke', ke, 'ke/D')
  out_tree.Branch('dedx_filter', dedx_filter, 'dedx_filter/I')
  out_tree.Branch('yz_corr_filter', yz_corr_filter, 'yz_corr_filter/I')
  out_tree.Branch('x_corr_filter', x_corr_filter, 'x_corr_filter/I')
  out_tree.Branch('Cx', Cx, 'Cx/I')
  out_tree.Branch('Cyz', Cyz, 'Cyz/I')
  
  #yz_tree = RT.TTree('yz_tree', '')
  yz_trees = [RT.TTree('yz_tree_%i'%i, '') for i in range(3)]
  for yz_tree in yz_trees:
    yz_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
    yz_tree.Branch('hit_plane', hit_plane_out, 'hit_plane/I')
    yz_tree.Branch('hit_x', hit_x, 'hit_x/D')
    yz_tree.Branch('hit_y', hit_y, 'hit_y/D')
    yz_tree.Branch('hit_z', hit_z, 'hit_z/D')
    yz_tree.Branch('efield', efield_out, 'efield/D')
    yz_tree.Branch('resrange', resrange_out, 'resrange/D')
    yz_tree.Branch('range_bin', range_bin, 'range_bin/I')
    yz_tree.Branch('ke', ke, 'ke/D')
  #x_tree = RT.TTree('x_tree', '')
  x_trees = [RT.TTree('x_tree_%i'%i, '') for i in range(3)]
  for x_tree in x_trees:
    x_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
    x_tree.Branch('hit_plane', hit_plane_out, 'hit_plane/I')
    x_tree.Branch('hit_x', hit_x, 'hit_x/D')
    x_tree.Branch('hit_y', hit_y, 'hit_y/D')
    x_tree.Branch('hit_z', hit_z, 'hit_z/D')
    x_tree.Branch('efield', efield_out, 'efield/D')
    x_tree.Branch('resrange', resrange_out, 'resrange/D')
    x_tree.Branch('range_bin', range_bin, 'range_bin/I')
    x_tree.Branch('ke', ke, 'ke/D')
  
  dedx_tree = RT.TTree('dedx_tree', '')
  dedx_tree.Branch('dq_dx', dq_dx, 'dq_dx/D')
  dedx_tree.Branch('hit_plane', hit_plane_out, 'hit_plane/I')
  dedx_tree.Branch('hit_x', hit_x, 'hit_x/D')
  dedx_tree.Branch('hit_y', hit_y, 'hit_y/D')
  dedx_tree.Branch('hit_z', hit_z, 'hit_z/D')
  dedx_tree.Branch('efield', efield_out, 'efield/D')
  dedx_tree.Branch('resrange', resrange_out, 'resrange/D')
  dedx_tree.Branch('range_bin', range_bin, 'range_bin/I')
  dedx_tree.Branch('ke', ke, 'ke/D')
  
  nevent = 0
  xmax = config['track_xmax']
  track_ymax = config['track_ymax']
  track_ymin = config['track_ymin']
  track_zmax = config['track_zmax']
  track_zmin = config['track_zmin']

  ediv_min = config['track_ediv_min']
  ediv_max = config['track_ediv_max']
  ediv2_min = config['track_ediv2_min']
  ediv2_max = config['track_ediv2_max']
  peakT_min = config['peakT_min']
  peakT_max = config['peakT_max']
  track_len_min = config['track_len_min']
  track_len_max = config['track_len_max']

  with cProfile.Profile() as pr:
    for e in tree:
      if nevent < args.nskip:
        nevent += 1
        print('skipping %i'%nevent)
        continue
      if args.n > 0 and nevent > args.n: break
      print(nevent, end='\r')
      nevent += 1
      for i in range(0, e.cross_trks):
        track_index = i*9000
        #dedx_event_filter = filter_event_for_dedx(e, i, config)
        dedx_event_filter = filter_event_for_dedx(e, i, ediv_min, ediv_max,
                                                  ediv2_min, ediv2_max,
                                                  peakT_min, peakT_max,
                                                  track_len_min, track_len_max)
        #yz_corr_filter[0] = filter_event_for_yz_corr(e, i, config)
        yz_corr_filter[0] = filter_event_for_yz_corr(e, i, xmax,
                                                     track_ymax, track_ymin,
                                                     track_zmax, track_zmin)
        x_cross_filter = e.trkstartx[i]*e.trkendx[i] < 0
    
        if (not dedx_event_filter) and (not yz_corr_filter[0]): continue
    
        testneg = (e.trkstartx[i]<-xmax or e.trkendx[i]<-xmax)
        testpos = (e.trkstartx[i]>xmax or e.trkendx[i]>xmax)
    
        theta_xz_deg = 180./RT.TMath.Pi()*e.trackthetaxz[i]
        theta_yz_deg = 180./RT.TMath.Pi()*e.trackthetayz[i]
        for hit_plane in range(0, 3):
          hit_index = hit_plane*3000
          hit_plane_out[0] = hit_plane
          if (hit_plane == 2 and
              ((abs(theta_xz_deg)>config['plane2_theta_xz_min'] and
                abs(theta_xz_deg)<config['plane2_theta_xz_max']) or
               (abs(theta_yz_deg)>config['plane2_theta_yz_min'] and
                abs(theta_yz_deg)<config['plane2_theta_yz_max']))): continue
    
          #res = []
          #dq = []
          nhits = e.ntrkhits[i*3 + hit_plane]
          start_index = track_index + hit_index
          res = np.array(e.trkresrange)[start_index:start_index + nhits]
          dq = np.array(e.trkdqdx)[start_index:start_index + nhits]
          #for k in range(track_index + hit_index, track_index + hit_index + nhits):
          #  res.append(e.trkresrange[k])
          #  dq.append(e.trkdqdx[k])
          #dq = np.array(dq)
          #res = np.array(res)
          dedx_plane_filter = filter_plane_for_dedx(e, i, hit_plane, track_index, hit_index, res, dq)
          hit_limit = min(nhits, 3000)
          for j in range(0, hit_limit):
            hit_x[0] = e.trkhitx[track_index + hit_index + j]
            hit_y[0] = e.trkhity[track_index + hit_index + j]
            hit_z[0] = e.trkhitz[track_index + hit_index + j]
    
            if (hit_y[0] < 0 or hit_y[0] > ymax or hit_z[0] < 0 or
                hit_z[0] > config['zmax']):
              continue
    
            ##negative X
            if (hit_x[0] > -360 and hit_x[0] < 0.):
              hit_pos = False
              if hit_plane == 1 and abs(theta_xz_deg) < 140: continue
              elif hit_plane == 0 and abs(theta_xz_deg) > 40: continue
            elif (hit_x[0] > 0. and hit_x[0] < 360.):
              hit_pos = True
              if hit_plane == 1 and abs(theta_xz_deg) > 40: continue
              elif hit_plane == 0 and abs(theta_xz_deg) < 140: continue
            else: continue
    
            ke[0] = ke_spline.Eval(e.trkresrange[track_index + hit_index + j]);
            resrange_out[0] = e.trkresrange[track_index + hit_index + j]
            range_bin[0] = int(resrange_out[0]/5.) 
            hit_trkpitch = e.trkpitch[track_index + hit_index + j]
    
            ke_pitch_filter = (ke[0] > 250. and ke[0] < 450. and
                               hit_trkpitch > 0.5 and hit_trkpitch < 0.8)
    
            dedx_filter[0] = dedx_event_filter and dedx_plane_filter and ke_pitch_filter
    
            if (not dedx_filter[0]) and (not yz_corr_filter[0]): continue
    
            dq_dx[0] = e.trkdqdx[track_index + hit_index + j]
            if dedx_filter[0]:
              efield_out[0] = get_efield(hit_x[0], hit_y[0], hit_z[0],
                                         pos_hists, neg_hists)
            else: efield_out[0] = -999.
            
            if (j < nhits-1) and (j > 0):
              if testpos and hit_pos:
                good_x_hit = (e.trkhitx[track_index + hit_index + j+1] > 0 and e.trkhitx[track_index + hit_index + j-1] > 0)
              elif testneg and not hit_pos:
                good_x_hit = (e.trkhitx[track_index + hit_index + j+1] < 0 and e.trkhitx[track_index + hit_index + j-1] < 0)
              else: good_x_hit = False
            else: good_x_hit = False
            x_corr_filter[0] = good_x_hit and x_cross_filter
    
            if (not dedx_filter[0]) and (not yz_corr_filter[0]) and (not x_corr_filter[0]): continue
    
            if dedx_filter[0]: dedx_tree.Fill()
            if yz_corr_filter[0]:
              #yz_tree.Fill()
              yz_trees[hit_plane].Fill()
              y_bin = int(hit_y[0]/args.ywidth)
              z_bin = int(hit_z[0]/args.zwidth)
              #if y_bin < nbins_y and z_bin < nbins_z:
              #  if hit_x[0] < 0.:
              #    dqdx_neg[hit_plane][z_bin][y_bin].push_back(dq_dx[0])
              #    all_dqdx_neg[hit_plane].append(dq_dx[0])
              #  else:
              #    dqdx_pos[hit_plane][z_bin][y_bin].push_back(dq_dx[0])
              #    all_dqdx_pos[hit_plane].append(dq_dx[0])
            #if x_corr_filter[0]: x_tree.Fill()
            if x_corr_filter[0]: x_trees[hit_plane].Fill()
            out_tree.Fill()
  
  stats = pstats.Stats(pr)
  stats.sort_stats(pstats.SortKey.TIME)
  stats.print_stats()
  stats.dump_stats(filename='profile.prof')

  fOut.cd()
  out_tree.Write()
  for t in yz_trees: t.Write()
  for t in x_trees: t.Write()
  #yz_tree.Write()
  #x_tree.Write()
  dedx_tree.Write()
  
  '''
  for i in range(len(dqdx_pos)):
    for j in range(len(dqdx_pos[i])):
      for k in range(len(dqdx_pos[i][j])):
        fOut.WriteObject(dqdx_pos[i][j][k], 'dqdx_pos_%i_%i_%i'%(i, j, k))
        fOut.WriteObject(dqdx_neg[i][j][k], 'dqdx_neg_%i_%i_%i'%(i, j, k))
  
  median_dqdx_pos = RT.TVectorD(3*nbins_y*nbins_z)
  median_dqdx_neg = RT.TVectorD(3*nbins_y*nbins_z)
  
  for i in range(3):
    global_median_dqdx_pos = RT.TVectorD(1)
    global_median_dqdx_neg = RT.TVectorD(1)
    global_median_dqdx_pos[0] = median(all_dqdx_pos[i])
    global_median_dqdx_neg[0] = median(all_dqdx_neg[i])
  
    global_median_dqdx_pos.Write('global_median_dqdx_pos_%i'%i)
    global_median_dqdx_neg.Write('global_median_dqdx_neg_%i'%i)
    for j in range(len(dqdx_pos[i])):
      for k in range(len(dqdx_pos[i][j])):
        median_dqdx_pos[i*nbins_y*nbins_z + j*nbins_y + k] = median(dqdx_pos[i][j][k]) if len(dqdx_pos) >= 5 else -999.
  median_dqdx_pos.Write('median_dqdx_pos')
  median_dqdx_neg.Write('median_dqdx_neg')
  '''
  
  fOut.Close()
