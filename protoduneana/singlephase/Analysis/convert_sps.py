import numpy as np, h5py as h5, ROOT as RT
from argparse import ArgumentParser as ap

class RootFile():
  def __init__(self, filename):
    self.filename = filename.replace(
      '/pnfs', 'root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr')
  def __enter__(self):
    self.file = RT.TFile.Open(self.filename)
    return self.file
  def __exit__(self, *args):
    self.file.Close()

def get_tree(f):
  return f.Get('savehits/space_points')

def get_node_features(e):
  return np.array([
     e.x, e.y, e.z,
     e.charge0, e.charge1, e.charge2,
     e.rms0, e.rms1, e.rms2
  ])

def get_unique_indices():
  return (np.array([0, 0, 0, 1, 1, 2]), np.array([0, 1, 2, 1, 2, 2]))

def get_formed_edges(e, knn=5):
  edge_indices = [np.array(e.edge_i)[:, :knn].flatten(),
                  np.array(e.edge_j)[:, :knn].flatten()]

  return edge_indices

def make_edges(nf, i, j):
  nfi = nf[:, i]
  nfj = nf[:, j]

  d = nfi - nfj
  outer_indices = [(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)]
  d_out = np.array([
    d[oi]*d[oj] for oi, oj in outer_indices
  ])

  dc_out = np.array([
    d[3+oi]*d[3+oj] for oi, oj in outer_indices
  ])
  #d_out0 = d[0]*d[0]
  #d_out1 = d[0]*d[1]
  #d_out2 = d[0]*d[2]
  #d_out3 = d[1]*d[1]
  #d_out4 = d[1]*d[2]
  #d_out5 = d[2]*d[2]

  #dc_out0 = d[3]*d[3]
  #dc_out1 = d[3]*d[4]
  #dc_out2 = d[3]*d[5]
  #dc_out3 = d[4]*d[4]
  #dc_out4 = d[4]*d[5]
  #dc_out5 = d[5]*d[5]

  ave_c = (nfi[3:6, :] + nfj[3:6, :]) / 2.

  edge_features = np.zeros((21, len(i)))
  edge_features[:3, :] = d[:3]
  edge_features[3:9 :] = d_out
  edge_features[9:12, :] = d[3:6]
  edge_features[12:18, :] = dc_out
  edge_features[18:21, :] = ave_c

  return edge_features

def get_edges(nf, knn=5):
  edge_indices = [[], []]
  edge_features = []
  unique = get_unique_indices()
  for i in range(nf.shape[1]):
    ni = nf[:, i]

    dists = []
    for j in range(nf.shape[1]):
      if j == i: continue
      nj = nf[:, j]

      dij = ni - nj
      dists.append([i, j, dij, np.sqrt(dij[:3].dot(dij[:3]))])
      #dists.append([j, np.sqrt(dij.dot(dij))])

    nn = sorted(dists, key=lambda x: x[3])[:knn]
    #print(nn)
    for n in nn:
      edge_indices[0].append(n[0])
      edge_indices[1].append(n[1])
      #edge_features.append(n[2])

      ds = n[2][:3]
      outer_space = np.outer(ds, ds)[unique]

      dc = n[2][3:6]
      outer_charge = np.outer(dc, dc)[unique]

      ci = nf[3:6, n[0]]
      cj = nf[3:6, n[1]]

      ave_charge = (ci + cj)/2.

      #Output len is 21
      edge_features.append(
        [*ds, *outer_space, *dc, *outer_charge, *ave_charge] 
      )
  return edge_indices, edge_features

#def get_edge_features(nf1, nf2):
#  np.array([
#    nf1[0] - nf2[0], nf1[0] - nf2[0],
#  ])

def setup_normal_dataset(f, init_entries, dname, nfeatures):
  dt = np.dtype('float32')
  dset = f.create_dataset(dname,
                          (init_entries, nfeatures),
                          dtype=dt,
                          maxshape=(None, nfeatures)
                         )
  return dset

def setup_dataset(f, init_entries, dname='node_features', nfeatures=9):
  dt = h5.vlen_dtype(np.dtype('float32'))
  dset = f.create_dataset(dname,
                          (init_entries, nfeatures, ),
                          dtype=dt,
                          maxshape=(None, nfeatures, ))
  return dset

def get_beam_fraction(e):
  return e.reco_beam_fraction
def get_truth(e): 
  results = np.zeros(6)
  if e.true_pdg == -13:
    results[3] = 1.
  elif e.true_pdg == 211:
    if e.true_end_z < -.49375:
      results[5] = 1.
    elif e.true_end_process == 'pi+Inelastic':
      if e.true_n_piplus == 0 and e.true_n_piminus == 0 and e.true_n_pi0 == 0:
        results[0] = 1.
      elif e.true_n_piplus == 0 and e.true_n_piminus == 0 and e.true_n_pi0 == 1: 
        results[1] = 1.
      else:
        results[2] = 1.
    else: 
      results[4] = 1.
  else:
    print('WARNING') 

  return results

if __name__ == '__main__':
  parser = ap()
  parser.add_argument('-i', required=True, type=str)
  parser.add_argument('-o', default='convert.hdf5', type=str)
  parser.add_argument('--max', type=int, default=None)
  args = parser.parse_args()

  with open(args.i, 'r') as f:
    files = [
      l.strip('\n') for l in f.readlines()
    ]
  print(files)

  with h5.File(args.o, 'w') as h5_file:
    total_entries = 0
    for iFile, filename in enumerate(files[:args.max]):
      with RootFile(filename) as f:
        print(f)
        t = get_tree(f)
        i = total_entries
        total_entries += t.GetEntries()
        if iFile == 0:
          output_nodes = setup_dataset(h5_file, t.GetEntries())
          output_edges = setup_dataset(h5_file, t.GetEntries(),
                                       dname='edge_features', nfeatures=21)
          output_edge_indices = setup_dataset(
              h5_file, t.GetEntries(),
              dname='edge_indices', nfeatures=2)
          output_truth = setup_normal_dataset(
              h5_file, t.GetEntries(),
              dname='truth', nfeatures=6)
          output_beam_fraction = setup_normal_dataset(
              h5_file, t.GetEntries(),
              dname='beam_fraction', nfeatures=1)

        else:
          print('Resizing')
          output_nodes.resize(total_entries, 0)
          output_edges.resize(total_entries, 0)
          output_edge_indices.resize(total_entries, 0)
          output_truth.resize(total_entries, 0)
          output_beam_fraction.resize(total_entries, 0)
          print('Done')
        for e in t:
          print(i)
          nf = get_node_features(e)
          ei, ej = get_formed_edges(e)
          ef = make_edges(nf, ei, ej)
          #ei, ef = get_edges(nf)

          edge_indices = np.zeros((2, len(ei)))
          edge_indices[0, :] = ei
          edge_indices[1, :] = ej

          output_nodes[i] = nf
          output_edges[i] = ef#np.array(ef).T
          output_edge_indices[i] = edge_indices#np.array(ei)
          output_truth[i] = get_truth(e)
          output_beam_fraction[i] = get_beam_fraction(e)

          i += 1
