import numpy as np, h5py as h5, ROOT as RT
from argparse import ArgumentParser as ap
from scipy.spatial import Delaunay

##NEED TO UPDATE SO EACH EVENT IS FLAT IN THE OUTPUT -- ALSO OPEN AND REWRITE?

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
  results = np.zeros(len(e.x)*9)
  results[0:len(e.x)*9:9] = e.x
  results[1:len(e.x)*9:9] = e.y
  results[2:len(e.x)*9:9] = e.z

  results[3:len(e.x)*9:9] = e.charge0
  results[4:len(e.x)*9:9] = e.charge1
  results[5:len(e.x)*9:9] = e.charge2

  results[6:len(e.x)*9:9] = e.rms0
  results[7:len(e.x)*9:9] = e.rms1
  results[8:len(e.x)*9:9] = e.rms2
  return results

def get_unique_indices():
  return (np.array([0, 0, 0, 1, 1, 2]), np.array([0, 1, 2, 1, 2, 2]))

def get_delaunay_edges(e, nf):
  points = np.zeros((len(e.x), 3))
  points[:, 0] = nf[0:len(e.x)*9:9]
  points[:, 1] = nf[1:len(e.x)*9:9]
  points[:, 2] = nf[2:len(e.x)*9:9]

  tri = Delaunay(points)

  tri_points = np.zeros((len(tri.simplices)*12, 2))
  start = 0
  for v in tri.simplices:
    pairs = np.array([[i,j] for i in v for j in v if i != j])
    tri_points[start:start+12, :] = pairs
    start += 12

  return np.unique(tri_points.astype(int), axis=0).T
  
def make_edge_i(n, knn=5):
  print('Making', n)
  return np.repeat(np.arange(n), knn, axis=0)

  
def get_formed_edges(e, knn=5, force_bidirection=False):
  edge_indices = [make_edge_i(len(e.edge_i), knn=knn),
                  np.array(e.edge_j)[:, :knn].flatten()]

  if force_bidirection:
    edge_indices = np.array(edge_indices)
    #print(edge_indices.shape)
    all_edges = np.zeros((2*len(edge_indices[0]), 2), dtype=int)
    all_edges[0::2] = edge_indices.T[:, ::-1]
    all_edges[1::2] = edge_indices.T[:]

    #print(all_edges.shape)
    #print(all_edges[:10])
    edge_indices = np.unique(all_edges, axis=0)
    #print(edge_indices.shape)
    edge_indices = edge_indices.T

  return edge_indices

#def get_formed_edges(e, knn=5):
#  edge_indices = [np.array( 
#  #edge_indices = [np.array(e.edge_i)[:, :knn].flatten(),
#  #                np.array(e.edge_j)[:, :knn].flatten()]
#
#  return edge_indices


def make_edges(nf, i, j):
  nf_reshaped = nf.reshape(-1, 9)
  nfi = nf_reshaped[i, :]
  nfj = nf_reshaped[j, :]
  print(nfi, nfj)

  d = nfi - nfj
  print(d.shape)

  ave = (nfi + nfj)/2.

  #ave_c = (nfi[:, 3:6] + nfj[:, 3:6]) / 2.
  #print(ave_c.shape)

  nfeatures = 14 # + 2 for edge indices
  nedges = ave.shape[0]

  edge_features = np.zeros((nedges, nfeatures))
  edge_features[:, :6] = d[:, :6] #diff in position (3), charge (3)
  edge_features[:, 6:12] = ave[:, :6] #ave in position (3), charge (3)
  edge_features[:, -2] = i
  edge_features[:, -1] = j

  return edge_features.flatten()

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
  ##Setting up dataset. Each entry will coincide with an event.
  ##Each entry will be nfeatures*nobs (i.e. nodes/edges)
  dt = h5.vlen_dtype(np.dtype('float32'))
  dset = f.create_dataset(dname,
                          (init_entries,),
                          dtype=dt,
                          maxshape=(None,))
  return dset

def get_beam_fraction(e):
  return e.reco_beam_fraction
def get_truth(e): 
  results = np.zeros(6)
  if e.true_pdg == -13: ##Is muon
    results[3] = 1.
  elif e.true_pdg == 211:
    if e.true_end_z < -.49375: ##Upstream int
      results[5] = 1.
    elif e.true_end_process == 'pi+Inelastic':
      if e.true_n_piplus == 0 and e.true_n_piminus == 0 and e.true_n_pi0 == 0: ##Abs
        results[0] = 1.
      elif e.true_n_piplus == 0 and e.true_n_piminus == 0 and e.true_n_pi0 == 1: ##Ch. Ex
        results[1] = 1.
      else: ##Other int
        results[2] = 1.
    else: ##Other/decay?
      results[4] = 1.
  else:
    print('WARNING') 

  return results

def get_truth_pdg(e): 
  results = np.zeros(4)
  if e.true_pdg == -11:
    results[0] = 1.
  elif e.true_pdg == 211:
    results[1] = 1.
  elif e.true_pdg == 2212:
    results[2] = 1.
  elif e.true_pdg == -13: ##Is muon
    results[3] = 1.
  else:
    print(f'WARNING: pdg is {e.true_pdg}')

  return results

if __name__ == '__main__':
  parser = ap()
  parser.add_argument('-i', required=True, type=str)
  parser.add_argument('-o', default='convert.hdf5', type=str)
  parser.add_argument('-n', type=int, default=1, help='Number of files')
  parser.add_argument('--nskip', type=int, default=0, help='Beginning file')
  parser.add_argument('--max', type=int, default=-1)
  parser.add_argument('--input_root', action='store_true')
  parser.add_argument('--delaunay', action='store_true')
  parser.add_argument('--bidirection', action='store_true')
  args = parser.parse_args()

  if args.input_root:
    filelist = [args.i]
  else:
    with open(args.i, 'r') as f:
      files = [
        l.strip('\n') for l in f.readlines()
      ]
    filelist = files[args.nskip:args.nskip+args.n]
  with h5.File(args.o, 'w') as h5_file:
    total_entries = 0
    #for iFile, filename in enumerate(files[:args.max]):
    for iFile, filename in enumerate(filelist):
      with RootFile(filename) as f:
        print(f)
        t = get_tree(f)
        i = total_entries
        total_entries += t.GetEntries()
        if iFile == 0:
          #TODO -- remove nfeatures
          output_nodes = setup_dataset(h5_file, t.GetEntries())
          output_edges = setup_dataset(h5_file, t.GetEntries(),
                                       dname='edge_features', nfeatures=14)
          output_edge_indices = setup_dataset(
              h5_file, t.GetEntries(),
              dname='edge_indices', nfeatures=2)
          output_truth = setup_normal_dataset(
              h5_file, t.GetEntries(),
              dname='truth', nfeatures=6)
          output_truth_pdg = setup_normal_dataset(
              h5_file, t.GetEntries(),
              dname='truth_pdg', nfeatures=4)
          output_beam_fraction = setup_normal_dataset(
              h5_file, t.GetEntries(),
              dname='beam_fraction', nfeatures=1)

        else:
          print('Resizing')
          output_nodes.resize(total_entries, axis=0)
          output_edges.resize(total_entries, axis=0)
          output_edge_indices.resize(total_entries, axis=0)
          output_truth.resize(total_entries, axis=0)
          output_truth_pdg.resize(total_entries, axis=0)
          output_beam_fraction.resize(total_entries, axis=0)
          print('Done')
        a = 0
        for e in t:
          print(i)
          if a >= args.max and args.max > 0: break
          nf = get_node_features(e)
          ei, ej = get_delaunay_edges(e, nf) if args.delaunay else  get_formed_edges(e, force_bidirection=args.bidirection)
          ef = make_edges(nf, ei, ej)

          #edge_indices = np.zeros((2, len(ei)))
          #edge_indices[0, :] = ei
          #edge_indices[1, :] = ej

          #print('nf shape', nf.shape)
          output_nodes[i] = nf
          output_edges[i] = ef#np.array(ef).T
          #output_edge_indices[i] = edge_indices#np.array(ei)
          output_truth[i] = get_truth(e)
          output_truth_pdg[i] = get_truth_pdg(e)
          output_beam_fraction[i] = get_beam_fraction(e)

          i += 1
          a += 1
