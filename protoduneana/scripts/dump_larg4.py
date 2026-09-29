from utils import read_header, provide_list
from classes import mcpartv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import h5py as h5
import numpy as np

"""
To read one particle (i):
    g = h5py.File('dump_larg4.h5')['largeant']
    off = g['offsets'][:]
    pos_i = g['positions'][off[i]:off[i+1]]      # (npts, 4) — x, y, z, t
    mom_i = g['momenta'][off[i]:off[i+1]]        # (npts, 4) — px, py, pz, E
    proc  = g['processes'][i].decode()
"""


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--tags', type=str, nargs='+', default=['largeant'])
    parser.add_argument('-o', type=str, default='dump_larg4.h5')
    parser.add_argument('--compression', type=str, default='gzip',
                        choices=['gzip', 'lzf', 'none'])
    parser.add_argument('--complevel', type=int, default=4)
    args = parser.parse_args()

    #HDF5 keeps variable-length data in the global heap, where the filters
    #can't reach it, so the trajectories are concatenated into flat (N, 4)
    #datasets instead. 'offsets' gives each particle its slice of them:
    #  pos_i = positions[offsets[i]:offsets[i+1]]
    comp = {}
    if args.compression != 'none':
        comp['compression'] = args.compression
        comp['shuffle'] = True
        if args.compression == 'gzip':
            comp['compression_opts'] = args.complevel

    read_header('gallery/ValidHandle.h')
    provide_list(mcpartv)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_mcpart = ev.getValidHandle[mcpartv]

    with h5.File(args.o, 'w') as f:
        for tag in args.tags:
            tag_group = f.create_group(tag)
            mcparts = get_mcpart(InputTag(tag))

            pdgs = []
            trackids = []
            mothers = []
            processes = []
            ntrajpoints = []
            positions = []
            momenta = []
            for part in mcparts.product():
                pdgs.append(part.PdgCode())
                trackids.append(part.TrackId())
                mothers.append(part.Mother())
                processes.append(str(part.Process()))

                npts = part.NumberTrajectoryPoints()
                ntrajpoints.append(npts)

                for i in range(npts):
                    pos = part.Position(i)
                    mom = part.Momentum(i)
                    positions.append([pos.X(), pos.Y(), pos.Z(), pos.T()])
                    momenta.append([mom.Px(), mom.Py(), mom.Pz(), mom.E()])

            ntrajpoints = np.asarray(ntrajpoints, dtype='int32')
            offsets = np.zeros(len(ntrajpoints) + 1, dtype='int64')
            np.cumsum(ntrajpoints, out=offsets[1:])

            #Process names are written as fixed-width bytes rather than
            #variable-length strings so they land in the compressed chunks too
            maxproc = max([len(p) for p in processes] + [1])

            def make(name, data):
                #A chunk can't be bigger than an empty dataset, and there's
                #nothing to compress in one either
                if len(data) == 0:
                    return tag_group.create_dataset(name, data=data)
                chunk = (min(4096, len(data)),) + data.shape[1:]
                return tag_group.create_dataset(name, data=data, chunks=chunk,
                                                **comp)

            make('PDGs', np.asarray(pdgs, dtype='int32'))
            make('trackids', np.asarray(trackids, dtype='int32'))
            make('mothers', np.asarray(mothers, dtype='int32'))
            make('processes', np.asarray(processes, dtype='S%d' % maxproc))
            make('ntrajpoints', ntrajpoints)
            make('offsets', offsets)
            make('positions', np.asarray(positions, dtype='float64').reshape(-1, 4))
            make('momenta', np.asarray(momenta, dtype='float64').reshape(-1, 4))
