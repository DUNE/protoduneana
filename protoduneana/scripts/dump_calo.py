from utils import handle_header, provide_list
from classes import calov
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import numpy as np


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--tag', type=str, default='pandoraGnocchiCalo')
    parser.add_argument('-o', type=str, default='dump_calo.npy')
    args = parser.parse_args()


    handle_header()
    provide_list(calov)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_calos = ev.getValidHandle[calov]

    calos = get_calos(InputTag(args.tag))

    print(len(calos.product()), 'calos')

    calo_out = []
    for i, calo in enumerate(calos.product()):
        xyz = calo.XYZ()
        dedx = calo.dEdx()
        rr = calo.ResidualRange()
        for j in range(len(xyz)):
            xyzj = xyz[j]
            calo_out.append([
                xyzj.X(), xyzj.Y(), xyzj.Z(), dedx[j], rr[j],
                calo.PlaneID().Plane, calo.PlaneID().TPC,
            ])

    # track_pts = []
    # for i, track in enumerate(tracks.product()):


    np.save(args.o, calo_out)