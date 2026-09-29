from utils import handle_header, provide_list
from classes import trackv, showerv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import numpy as np


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--shower-tag', type=str, default='pandoraShower')
    parser.add_argument('--track-tag', type=str, default='pandoraTrack')
    parser.add_argument('-o', type=str, default='dump_reco.npz')
    args = parser.parse_args()


    handle_header()
    provide_list(trackv, showerv)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_showers = ev.getValidHandle[showerv]
    get_tracks = ev.getValidHandle[trackv]

    tracks = get_tracks(InputTag(args.track_tag))
    showers = get_showers(InputTag(args.shower_tag))

    track_pts = []
    for i, track in enumerate(tracks.product()):
        # print(f'{n}/{len(beams.product())}', end='\r')
        for j in range(track.NumberTrajectoryPoints()):
            pt = track.TrajectoryPoint(j).position
            track_pts.append([pt.X(), pt.Y(), pt.Z()])

    shower_infos = []
    for i, shower in enumerate(showers.product()):

        shower_infos.append([
            shower.ShowerStart().X(), shower.ShowerStart().Y(), shower.ShowerStart().Z(),
            shower.Direction().X(), shower.Direction().Y(), shower.Direction().Z(),
             shower.Length(), shower.OpenAngle(),
        ])
    np.savez(args.o, tracks=track_pts, showers=shower_infos)