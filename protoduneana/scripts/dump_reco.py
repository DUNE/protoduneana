from utils import handle_header, provide_list
from classes import trackv, showerv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import numpy as np
import cppyy

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


    track_pts = []
    try:
        tracks = get_tracks(InputTag(args.track_tag))
        for i, track in enumerate(tracks.product()):
            # print(f'{n}/{len(beams.product())}', end='\r')
            for j in range(track.NumberTrajectoryPoints()):
                pt = track.TrajectoryPoint(j).position
                track_pts.append([pt.X(), pt.Y(), pt.Z()])
    except cppyy.gbl.std.exception:
        print('Probably could not get the tracks from event')
        

    shower_infos = []
    try:
        showers = get_showers(InputTag(args.shower_tag))
        for i, shower in enumerate(showers.product()):

            shower_infos.append([
                shower.ShowerStart().X(), shower.ShowerStart().Y(), shower.ShowerStart().Z(),
                shower.Direction().X(), shower.Direction().Y(), shower.Direction().Z(),
                shower.Length(), shower.OpenAngle(),
            ] + [e for e in shower.Energy()])
            print(shower_infos[-1][8:])
    except cppyy.gbl.std.exception:
            print('Probably could not get the showers from event')
    np.savez(args.o, tracks=track_pts, showers=shower_infos)