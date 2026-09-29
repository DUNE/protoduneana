from utils import handle_header, provide_list
from classes import beamv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import numpy as np


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--tag', type=str, default='generator')
    parser.add_argument('-o', type=str, default='dump_beam.npy')
    args = parser.parse_args()


    handle_header()
    provide_list(beamv)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_beams = ev.getValidHandle[beamv]

    beams = get_beams(InputTag(args.tag))

    output = {}
    for n, beam in enumerate(beams.product()):
        # print(f'{n}/{len(beams.product())}', end='\r')
        track = beam.GetBeamTrack(0)
        for p in range(track.NumberTrajectoryPoints()):
            pos = track.TrajectoryPoint(p).position
            print('\t', pos.X(), pos.Y(), pos.Z())
        print(beam.GetRecoBeamMomenta(),)
        # output[beam.Channel()] = np.frombuffer(
        #     adcs.data(), dtype=np.int16, count=adcs.size()
        # ).copy()
    print()
    # output = [output[k] for k in sorted(output)]
    # np.save(args.o, output)
