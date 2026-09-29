from utils import handle_header, provide_list
from classes import wirev
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import numpy as np


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--tag', type=str, default='wclsdatahd:gauss')
    parser.add_argument('-o', type=str, default='dump_wires.npy')
    args = parser.parse_args()


    handle_header()
    provide_list(wirev)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_wires = ev.getValidHandle[wirev]

    wires = get_wires(InputTag(args.tag))

    output = {}
    for n, wire in enumerate(wires.product()):
        print(f'{n}/{len(wires.product())}', end='\r')
        signals = wire.Signal()
        output[wire.Channel()] = np.frombuffer(
            signals.data(), dtype=np.float32, count=signals.size()
        ).copy()
    print()
    output = [output[k] for k in sorted(output)]
    np.save(args.o, output)