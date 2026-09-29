from utils import handle_header, provide_list
from classes import digitv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import numpy as np


def plot(frame):
    import matplotlib.pyplot as plt
    # plt.ion()
    plt.imshow((frame.T - np.median(frame, axis=1)).T, aspect='auto')
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--tag', type=str, default='tpcrawdecoder:daq')
    parser.add_argument('-o', type=str, default='dump_digits.npy')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()

    if args.plot:
        t = np.load(args.i)
        plot(t)
        exit()

    handle_header()
    provide_list(digitv)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_digits = ev.getValidHandle[digitv]

    digits = get_digits(InputTag(args.tag))

    output = {}
    for n, digit in enumerate(digits.product()):
        print(f'{n}/{len(digits.product())}', end='\r')
        adcs = digit.ADCs()
        output[digit.Channel()] = np.frombuffer(
            adcs.data(), dtype=np.int16, count=adcs.size()
        ).copy()
    print()
    output = [output[k] for k in sorted(output)]
    np.save(args.o, output)