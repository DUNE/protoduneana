from utils import handle_header, provide_list
from classes import tsv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import numpy as np


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--tag', type=str, default='timingrawdecoder:daq')
    parser.add_argument('-o', type=str, default='')
    args = parser.parse_args()


    handle_header()
    provide_list(tsv)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_tss = ev.getValidHandle[tsv]

    tss = get_tss(InputTag(args.tag))

    output = {}
    for n, ts in enumerate(tss.product()):
        print(ts.GetTimeStamp()*16e-9, ts.GetTimeStamp_Low(), ts.GetTimeStamp_High())
    print()
