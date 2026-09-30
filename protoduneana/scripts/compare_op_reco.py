from utils import handle_header, provide_list
from classes import ophitv, opwfv, opflashv, opdetwfv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import numpy as np
import h5py as h5



class Geom:
    def __init__(
        self,
        self_triggered_size=1024,
        streaming_size=343808,
        streaming_apas=[0],
        apas={
            0:(120,159),
            1:(80,119),
            2:(40,79),
            3:(0, 39)
        }
    ):
        self.apas=apas
        self.self_triggered_size=self_triggered_size
        self.streaming_size=streaming_size
        self.streaming_apas=streaming_apas

if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--digits', type=str, default=None)
    parser.add_argument('--deco', type=str, default=None)
    parser.add_argument('--hits', type=str, default=None)
    parser.add_argument('--flashes', type=str, default=None)
    parser.add_argument('-o', type=str, default='dump_op_reco.h5')
    # parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()

    # if args.plot:
    #     t = np.load(args.i)
    #     plot(t)
    #     exit()

    handle_header()
    provide_list(ophitv, opwfv, opflashv, opdetwfv)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    with h5.File(args.o, 'w') as fout:
        if args.digits is not None:
            get_opdetwf = ev.getValidHandle[opdetwfv]    
            opdetwfs = get_opdetwf(InputTag(args.digits))
            output_wfs = {}
            output_tss = {}
            odwf_group = fout.create_group('opdetwf')
            for n, opdetwf in enumerate(opdetwfs.product()):
                adcs = opdetwf.Waveform()
                chan = opdetwf.ChannelNumber()
                if chan not in output_wfs:
                    output_wfs[chan] = []
                    output_tss[chan] = []
                output_wfs[chan].append(
                    np.frombuffer(
                                    adcs.data(), dtype=np.int16, count=adcs.size()
                                ).copy()
                )
                output_tss[chan].append(
                    opdetwf.TimeStamp()
                )

                    
                print(f'{n}/{len(opdetwfs.product())}', end='\r')
            print()

            for chan, wfs in output_wfs.items():
                chan_group = odwf_group.create_group(str(chan))
                chan_group.create_dataset('timestamps', data=output_tss[chan])
                chan_group.create_dataset('waveforms', data=output_wfs[chan])
                # print('Wrote', chan, len(output_wfs[chan]))
            del output_wfs
            del output_tss


        if args.deco is not None:
            
            get_opwf = ev.getValidHandle[opwfv]
            opwfs = get_opwf(InputTag(args.deco))
            output_wfs = {}
            output_tss = {}

            wf_group = fout.create_group('opwf')
            for n, opwf in enumerate(opwfs.product()):
                vals = opwf.Signal()
                chan = opwf.Channel()
                if chan not in output_wfs:
                    output_wfs[chan] = []
                    output_tss[chan] = []
                output_wfs[chan].append(
                    np.frombuffer(
                                    vals.data(), dtype=np.float32, count=vals.size()
                                ).copy()
                )
                output_tss[chan].append(
                    opwf.TimeStamp()
                )

                    
                print(f'{n}/{len(opwfs.product())}', end='\r')
            print()

            for chan, wfs in output_wfs.items():
                chan_group = wf_group.create_group(str(chan))
                chan_group.create_dataset('timestamps', data=output_tss[chan])
                chan_group.create_dataset('waveforms', data=output_wfs[chan])
                # print('Wrote', chan, len(output_wfs[chan]))
            del output_wfs, output_tss

        if args.hits is not None:
            
            get_hits = ev.getValidHandle[ophitv]
            hits = get_hits(InputTag(args.hits))
            output = {}

            hit_group = fout.create_group('ophit')
            for n, hit in enumerate(hits.product()):
                chan = hit.OpChannel()
                if chan not in output:
                    output[chan] = []

                output[chan].append(
                    [hit.PeakTime(), hit.StartTime(), hit.Width(), hit.Amplitude(), hit.Area(), hit.PE()]
                )
                    
                print(f'{n}/{len(hits.product())}', end='\r')
            print()
            for chan, hits in output.items():
                chan_group = hit_group.create_group(str(chan))
                chan_group.create_dataset('hits', data=hits)
            del output
