from utils import read_header, provide_list
from classes import mctruthv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import h5py as h5

# def create_datasets(group):
    # group.create_dataset()

if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--tags', type=str, nargs='+')
    parser.add_argument('-o', type=str, default='dump_mctruth.h5')
    args = parser.parse_args()


    header_handle()
    provide_list(mctruthv)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_mctruth = ev.getValidHandle[mctruthv]

    with h5.File(args.o, 'w') as f:

        for tag in args.tags:
            tag_group = f.create_group(tag)
            mctruths = get_mctruth(InputTag(tag))
            mct = mctruths.product()[0]

            tag_group.create_dataset('origins', data=[mct.Origin()])
            tag_group.create_dataset('nparts', data=[mct.NParticles()])
            pdgs = []
            positions = []
            momenta = []
            processes = []
            trackids = []
            for i in range(mct.NParticles()):
                part = mct.GetParticle(i)
                pdgs.append(part.PdgCode())
                positions.append([part.Position(0).X(), part.Position(0).Y(), part.Position(0).Z(), part.T()])
                momenta.append([part.Px(), part.Py(), part.Pz(), part.E()])
                processes.append(part.Process())
                trackids.append(part.TrackId())
            tag_group.create_dataset('PDGs', data=pdgs)
            tag_group.create_dataset('positions', data=positions)
            tag_group.create_dataset('momenta', data=momenta)
            tag_group.create_dataset('processes', data=processes)
            tag_group.create_dataset('trackids', data=trackids)

