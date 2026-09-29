from utils import handle_header, provide_list
from classes import edepv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import h5py as h5


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--tags', type=str, nargs='+')
    parser.add_argument('-o', type=str, default='dump_edeps.h5')
    args = parser.parse_args()


    handle_header()
    provide_list(edepv)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_edeps = ev.getValidHandle[edepv]

    with h5.File(args.o, 'w') as f:
        for tag in args.tags:
            print('Processing', tag)
            tag_group = f.create_group(tag)

            edeps = get_edeps(InputTag(tag))
            photons = []
            electrons = []
            energy = []
            positions = []

            print(f'Iterating over {len(edeps.product())} edeps')
            for edep in edeps.product():
                photons.append([edep.NumFPhotons(), edep.NumSPhotons()])
                electrons.append(edep.NumElectrons())
                energy.append(edep.Energy())
                positions.append([
                    edep.MidPointX(),
                    edep.MidPointY(),
                    edep.MidPointZ(),
                    edep.T(),
                ])
                
            tag_group.create_dataset('photons', data=photons)
            tag_group.create_dataset('electrons', data=electrons)
            tag_group.create_dataset('energy', data=energy)
            tag_group.create_dataset('positions', data=positions)
