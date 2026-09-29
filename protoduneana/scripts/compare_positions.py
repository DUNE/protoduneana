from argparse import ArgumentParser as ap
import os
import h5py as h5
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#(horizontal column, vertical column) into the (N, 4) x/y/z/t position arrays
VIEWS = [(0, 1), (2, 1), (2, 0)]
AXNAMES = ['X', 'Y', 'Z']


def split_spec(spec, kind):
    """Split a 'file|tag|label' spec. The label is optional. Note that art
       tags can contain colons ('IonAndScint:priorSCE'), hence the pipe."""
    parts = spec.split('|')
    if len(parts) == 2:
        fname, tag = parts
        label = '%s %s' % (os.path.basename(fname), tag)
    elif len(parts) == 3:
        fname, tag, label = parts
    else:
        raise ValueError('%s spec %r is not of the form file|tag[|label]' %
                         (kind, spec))

    with h5.File(fname, 'r') as f:
        if tag not in f:
            raise KeyError('Tag %r not found in %s -- has %s' %
                           (tag, fname, list(f.keys())))
    return fname, tag, label


def read_positions(fname, tag):
    """Positions as a flat (Npoints, 4) array of x, y, z, t.

       dump_edeps.py writes one row per deposit. dump_larg4.py writes one row
       per trajectory point of every particle -- either already flattened, or
       as a variable-length row per particle in older dumps, which is what the
       object dtype below picks up."""
    with h5.File(fname, 'r') as f:
        pos = f[tag]['positions'][:]

    if pos.dtype == object:
        pos = ([np.asarray(p, dtype='float64').reshape(-1, 4) for p in pos]
               if len(pos) else [np.empty((0, 4))])
        pos = np.concatenate(pos)
    return np.asarray(pos, dtype='float64').reshape(-1, 4)


def flatten(specs):
    #The inputs take both repeated flags and several specs per flag
    return [s for group in (specs if specs else []) for s in group]


def in_box(pos, lims):
    """Mask of the points inside the box the limits make. The limits given for
       any one coordinate cut every view, not just the ones that draw it, so
       e.g. a --zlim keeps the Y vs X view to that same slab."""
    mask = np.ones(len(pos), dtype=bool)
    for col, lim in enumerate(lims):
        if lim is not None:
            mask &= ((pos[:, col] >= min(lim)) & (pos[:, col] <= max(lim)))
    return mask


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('--g4-input', type=str, nargs='+', action='append',
                        metavar='FILE|TAG[|LABEL]',
                        help='dump_larg4.py output, repeatable (quote it)')
    parser.add_argument('--edep-input', type=str, nargs='+', action='append',
                        metavar='FILE|TAG[|LABEL]',
                        help='dump_edeps.py output, repeatable (quote it)')
    parser.add_argument('-o', type=str, default='compare_positions.png')
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--alpha', type=float, default=1.)
    parser.add_argument('--xlim', type=float, nargs=2, metavar=('MIN', 'MAX'))
    parser.add_argument('--ylim', type=float, nargs=2, metavar=('MIN', 'MAX'))
    parser.add_argument('--zlim', type=float, nargs=2, metavar=('MIN', 'MAX'))
    args = parser.parse_args()

    #Indexed like the position columns, so each view can look up the range of
    #whichever coordinate it puts on which axis
    lims = [args.xlim, args.ylim, args.zlim]

    specs = ([(s, 'g4') for s in flatten(args.g4_input)] +
             [(s, 'edep') for s in flatten(args.edep_input)])
    if not specs:
        parser.error('Provide at least one --g4-input or --edep-input')

    if not args.show:
        matplotlib.use('Agg')

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for spec, kind in specs:
        fname, tag, label = split_spec(spec, kind)
        pos = read_positions(fname, tag)
        npoints = len(pos)
        pos = pos[in_box(pos, lims)]
        print('%-5s %s [%s] -- %d of %d points in box' %
              (kind, fname, tag, len(pos), npoints))

        for ax, (h, v) in zip(axes, VIEWS):
            ax.scatter(pos[:, h], pos[:, v], s=.5, alpha=args.alpha,
                       label=label, rasterized=True)

    for ax, (h, v) in zip(axes, VIEWS):
        ax.set_xlabel('%s [cm]' % AXNAMES[h])
        ax.set_ylabel('%s [cm]' % AXNAMES[v])
        ax.set_aspect('equal')
        if lims[h] is not None:
            ax.set_xlim(*lims[h])
        if lims[v] is not None:
            ax.set_ylim(*lims[v])

    handles, labels = axes[0].get_legend_handles_labels()
    #Keep the key readable no matter how transparent the points themselves are
    for handle in handles:
        handle.set_alpha(1.)
    fig.legend(handles, labels, loc='upper center', markerscale=20,
               ncol=min(len(labels), 4))
    fig.tight_layout(rect=(0, 0, 1, .92))

    #bbox_inches keeps the axis labels when equal aspect leaves the panels
    #different heights than tight_layout planned for
    fig.savefig(args.o, dpi=150, bbox_inches='tight')
    print('Wrote', args.o)
    if args.show:
        plt.show()
