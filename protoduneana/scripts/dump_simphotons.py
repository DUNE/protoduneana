from utils import handle_header, provide_list
from classes import simphotv
from ROOT import vector, string
from ROOT.gallery import Event
from ROOT.art import InputTag
from argparse import ArgumentParser as ap
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox


class WaveformDisplay:
    """Asks for an optical channel and draws that channel's photon waveform
       from each tag. Stays open, so a new channel can be asked for at any
       time -- but a channel that's already on screen is left alone."""

    def __init__(self, get_phots, tags, tick_ns=512):
        self.get_phots = get_phots
        self.tags = tags
        self.tick_ns = tick_ns
        self.channel = None

        self.fig, self.ax = plt.subplots(figsize=(11, 6))
        self.fig.subplots_adjust(bottom=.2)

        self.box = TextBox(self.fig.add_axes([.12, .04, .1, .06]), 'Channel ')
        self.button = Button(self.fig.add_axes([.26, .04, .1, .06]), 'Confirm')
        self.status = self.fig.text(.4, .06, '')

        #Confirm and hitting enter in the box do the same thing. Clicking the
        #button also drops focus from the box, which fires on_submit as well,
        #but the second call asks for the channel already drawn and does nothing
        self.box.on_submit(self.submit)
        self.button.on_clicked(lambda event: self.submit(self.box.text))

        self.blank('Enter a channel number and hit Confirm')

    def set_status(self, message, color='k'):
        self.status.set_text(message)
        self.status.set_color(color)
        self.fig.canvas.draw_idle()

    def blank(self, message):
        self.ax.clear()
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.text(.5, .5, message, transform=self.ax.transAxes,
                     ha='center', va='center', color='grey')
        self.fig.canvas.draw_idle()

    def submit(self, text):
        text = text.strip()
        if not text:
            return
        try:
            channel = int(text)
        except ValueError:
            self.set_status('%r is not a channel number' % text, 'firebrick')
            return

        if channel == self.channel:
            #Already on screen, so there's nothing to redo
            return
        self.draw(channel)

    def read_channel(self, tag, channel):
        """Ticks and photon counts for one channel of one tag. Only the map of
           the requested channel gets walked -- the rest of the tag's
           SimPhotonsLite are only asked for their channel number."""
        for phot in self.get_phots(InputTag(tag)).product():
            if phot.OpChannel != channel:
                continue

            detected = phot.DetectedPhotons
            ticks = np.empty(detected.size(), dtype='int64')
            counts = np.empty(detected.size(), dtype='int64')
            for i, (tick, nphotons) in enumerate(detected):
                ticks[i], counts[i] = tick, nphotons
            return ticks, counts
        return np.empty(0, dtype='int64'), np.empty(0, dtype='int64')

    def draw(self, channel):
        self.set_status('Reading channel %d ...' % channel)

        waveforms = {}
        for tag in self.tags:
            ticks, counts = self.read_channel(tag, channel)
            print('%-20s channel %d -- %d ticks, %d photons' %
                  (tag, channel, len(ticks), counts.sum()))
            if len(ticks):
                waveforms[tag] = (ticks // self.tick_ns, counts)

        self.channel = channel
        if not waveforms:
            self.blank('No photons on channel %d' % channel)
            self.set_status('Channel %d has no photons in any tag' % channel,
                            'firebrick')
            return

        #One grid across every tag so the waveforms line up and can be summed.
        #The photon times span the whole readout window, so this is only cheap
        #to build and draw because they've been binned into readout ticks
        lo = min(ticks.min() for ticks, _ in waveforms.values())
        hi = max(ticks.max() for ticks, _ in waveforms.values())
        edges = np.arange(lo, hi + 1) * self.tick_ns

        self.ax.clear()
        total = np.zeros(len(edges))
        for tag, (ticks, counts) in waveforms.items():
            waveform = np.bincount(ticks - lo, weights=counts,
                                   minlength=len(edges))
            total += waveform
            self.ax.step(edges, waveform, where='post', lw=1, label=tag)

        if len(waveforms) > 1:
            self.ax.step(edges, total, where='post', lw=1, ls='--', color='k',
                         label='Total')

        self.ax.set_xlabel('Time [ns]')
        self.ax.set_ylabel('Detected photons' if self.tick_ns == 1 else
                           'Detected photons / %d ns tick' % self.tick_ns)
        self.ax.set_title('Channel %d' % channel)
        #'best' has to look at every point drawn, which is slow on a waveform
        self.ax.legend(loc='upper right')
        self.set_status('Channel %d -- %d photons over %d tags' %
                        (channel, total.sum(), len(waveforms)))


if __name__ == '__main__':
    parser = ap()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-e', type=int, default=0)
    parser.add_argument('--tags', type=str, nargs='+', required=True)
    parser.add_argument('--tick-resolution', type=int, default=512,
                        help='ns per readout tick to sum the photons into')
    parser.add_argument('--channel', type=int, default=None,
                        help='draw this channel on startup')
    args = parser.parse_args()

    if args.tick_resolution < 1:
        parser.error('--tick-resolution must be at least 1 ns')

    handle_header()
    provide_list(simphotv)

    ev = Event(vector(string)(1, args.i))
    ev.goToEntry(args.e)
    get_phots = ev.getValidHandle[simphotv]

    if matplotlib.get_backend().lower() == 'agg':
        print('matplotlib is on the non-interactive Agg backend -- the display '
              'needs a GUI backend to take input (is DISPLAY set?)')

    display = WaveformDisplay(get_phots, args.tags, args.tick_resolution)
    if args.channel is not None:
        display.submit(str(args.channel))
    plt.show()
