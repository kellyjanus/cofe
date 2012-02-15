import os.path
import glob
from enthought.traits.api import HasTraits, Instance, Int, Enum, File, Directory
from enthought.traits.ui.api import View, Group, Item
from enthought.traits.ui.menu import ToolBar, Action
from enthought.enable.api import ColorTrait
from enthought.enable.component_editor import ComponentEditor
from enthought.chaco.api import marker_trait, Plot, ArrayPlotData, Legend
from enthought.chaco.tools.api import PanTool, ZoomTool, LegendTool, \
        TraitsTool, DragZoom
from enthought.pyface.api import ImageResource

from numpy import linspace, sin
import pyfits

def read_demod_file(filename):
    f = pyfits.open(filename)
    return f

class GroundStation(HasTraits):

    plot = Instance(Plot)
    plot_P = Instance(Plot)
    file = File(filter=[u'*.fits'])
    root_folder = Directory()

    ch = Enum(*['ch%d' % i for i in range(16)])
    latest = Action(name = "Second last file",
                action = "get_last_file",
                toolip = "Get last file",
                )

    prev = Action(name = "Previous",
                action = "get_prev",
                toolip = "Get prev file",
                )
    next = Action(name = "Next",
                action = "get_next",
                toolip = "Get next file",
                )
    traits_view = View(
        Group(
              Item('file'),
              Item('root_folder'),
              Item('ch'),
              Item('plot', editor=ComponentEditor(), show_label=False),
              Item('plot_P', editor=ComponentEditor(), show_label=False)
             ),
              width=800, height=1000, resizable=True, title="COFE Ground Station",
              toolbar=ToolBar(latest,prev,next),
              )
        
    @property
    def last_day_files(self):
       last_day = sorted(glob.glob(os.path.join(self.root_folder ,'*')))[-1]
       return sorted(glob.glob(os.path.join(last_day ,'*.fits')))

    def get_last_file(self):
       self.file = self.last_day_files[-2]

    def get_prev(self):
            self.file = self.last_day_files[max(0,self.last_day_files.index(self.file) - 1)]
    def get_next(self):
        try:
            self.file = self.last_day_files[self.last_day_files.index(self.file) + 1]
        except:
            self.file = self.last_day_files[-1] 

    def __init__(self):
        super(GroundStation, self).__init__()

    def _file_changed(self):
        self.data = read_demod_file(self.file)
        self.create_plot()

    def _ch_changed(self):
        self.create_plot()

    def get_rev(self):
        return self.data['rev'].data.field(0)

    def get_sci(self, tqu='T'):
        return self.data[self.ch].data.field(tqu)

    def create_plot(self):
        plot_data = ArrayPlotData(rev = self.get_rev())
        self.plot = Plot(plot_data,title='Temperature')
        
        plot_data.set_data('T', self.get_sci('T'))
        self.plot.plot(('rev','T'))

        color = {'Q':'red','U':'green'}
        self.plot_P = Plot(plot_data,title='Polarization %s' % color)

        for qu in ['Q','U']:
            plot_data.set_data(qu, self.get_sci(qu))
            self.plot_P.plot(('rev',qu),color=color[qu])

        for p in [self.plot,self.plot_P]:
            zoom = ZoomTool(p, tool_mode="box", always_on=False)
            p.overlays.append(zoom)
            p.tools.append(PanTool(p))



if __name__ == "__main__":
    g = GroundStation()
    g.configure_traits()
