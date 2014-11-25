# makes an animation

from yt.mods import *
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
from JSAnimation import IPython_display
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

prj = ProjectionPlot(load('/media/DATA/Simulations/run1.e-6/DD0000/CE0000'), 0, 'Density', weight_field='Density')
prj.set_zlim('Density',1e-32,1e-26)
fig = prj.plots['Density'].figure
fig.canvas = FigureCanvasAgg(fig)

# animation function.  This is called sequentially
def animate(i):
    pf = load('/media/DATA/Simulations/run1.e-6/DD%04i/CE%04i' % (i,i))
    prj._switch_pf(pf)

# call the animator.  blit=True means only re-draw the parts that have changed.
animation.FuncAnimation(fig, animate, frames=44, interval=200, blit=False)
