from yt.mods import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def _Temperature(field, data):
    #This find the temperature, assuming the gas is idea (ie gamma = 5/3)
    #Thermal energy is given in erg/g,
    #Density is given in g/cm^3
    #Number density is given in 1/cm^3
    gamma = data.pf['Gamma']
    Boltzmann_constant = 1.3806488*(10**(-16)) #erg/K
    top= data["ThermalEnergy"] * data["Density"] * (gamma - 1.)
    bottom = data["NumberDensity"] * Boltzmann_constant
    temperature = top/bottom
    inf = np.all(np.isfinite(data["NumberDensity"]))
    return temperature
    
add_field("Temperature", function=_Temperature, units=r"K")

fns = ['/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000','/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0001/CE0001', '/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0002/CE0002', '/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0003/CE0003']

fig = plt.figure()

# See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
# These choices of keyword arguments produce a four panel plot with a single
# shared narrow colorbar on the right hand side of the multipanel plot. Axes
# labels are drawn for all plots since we're slicing along different directions
# for each plot.
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 2),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")

for i, fn in enumerate(fns):
    # Load the data and create a single plot
    pf = load(fn) # load data
    p = ProjectionPlot(pf, 'z', 'Density')

    # Ensure the colorbar limits match for all plots
    p.set_zlim('Density', 1e1, 1e7)

    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots['Density']
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Finally, this actually redraws the plot.
    p._setup_plots()

plt.savefig('multiplot_2x2_time_series.png')
