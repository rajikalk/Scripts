
from yt.mods import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def _KeplerianVelocity(field, data):
    #Calculates the Keplerian Velocity is at point in cgs
    Msun = 1.9891*(10**33.)
    lu = 10000000000000.     #length units in cm
    gl = lu/256.             #cm in grid length
    r1 = ((data['x'] - data['particle_position_x'][0])**2. + (data['y'] - data['particle_position_y'][0])**2. + (data['z'] - data['particle_position_z'][0])**2.)**0.5
    M1 = data['ParticleMass'][0]
    G = 6.67259*(10**(-8.))
    GPE1 = (G*M1)/(r1*lu)
    r2 = ((data['x'] - data['particle_position_x'][1])**2. + (data['y'] - data['particle_position_y'][1])**2. + (data['z'] - data['particle_position_z'][1])**2.)**0.5
    M2 = data['ParticleMass'][1]
    GPE2 = (G*M2)/(r2*lu)
    GPE = GPE1 + GPE2
    v_k = (GPE)**(0.5)
    v_g = ((data['x-velocity']**2) + (data['y-velocity']**2) + (data['z-velocity']**2))**(0.5)
    v = v_g/v_k
    if v.any() > 1.0:
        v = 10.0
    return v
    
add_field("KeplerianVelocity", function=_KeplerianVelocity, units=r"v/v_k")

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

'''
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
    p.set_zlim('Density', 1e-4, 1e-2)

    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots['Density']
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]
    plt.xticks(rotation='vertical')
    # Finally, this actually redraws the plot.
    p._setup_plots()

plt.savefig('multiplot_2x2_time_series.png')

'''
import yt
from yt.mods import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (4, 2),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")

i = 0
for fn in fns:
    pf = load(fn) # load data
    #p1 = SlicePlot(pf, "z", "Density")
    #p2 = SlicePlot(pf, 'y', "Density")
    p1 = SlicePlot(pf, 'z', 'Temperature')
    p2 = SlicePlot(pf, 'y', 'Temperature')
    #p1.set_zlim('Density', 1.e-9, 1.e-4)
    #p2.set_zlim('Density', 1.e-9, 1.e-4)
    p1.set_zlim('Temperature', 1e1, 1e7)
    p2.set_zlim('Temperature', 1e1, 1e7)
    p1.annotate_velocity()
    p1.annotate_particles(10.0, p_size = 10)
    p2.annotate_velocity()
    p2.annotate_particles(10.0, p_size = 10)
    #p3.annotate_velocity()
    #p3.annotate_particles(10.0, p_size = 10)
    #plot1 = p1.plots['Density']
    plot1 = p1.plots['Temperature']
    plot1.figure = fig
    plot1.axes = grid[i].axes
    plot1.cax = grid.cbar_axes[i]
    p1._setup_plots()
    i = i + 1
    if i == 1:
	p2.annotate_contour("KeplerianVelocity", ncont=1, clim=(1.0, 1.0))
    #plot2 = p2.plots['Density']
    plot2 = p2.plots['Temperature']
    plot2.figure = fig
    plot2.axes = grid[i].axes
    plot2.cax = grid.cbar_axes[i]
    p2._setup_plots()
    i = i + 1
    #plot3 = p3.plots['Temperature']
    #plot3.figure = fig
    #plot3.axes = grid[i].axes
    #plot3.cax = grid.cbar_axes[i]
    #p3._setup_plots()
    #i = i + 1

plt.savefig('Temperature_fallback.png')

