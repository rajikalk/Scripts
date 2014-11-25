#Creates projection plots of each data dump.

from yt.mods import *
import matplotlib.pyplot as plt

def _com(field, data):
    lu = 10000000000000.
    gl = lu/256.
    x_top = data['particle_mass'][0]*data['particle_position_x'][0] + data['particle_mass'][1]*data['particle_position_x'][1]
    y_top = data['particle_mass'][0]*data['particle_position_y'][0] + data['particle_mass'][1]*data['particle_position_y'][1]
    z_top = data['particle_mass'][0]*data['particle_position_z'][0] + data['particle_mass'][1]*data['particle_position_z'][1]
    bot = data['particle_mass'][0] + data['particle_mass'][1]
    x = x_top/bot
    y = y_top/bot
    z = z_top/bot
    com = [x,y,z]
    return com

add_field("com", function=_com, units=r"cm")

def _Lz(field, data):
    lu = 10000000000000.
    gl = lu/256.
    Msun = 1.9891e+33
    x_top = data['particle_mass'][0]*data['particle_position_x'][0] + data['particle_mass'][1]*data['particle_position_x'][1]
    y_top = data['particle_mass'][0]*data['particle_position_y'][0] + data['particle_mass'][1]*data['particle_position_y'][1]
    bot = data['particle_mass'][0] + data['particle_mass'][1]
    x = x_top/bot
    y = y_top/bot
    xpos = (data['x']-x)*gl
    ypos = (data['y']-y)*gl
    mass = data['TotalMassMsun']
    xvel = data['x-velocity']
    yvel = data['y-velocity']
    Lz = xpos*mass*yvel - ypos*mass*xvel
    return Lz
    
add_field("Lz", function=_Lz, units=r"gcm^2s")

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/media/DATA/Simulations/smallbox/hot2/DD00*/CE00*.hierarchy")

save_directory = '/media/DATA/YT_output/hot2/xz-plane-Lz/'
it = 0

for pf in ts:
    p = SlicePlot(pf, "y", "Lz")
    p.annotate_velocity()
    p.annotate_particles(10.0, p_size = 50)
    #p.annotate_contour("Density", ncont = 10)
    time = str(pf.current_time*pf["years"])
    title = "Current Time:" + time + "years"
    p.annotate_title(title)
    #p.set_zlim('all', zmin = 1e5, zmax = 1e15)
    filename = save_directory + "frame_" + ("%03d" % it) + ".png"
    p.save(filename)
    it = it + 1
