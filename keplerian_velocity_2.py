#relative keplering velocity

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

def center_of_mass(mass1, mass2, position1, position2):
    #Everything in cgs
    x_top = mass1*position1[0] + mass2*position2[0]
    y_top = mass1*position1[1] + mass2*position2[1]
    z_top = mass1*position1[2] + mass2*position2[2]
    bot = mass1 + mass2
    x = x_top/bot
    y = y_top/bot
    z = z_top/bot
    com = [x,y,z]
    return com
    
def _com(field, data):
    x_top = data['particle_mass'][0]*data['particle_position_x'][0] + data['particle_mass'][1]*data['particle_position_x'][1]
    y_top = data['particle_mass'][0]*data['particle_position_y'][0] + data['particle_mass'][1]*data['particle_position_y'][1]
    z_top = data['particle_mass'][0]*data['particle_position_z'][0] + data['particle_mass'][1]*data['particle_position_z'][1]
    bot = data['particle_mass'][0] + data['particle_mass'][1]
    x = x_top/bot
    y = y_top/bot
    z = z_top/bot
    com = [x,y,z]
    return com

add_field("com", function=_com, units=r"gridunits")

def _RadfromCom(field, data):
    position = []
    for val in range(256^3):
        pos = [data["x"][val], data["y"][val], data["z"][val]]
        position.append(pos)
    rel_pos = []
    for i in range(256^3):
        x = position[i][0] - data["com"][0]
        y = position[i][1] - data["com"][1]
        z = position[i][2] - data["com"][2]
        r = (x**2 + y**2 + z**2)**0.5
        rel_pos.append(r)
        print rel_pos
    return rel_pos
    
add_field("RadfromCom", function=_RadfromCom, units=r"cm")
    
def _KeplerianVelocity(field, data):
    #Calculates the Keplerian Velocity is at point in cgs
    Msun = 1.9891*(10**33.)
    lu = 10000000000000.     #length units in cm
    gl = lu/256.             #cm in grid length
    r = data['Radius']*gl
    M = data['particle_mass'][0]*Msun + data['particle_mass'][1]*Msun
    G = 6.67259*(10**(-8.))
    v_k = ((G*M)/r)**(0.5)
    v_g = ((data['x-velocity']**2) + (data['y-velocity']**2) + (data['z-velocity']**2))**(0.5)
    v = v_g/v_k
    return v
    
add_field("KeplerianVelocity", function=_KeplerianVelocity, units=r"v/v_k")
    
def _EscapeVelocity(field,data):
    Msun = 1.9891*(10**33.)
    lu = 10000000000000.     #length units in cm
    gl = lu/256.             #cm in grid length
    r = data['Radius']*gl
    M = data['particle_mass'][0]*Msun + data['particle_mass'][1]*Msun
    G = 6.67259*(10**(-8.))
    v_e = ((2.*G*M)/r)**(0.5)
    v_g = ((data['x-velocity']**2) + (data['y-velocity']**2) + (data['z-velocity']**2))**(0.5)
    v = v_g/v_e
    return v

add_field("EscapeVelocity", function=_EscapeVelocity, units=r"v/v_e")

def _RelativeKeplerianVelocity(field, data):
    v = (data['x-velocity']**2 + data['y-velocity']**2 + data['z-velocity']**2)**(0.5)
    v_rel = v/data['KeplerianVelocity']
    return v_rel
    
add_field("RelativeKeplerianVelocity", function=_RelativeKeplerianVelocity, units=r"#")

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/media/DATA/Simulations/smallbox/run1.e-6lessetot/longer/DD00*/CE00*.hierarchy")

#save directory
save_directory = '~/../../../media/DATA/YT_output/run1.e-6lessetot/longer/xz-plane-Escape/'
it = 0

#create plots
for pf in ts:
    p = SlicePlot(pf, "y", "EscapeVelocity")
    p.annotate_velocity()
    p.annotate_particles(10.0, p_size = 50)
    time = str(pf.current_time*pf["years"])
    title = "Current Time:" + time + "years"
    p.annotate_title(title)
    #p.set_zlim('all', zmin = 1e-1, zmax = 1e2)
    filename = save_directory + "frame_" + str(it) + "_time_" + time + ".png"
    p.save(filename)
    it = it + 1
