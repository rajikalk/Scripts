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

add_field("com", function=_com, units=r"gridunits")

def _RadfromCom(field, data):
    lu = 10000000000000.
    gl = lu/256.
    position = []
    for val_x in range(254):
        for val_y in range(254):
            for val_z in range(254):        
                pos = [data["x"][val_x], data["y"][val_y], data["z"][val_z]]
                position.append(pos)
    rel_pos = []
    for i in range(255):
        for j in range(255):
            for k in range(255):
                x = (position[i][0] - data["com"][0])
                y = (position[j][1] - data["com"][1])
                z = (position[k][2] - data["com"][2])
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
    
def _EscapeVelocity(field,data):
    Msun = 1.9891*(10**33.)
    lu = 10000000000000.     #length units in cm
    gl = lu/256.             #cm in grid length
    r = data['Radius']       #in cm
    M = data['ParticleMassMsun'][0]*Msun + data['ParticleMassMsun'][1]*Msun
    G = 6.67259*(10**(-8.))
    v_e = ((2.*G*M)/r)**(0.5)
    v_g = ((data['x-velocity']**2) + (data['y-velocity']**2) + (data['z-velocity']**2))**(0.5)
    v = v_g/v_e
    if v_g.any() > v_e.any():
        v = -1.e-46
    else:
        v = v_g/v_e
    return v

add_field("EscapeVelocity", function=_EscapeVelocity, units=r"v/v_e")

def _RelativeKeplerianVelocity(field, data):
    v = (data['x-velocity']**2 + data['y-velocity']**2 + data['z-velocity']**2)**(0.5)
    v_rel = v/data['KeplerianVelocity']
    return v_rel
    
add_field("RelativeKeplerianVelocity", function=_RelativeKeplerianVelocity, units=r"#")
    
def _Bound(field, data):
    Msun = 1.9891*(10**33.)
    mass = data["CellMass"]
    #get thermal energy:
    TE = data["ThermalEnergy"]*mass
    #get kinetic energy:
    v = ((data['x-velocity']**2) + (data['y-velocity']**2) + (data['z-velocity']**2))**(0.5)
    KE = 0.5*mass*(v**2.)
    #get potential energy:
    GPEg = data["Grav_Potential"]*(mass)
    pp = data['ParticleMassMsun'][0]*Msun + data['ParticleMassMsun'][1]*Msun
    top = -6.67259*(10**(-8.))*pp*mass
    bottom = data["Radius"]
    GPEp = top/bottom
    total = (-1)*(KE + GPEg + GPEp + TE)
    '''
    if total.all() > 0:
        total = data["Density"]
    else:
        total = total
    '''
    return total
    
add_field("Bound", function=_Bound, units=r"g/cm^3")

def _DensityBound(field, data):
    if data["Bound"].any() < 0:
        result = data["Bound"]
    else:
        result = data["Density"]
    return result

add_field("DensityBound", function=_DensityBound, units=r"g/cm^3")

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6vkepfix_0.75k/DD00*/CE00*.hierarchy")

#save directory
#save_directory = 'xz-plane/'
it = 0

#create plots
for pf in ts:
    p = SlicePlot(pf, "y", "Density")
    p.annotate_velocity()
    p.annotate_particles(10.0, p_size = 50)
    p.annotate_contour("KeplerianVelocity", ncont=1, clim=(1.0, 1.0))
    time = str(pf.current_time*pf["years"])
    title = "Current Time:" + time + "years"
    p.annotate_title(title)
    p.set_zlim('all', zmin = 1.e-9, zmax = 1.e-4)
    filename = "frame_y_" + ("%03d" % it) + ".png"
    p.save(filename)
    it = it + 1
