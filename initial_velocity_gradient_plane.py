#relative keplering velocity

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

def distance(point1, point2):
#Takes in the position as grid points and give the separation in grid points
    x_diff_sq = (point1[0] - point2[0])**2.
    y_diff_sq = (point1[1] - point2[1])**2.
    z_diff_sq = (point1[2] - point2[2])**2.
    result = (((x_diff_sq) + (y_diff_sq) + (z_diff_sq))**(0.5))
    return result
    
def KeplerianVelocity(mass, radius):
#All units are in CGS
    v = ((6.67259e-8*mass)/radius)**0.5
    return v
    
#Define GPE function
def GPE(Mass, distance):
#Takes in the masses in g and the separation in cm
#gives energy in erg
    top = -6.67259*(10**(-8.))*Mass
    bottom = distance
    result = top/bottom
    return result

#calculate gradient:
def Grad(GPE2, GPE1, dx):
    grad = (GPE2 - GPE1)/dx
    return grad    
    
#load data cube:
pf = load ("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")
g = pf.h.grids[0]

#Find size of domain
dim = pf.domain_dimensions[0]
dim_rsun = pf.domain_width*pf['rsun']

#Define values:
Msun = pf['MassUnits']       #grams/Msun
lu = pf['cm']     #length units in cm
tu = pf['TimeUnits']           #time units in sec
DensityUnits = pf['DensityUnits']
gl = lu/dim             #cm in grid length

#Set x and y position you wish to find vertical velocity distribution for.
#x_pos = dim/2.
y_pos = dim/2.
z_pos = dim/2.

#Find the mass of the particles
dd = pf.h.all_data()
pm1 = dd["ParticleMassMsun"][0]*Msun
pm2 = dd["ParticleMassMsun"][1]*Msun
mass = pm1 + pm2

#define arrays:
#z_pos = []  #z position relative to orbital plane
x_pos = []  #x position along orbital plane
kep_v= []   #keplerian velocity at that height above adn below orbital plane
vel = []    #Initial velocity of gas above and below the orbital plane
rel_v = []  #Relative velocity of gas compared to keplerian velocity.

for x in range(dim):
    g = pf.h.grids[0] 
    x_pos.append(((x-dim/2.)/dim)*dim_rsun) #append z position. Take away 128 to have it relative to orbital plane.
    
    #get particle position:
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    
    #calculate distances to each particle
    radiusp1 = distance([(x+0.5)*gl, (y_pos+0.5)*gl, (z_pos+0.5)*gl],pp1)
    radiusp1x = distance([(x+1.5)*gl, (y_pos+0.5)*gl, (z_pos+0.5)*gl],pp1)
    radiusp1y = distance([(x+0.5)*gl, (y_pos+1.5)*gl, (z_pos+0.5)*gl],pp1)
    radiusp1z = distance([(x+0.5)*gl, (y_pos+0.5)*gl, (z_pos+1.5)*gl],pp1)
    radiusp2 = distance([(x+0.5)*gl, (y_pos+0.5)*gl, (z_pos+0.5)*gl],pp2)
    radiusp2x = distance([(x+1.5)*gl, (y_pos+0.5)*gl, (z_pos+0.5)*gl],pp2)
    radiusp2y = distance([(x+0.5)*gl, (y_pos+1.5)*gl, (z_pos+0.5)*gl],pp2)
    radiusp2z = distance([(x+0.5)*gl, (y_pos+0.5)*gl, (z_pos+1.5)*gl],pp2)
    radius = distance([(x+0.5)*gl, (y_pos+0.5)*gl, (z_pos+0.5)*gl], [0.5*lu, 0.5*lu, 0.5*lu])

    #cmass = g["CellMassMsun"][x,y_pos,z_pos]*Msun #in grams
    GPE_val = GPE(pm1, radiusp1) + GPE(pm2, radiusp2)
    GPEx = GPE(pm1, radiusp1x) + GPE(pm2, radiusp2x)
    GPEy = GPE(pm1, radiusp1y) + GPE(pm2, radiusp2y)
    GPEz = GPE(pm1, radiusp1z) + GPE(pm2, radiusp2z)
    
    #Calculate grad:
    dx = Grad(GPEx, GPE_val, gl)
    dy = Grad(GPEy, GPE_val, gl)
    dz = Grad(GPEz, GPE_val, gl)

    #calculate v_kep:
    mag = (dx**2. + dy**2. + dz**2.)**0.5
    v_kep = (radius*mag)**0.5
    
    kep_v.append(v_kep)
    v_x = g["x-velocity"][x,y_pos,z_pos]
    v_y = g["y-velocity"][x,y_pos,z_pos]
    v_z = g["z-velocity"][x,y_pos,z_pos]
    v_mag = (v_x**2 + v_y**2 + v_z**2)**0.5
    vel.append(v_mag)
    v_rel = v_mag/v_kep
    rel_v.append(v_rel)
    print "rel_v", v_rel
    
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
#f.clf()
f.set_size_inches(7,7)
ax1.plot(x_pos, kep_v, 'b')
ax1.set_ylabel('$v_{keplerian}$ ($cms^{-1}$)')
ax2.plot(x_pos, vel, 'b')
ax2.set_ylabel('$v$ ($cms^{-1}$)')
ax3.plot(x_pos, rel_v, 'b')
ax3.axhline(y=1.0, color='r', ls='--')
ax3.set_xlabel('x position ($R_\odot$)')
ax3.set_ylabel('$v/v_{keplerian}$')
f.subplots_adjust(hspace=0.1)

filename = 'initial_velocity_z_'+str(z_pos)+'_y_'+str(y_pos)+'.png'
plt.savefig(filename)
