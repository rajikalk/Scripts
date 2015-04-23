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
x_pos = dim/2.
y_pos = dim/2.
#z_pos = dim/2.

#Find the mass of the particles
dd = pf.h.all_data()
pm1 = dd["ParticleMassMsun"][0]*Msun
pm2 = dd["ParticleMassMsun"][1]*Msun
mass = pm1 + pm2

#define arrays:
z_pos = []  #z position relative to orbital plane
#x_pos = []  #x position along orbital plane
kep_v= []   #keplerian velocity at that height above adn below orbital plane
vel = []    #Initial velocity of gas above and below the orbital plane
rel_v = []  #Relative velocity of gas compared to keplerian velocity.

for z in range(dim):
    g = pf.h.grids[0] 
    z_pos.append(((z-dim/2.)/dim)*dim_rsun) #append z position. Take away 128 to have it relative to orbital plane.
    
    #get particle position:
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    
    #calculate distances to each particle
    radiusp1 = distance([(x_pos+0.5)*gl, (y_pos+0.5)*gl, (z+0.5)*gl],pp1)
    radiusp1x = distance([(x_pos+1.5)*gl, (y_pos+0.5)*gl, (z+0.5)*gl],pp1)
    radiusp1y = distance([(x_pos+0.5)*gl, (y_pos+1.5)*gl, (z+0.5)*gl],pp1)
    radiusp1z = distance([(x_pos+0.5)*gl, (y_pos+0.5)*gl, (z+1.5)*gl],pp1)
    radiusp2 = distance([(x_pos+0.5)*gl, (y_pos+0.5)*gl, (z+0.5)*gl],pp2)
    radiusp2x = distance([(x_pos+1.5)*gl, (y_pos+0.5)*gl, (z+0.5)*gl],pp2)
    radiusp2y = distance([(x_pos+0.5)*gl, (y_pos+1.5)*gl, (z+0.5)*gl],pp2)
    radiusp2z = distance([(x_pos+0.5)*gl, (y_pos+0.5)*gl, (z+1.5)*gl],pp2)
    radius = distance([(x_pos+0.5)*gl, (y_pos+0.5)*gl, (z+0.5)*gl], [0.5*lu, 0.5*lu, 0.5*lu])

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
    v_x = g["x-velocity"][x_pos,y_pos,z]
    v_y = g["y-velocity"][x_pos,y_pos,z]
    v_z = g["z-velocity"][x_pos,y_pos,z]
    v_mag = (v_x**2 + v_y**2 + v_z**2)**0.5
    vel.append(v_mag)
    v_rel = v_mag/v_kep
    rel_v.append(v_rel)
    print "rel_v", v_rel
    
f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
#f.clf()
f.set_size_inches(7,7)
ax1.plot(kep_v, z_pos, 'b')
ax1.set_xlabel('$v_{keplerian}$ ($cms^{-1}$)')
ax1.set_ylabel('z position ($R_\odot$)')
ax2.plot(vel, z_pos, 'b')
ax2.set_xlabel('$v$ ($cms^{-1}$)')
ax3.plot(rel_v, z_pos, 'b')
ax3.axvline(x=1.0, color='r', ls='--')
ax3.set_xlabel('$v/v_{keplerian}$')
f.autofmt_xdate(rotation=90)
f.subplots_adjust(wspace=0.1)

filename = 'initial_velocity_x_'+str(x_pos)+'_y_'+str(y_pos)+'.png'
plt.savefig(filename)
