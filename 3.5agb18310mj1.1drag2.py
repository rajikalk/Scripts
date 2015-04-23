from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np
import math
import os

incr = 1
for i in range(1,49,incr):

  f=open('Cool_drag.csv','a')
  print str(i)
  if i > 99:
    DD0="/DD0"+str(i-incr)
    CE0="/CE0"+str(i-incr)
    DD="/DD0"+str(i)
    CE="/CE0"+str(i)
    DD2="/DD0"+str(i+incr)
    CE2="/CE0"+str(i+incr)
    if i-incr < 100:
      DD0="/DD00"+str(i-incr)
      CE0="/CE00"+str(i-incr)
    if i+incr == 1000:
      DD2="/DD1000"
      CE2="/CE1000"
  if i < 100:
    DD0="/DD00"+str(i-incr)
    CE0="/CE00"+str(i-incr)
    DD="/DD00"+str(i)
    CE="/CE00"+str(i)
    DD2="/DD00"+str(i+incr)
    CE2="/CE00"+str(i+incr)
    if i-incr < 10:
      DD0="/DD000"+str(i-incr)
      CE0="/CE000"+str(i-incr)
    if i+incr == 100:
      DD2="/DD0100"
      CE2="/CE0100"
  if i < 10:
    DD0="/DD000"+str(i-incr)
    CE0="/CE000"+str(i-incr)
    DD="/DD000"+str(i)
    CE="/CE000"+str(i)
    DD2="/DD000"+str(i+incr)
    CE2="/CE000"+str(i+incr)
    if i+incr == 10:
      DD2="/DD0010"
      CE2="/CE0010"

  pf=load('/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/'+DD+CE)
  pf0=load('/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/'+DD0+CE0)
  pf2=load('/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/'+DD2+CE2)
  G=6.67e-8
  time1=pf.current_time
  time2=pf2.current_time
  ce=pf.h.all_data()
  ce0=pf0.h.all_data()
  ce2=pf2.h.all_data()
  timeunit = ce.pf.parameters["TimeUnits"]
  lengthunit = ce.pf.parameters["LengthUnits"]
#  sl=ce.pf.parameters["PointMassSmoothingLength"]
  ppx=ce['particle_position_x']
  ppy=ce['particle_position_y']
  ppz=ce['particle_position_z']
  pvx=ce['particle_velocity_x']
  pvy=ce['particle_velocity_y']
  pvz=ce['particle_velocity_z']
  pm=ce['ParticleMass']
  phi=ce['Grav_Potential']*(lengthunit/timeunit)**2
  ppx2=ce2['particle_position_x']
  ppy2=ce2['particle_position_y']
  ppz2=ce2['particle_position_z']
  pvx0=ce0['particle_velocity_x']
  pvy0=ce0['particle_velocity_y']
  pvz0=ce0['particle_velocity_z']
  pvx2=ce2['particle_velocity_x']
  pvy2=ce2['particle_velocity_y']
  pvz2=ce2['particle_velocity_z']
  pm2=ce2['ParticleMass']
  phi2=ce2['Grav_Potential']*(lengthunit/timeunit)**2
  phiatp11=pf.h.find_field_value_at_point(['Grav_Potential'],[ppx[1],ppy[1],ppz[1]])[0]*(lengthunit/timeunit)**2
  phiatp12=pf2.h.find_field_value_at_point(['Grav_Potential'],[ppx2[1],ppy2[1],ppz2[1]])[0]*(lengthunit/timeunit)**2


  pv=(pvx[1]**2+pvy[1]**2+pvz[1]**2)**0.5
  pv0=(pvx0[1]**2+pvy0[1]**2+pvz0[1]**2)**0.5
  pv2=(pvx2[1]**2+pvy2[1]**2+pvz2[1]**2)**0.5

  psep1 = ((ppx[1]-ppx[0])**2+(ppy[1]-ppy[0])**2+(ppz[1]-ppz[0])**2)**0.5*lengthunit
  psep2 = ((ppx2[1]-ppx2[0])**2+(ppy2[1]-ppy2[0])**2+(ppz2[1]-ppz2[0])**2)**0.5*lengthunit

  
  Pe1=-G*pm[0]*pm[1]/psep1+phiatp11*pm[1]
  Ke1=0.5*pm[1]*pv**2
  Pe2=-G*pm2[0]*pm2[1]/psep2+phiatp12*pm2[1]
  Ke2=0.5*pm2[1]*pv2**2

  energydiff = (Pe1+Ke1)-(Pe2+Ke2)
  timediff = (time2-time1)*365*24*3600.0
  disttraveled = pv*timediff
  disttraveledavg = (pv+pv0+pv2)/3.0*timediff
  force = energydiff/disttraveled
  forceavg = energydiff/disttraveledavg
#  f.write(str(i)+"  "+str(force)+"  "+str(force/pm[1])+"  "+str(forceavg)+"  "+str(Pe1)+"\n")
  f.write(str(i)+"  "+str(force)+"  "+str(force/pm[1])+"  "+str(forceavg)+"\n")
#  print "force= " +str(force)+", acceleration= "+str(force/pm[1])
