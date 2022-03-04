#!/usr/bin/env python
# -*- coding: utf-8 -*-

# so that python 2.x and 3.x will work the same
from __future__ import print_function, absolute_import, division

import astropy.units as unit

import numpy as np

class Nbody(object):

  def __init__(self,*args,**kw):
    
    # ------- Handle keywords
    self.Grho  = kw.get('Grho',100.)
    self.rsoft = kw.get('rsoft',0.3/2**14)
    
    mass1 = kw.get('mass1',1.)   # Msun
    mass2 = kw.get('mass2',1.)   # Msun
    rmin  = kw.get('rmin',None)  # AU
    rmax  = kw.get('rmax',None)  # AU
    e     = kw.get('e',None)
    a     = kw.get('a',None)     # semi major axis (AU)
    j     = kw.get('j',None)     # specific relative angular momentum (cm2/s)
    J     = kw.get('J',None)     # total relative angular momentum (g cm2 / s)
    
    self.scale_t = kw.get('scale_t',711482453942051.0)
    self.scale_l = kw.get('scale_l',1.234271032e+19)
    self.scale_m = kw.get('scale_m',6.18296319475e+36)
    
    mass1, mass2 = mass1*unit.Msun.to(self.code_mass), mass2*unit.Msun.to(self.code_mass)

    self.currentPosition = np.zeros((2,3),dtype=np.float64)
    self.currentVelocity = np.zeros((2,3),dtype=np.float64)
    self.currentMass     = np.array([mass1,mass2],dtype=np.float64)
    self.currentTime     = 0.

    if rmin is not None and e is not None:
      rmin = rmin*unit.AU.to(self.code_length)
      c    = rmin * (1 + e) # semi lactus rectum
      J    = np.sqrt(c*self.Grho*mass1*mass2*self.reducedMass(self.currentMass))
      rmax = c/(1 - e)
    else:
      raise ValueError("Not enough input parameters has been set to calculate orbit")
    
    # initialize
    rx01 =  rmax*mass2/(mass1+mass2)
    rx02 = -rmax*mass1/(mass1+mass2)
    vy01 =  J/(mass1*rmax)
    vy02 = -J/(mass2*rmax)
    
    self.currentPosition[0,0] = rx01
    self.currentPosition[1,0] = rx02
    self.currentVelocity[0,1] = vy01
    self.currentVelocity[1,1] = vy02
    
    self.massList = []
    self.positionList = []
    self.velocityList = []
    self.timeList = []
  
  # ------- units
  @property
  def code_time(self):
    return unit.def_unit('code_time', self.scale_t * unit.s)

  @property
  def code_mass(self):
    return unit.def_unit('code_mass', self.scale_m * unit.g)

  @property
  def code_length(self):
    return unit.def_unit('code_length', self.scale_l * unit.cm)

  @property
  def code_velocity(self):
    return unit.def_unit('code_velocity', self.scale_l/self.scale_t * unit.cm/unit.s)

  @property
  def code_accretion(self):
    return unit.def_unit('code_accretion', self.scale_m/self.scale_t * unit.g/unit.s)
  
  # ------- Force and potential functions
  def newtonianGravitationalPotential(self,r):
    """
    Newtonian gravitational potential
    """
    return -1./r
  
  def softenedGravitationalPotential(self,r,rsoft):
    """
    Softened gravitational potential
    """
    
    eps = r/rsoft
    
    cnst1 = -(14/15 + 13/15 + 1)/rsoft
    cnst2 = -(11/5 + 1)/rsoft
    
    f1 = 4./rsoft*(8/5*eps**5 - 12/5*eps**4 + 4/3*eps**2) + cnst1
    f2 = 4./rsoft*(-8/15*eps**5 + 12/5*eps**4 - 4*eps**3 + 8/3*eps**2 + 1/60/eps) + cnst2
    f3 = -1./r
    
    return np.select([eps < 0.5, eps < 1., eps >= 1],[f1,f2,f3])

  def gravitationalPotential(self,position,mass1,mass2,rsoft=0.,cyclic=False,period=1.):
    """
    gravitational potential between particles
    Input:
      position : Cartesian particle positions. numpy array of shape (nparticle,3)
      position : Particle masses. numpy array of shape (nparticle,)
    keywords:
      rsoft    : softening length. If 0 newtonian acceleration is used
    """
    
    relPos = self.relativePosition(position,cyclic=cyclic,period=period)
    distance = np.linalg.norm(relPos,axis=-1)
    
    if rsoft == 0.:
      f = self.newtonianPotential(distance)
    else:
      f = self.softenedGravitationalPotential(distance,rsoft)
    
    return self.Grho * f * mass1*mass2

  def newtonianAcceleration(self,r):
    """
    Newtonian acceleration given distance, r
    """
    return r**(-2)
  
  def softenedAcceleration(self,r,rsoft):
    """
    Softened accelration given distance r and softening length rsoft
    See Federath et al. (2010)
    """
    eps = r/rsoft
  
    f1 = 4./rsoft**2*(8/3*eps - 48/5*eps**3 + 8*eps**4)
    f2 = 4./rsoft**2*(16/3*eps - 12*eps**2 + 48/5*eps**3 - 8/3*eps**4 - 1/60/eps**2)
    f3 = 1/r**2
  
    return np.select([eps < 0.5, eps < 1., eps >= 1],[f1,f2,f3])
    
  def gravitationalAcceleration(self,position,mass,rsoft=0.,cyclic=False,period=1.):
    """
    gravitational acceleration on particles
    Input:
      position : Cartesian particle positions. numpy array of shape (nparticle,3)
      position : Particle masses. numpy array of shape (nparticle,)
    keywords:
      rsoft    : softening length. If 0 newtonian acceleration is used
    """
    
    relPos = self.relativePosition(position,cyclic=cyclic,period=period)
    distance = np.linalg.norm(relPos,axis=-1)
    
    if rsoft == 0.:
      f = self.newtonianAcceleration(distance)
    else:
      f = self.softenedAcceleration(distance,rsoft)
    
    return self.Grho * f[:,:,np.newaxis] * mass[np.newaxis,:,np.newaxis] * relPos / distance[:,:,np.newaxis]
  
  # ------- Helper functions
  def relativePosition(self,x,cyclic=False,period=1.):
    """
    Relative position between two points
    
    Input:
      x1, x2 : The two points for which to calculate the relative position.
    Keywords:
      cyclic : Use cyclic position cordinates (Default: False)
      period : period of cyclic coordinates. Only relevant if cyclic is True (Default: 1)
    """
  
    # ------- Sanity checks
    x = np.array(x)
    
    x1, x2 = x[:,np.newaxis], x[np.newaxis,:]
  
    # ------- Relative Euclidian distances
    if cyclic:
      return (((x2 - x1 + period/2.) % period) - period/2.)
    else:
      return x2 - x1
  
  # ------- Quantities
  @property
  def time(self):
    """
    Time array
    """
    if len(self.timeList) == 0:
      return None
    else:
      return np.hstack(self.timeList) * self.code_time
  
  @property
  def position(self):
    """
    position array
    """
    if len(self.positionList) == 0:
      return None
    else:
      return np.transpose(np.dstack(self.positionList),axes=(2,0,1)) * self.code_length
  
  @property
  def velocity(self):
    """
    velocity array
    """
    if len(self.velocityList) == 0:
      return None
    else:
      return np.transpose(np.dstack(self.velocityList),axes=(2,0,1)) * self.code_velocity
  
  @property
  def mass(self):
    """
    mass array
    """
    if len(self.massList) == 0:
      return None
    else:
      return np.vstack(self.massList) * self.code_mass
  
  # ------- Derived quantities
  def reducedMass(self,*args):
    """
    reduced mass
    """
    
    if len(args) == 0:
      mass = self.mass
    else:
      mass, = args
    
    if mass.ndim == 2:
      m1, m2 = mass[:,0], mass[:,1]
    elif mass.ndim == 1:
      m1, m2 = mass
    else:
      raise ValueError("mass array must be one- or two-dimensional")
    
    return m1*m2/(m1+m2)
  
  def q(self,*args):
    """
    Mass ratio
    """
    
    if len(args) == 0:
      mass = self.mass
    else:
      mass, = args
    
    return mass.min(axis=-1)/mass.max(axis=-1)
  
  def separation(self,*args):
    """
    calculate separation
    """
    
    if len(args) == 0:
      pos = self.position
    else:
      pos, = args
    
    if pos.ndim == 3:
      p1, p2 = pos[:,0,:], pos[:,1,:]
    elif pos.ndim == 2:
      p1, p2 = pos[0,:], pos[1,:]
    else:
      raise ValueError("Position array must be two- or three-dimensional")
    
    from cyclic import dist

    if isinstance(p1,unit.quantity.Quantity):
      return dist(p1,p2,cyclic=False) * self.code_length
    else:
      return dist(p1,p2,cyclic=False)
  
  def centreOfMassVelocity(self,*args):
    """
    centre of mass velocity
    """
    if len(args) == 0:
      vel, mass = self.velocity, self.mass
    else:
      vel, mass = args
    
    if isinstance(mass,unit.quantity.Quantity):
      mass = mass.value
    
    if vel.ndim == 3:
      return (vel*mass[:,:,np.newaxis]).sum(axis=1)/mass.sum(axis=1)[:,np.newaxis]
    else:
      return (vel*mass[:,np.newaxis]).sum(axis=0)/mass.sum()
  
  def centreOfMassPosition(self,*args):
    """
    centre of mass
    """
    if len(args) == 0:
      pos, mass = self.position, self.mass
    else:
      pos, mass = args
    
    if isinstance(mass,unit.quantity.Quantity):
      mass = mass.value
    
    return (pos*mass[:,:,np.newaxis]).sum(axis=1)/mass.sum(axis=1)[:,np.newaxis]
  
  def Ekin(self,*args,absolute=False):
    """
    kinetic energy
    """
    
    if len(args) == 0:
      vel, mass = self.velocity, self.mass
    else:
      vel, mass = args
    
    if absolute:
      speed2 = (vel**2).sum(axis=-1)
    else:
      if vel.ndim == 3:
        speed2 = ((vel - self.centreOfMassVelocity(vel,mass)[:,np.newaxis,:])**2).sum(axis=-1)
      else:
        speed2 = ((vel - self.centreOfMassVelocity(vel,mass))**2).sum(axis=-1)
    
    return 0.5*(mass*speed2).sum(axis=-1)
  
  def Epot(self,*args):
    """
    Gravitational potential energy
    """
    
    if len(args) == 0:
      pos, mass = self.position, self.mass
    else:
      pos, mass = args
    
    sep = self.separation(pos)
    
    if isinstance(sep,unit.quantity.Quantity):
      f = self.softenedGravitationalPotential(sep.value,rsoft=self.rsoft) if self.rsoft > 0 else self.newtonianGravitationalPotential(sep.value)
      f = f/self.code_length
    else:
      f = self.softenedGravitationalPotential(sep,rsoft=self.rsoft) if self.rsoft > 0 else self.newtonianGravitationalPotential(sep)
    
    if isinstance(mass,unit.quantity.Quantity):
      return self.Grho * self.code_length**3/self.code_mass/self.code_time**2 * np.product(mass,axis=-1) * self.code_mass * f
    else:
      return self.Grho * np.product(mass,axis=-1) * f

  def Etot(self,*args):
    """
    Total kinetic plus gravitational potential energy
    """
    
    if len(args) == 0:
      pos, vel, mass = self.position, self.velocity, self.mass
    else:
      pos, vel, mass = args
    
    return self.Ekin(vel,mass) + self.Epot(pos,mass)

  def specificEkin(self,*args,absolute=False):
    """
    specific kinetic energy
    """
    
    if len(args) == 0:
      vel, mass = self.velocity, self.mass
    else:
      vel, mass = args
    
    return self.Ekin(vel,mass,absolute=absolute)/self.reducedMass(mass)
  
  def specificEpot(self,*args):
    """
    specific gravitational potential energy
    """
    
    if len(args) == 0:
      pos, mass = self.position, self.mass
    else:
      pos, mass = args
    
    return self.Epot(pos,mass)/self.reducedMass(mass)
  
  def specificEtot(self,*args):
    """
    Total specific kinetic plus specific gravitational potential energy
    """
    
    if len(args) == 0:
      pos, vel, mass = self.position, self.velocity, self.mass
    else:
      pos, vel, mass = args
    
    return (self.Ekin(vel,mass) + self.Epot(pos,mass))/self.reducedMass(mass)
  
  def relativeAngularMomentum(self,*args):
    """
    Relative angular momentum
    """
    if len(args) == 0:
      pos, vel, mass = self.position, self.velocity, self.mass
    else:
      pos, vel, mass = args
    
    if pos.ndim == 3:
      relPos = pos[:,1,:] - pos[:,0,:]
      relVel = mass[:,1][:,np.newaxis]*vel[:,1,:] - mass[:,0][:,np.newaxis]*vel[:,0,:]
    else:
      relPos = pos[1,:] = pos[0,:]
      relVel = mass[1]*vel[1,:] - mass[0]*vel[0,:]
    
    if isinstance(relPos,unit.quantity.Quantity):
      return np.cross(relPos.value,relVel.value,axis=-1) * self.code_length * self.code_mass * self.code_velocity
    else:
      return np.cross(relPos,relVel,axis=-1)
  
  def relativeAngularMomentumMagnitude(self,*args):
    """
    Relative angular momentum magnitude
    """
    if len(args) == 0:
      pos, vel, mass = self.position, self.velocity, self.mass
    else:
      pos, vel, mass = args
    
    return np.sqrt((self.relativeAngularMomentum(pos,vel,mass)**2).sum(axis=-1))

  def specificRelativeAngularMomentum(self,*args):
    """
    Specific Relative angular momentum
    """
    if len(args) == 0:
      pos, vel = self.position, self.velocity
    else:
      pos, vel = args
    
    if pos.ndim == 3:
      relPos = pos[:,1,:] - pos[:,0,:]
      relVel = vel[:,1,:] - vel[:,0,:]
    else:
      relPos = pos[1,:] = pos[0,:]
      relVel = vel[1,:] - vel[0,:]
    
    if isinstance(relPos,unit.quantity.Quantity):
      return np.cross(relPos.value,relVel.value,axis=-1) * self.code_length * self.code_velocity
    else:
      return np.cross(relPos,relVel,axis=-1)
  
  def specificRelativeAngularMomentumMagnitude(self,*args):
    """
    Specific Relative angular momentum magnitude
    """
    if len(args) == 0:
      pos, vel = self.position, self.velocity
    else:
      pos, vel = args
    
    return np.sqrt((self.specificRelativeAngularMomentum(pos,vel)**2).sum(axis=-1))
  
  def semiMajorAxis(self,*args):
    """
    semi major axis
    """
    
    if len(args) == 0:
      pos, vel, mass = self.position, self.velocity, self.mass
    else:
      pos, vel, mass = args
    
    etot = self.specificEtot(pos,vel,mass)
    
    if isinstance(mass,unit.quantity.Quantity):
      mu = self.Grho * self.code_length**3/self.code_mass/self.code_time**2 * mass.sum(axis=-1)
    else:
      mu = self.Grho * mass.sum(axis=-1)
    
    return -mu/(2*etot)
  
  def eccentricity(self,*args):
    """
    orbital eccentricity
    """
    
    if len(args) == 0:
      pos, vel, mass = self.position, self.velocity, self.mass
    else:
      pos, vel, mass = args
    
    etot = self.specificEtot(pos,vel,mass)
    h    = self.specificRelativeAngularMomentumMagnitude(pos,vel)
    mu   = self.Grho * self.code_length**3/self.code_mass/self.code_time**2 * mass.sum(axis=-1) if isinstance(mass,unit.quantity.Quantity) else self.Grho * mass.sum(axis=-1)
    
    return np.sqrt(1 + 2*etot*h**2/mu**2).value
  
  def orbitalPeriod(self,*args):
    """
    orbital Period
    """
    if len(args) == 0:
      pos, vel, mass = self.position, self.velocity, self.mass
    else:
      pos, vel, mass = args
    
    a = self.semiMajorAxis(pos,vel,mass)
    
    if isinstance(mass,unit.quantity.Quantity):
      mu = self.Grho * self.code_length**3/self.code_mass/self.code_time**2 * mass.sum(axis=-1)
    else:
      mu = self.Grho * mass.sum(axis=-1)
    
    return 2*np.pi * np.sqrt(a**3/mu)

  
  # ------- Integrate
  def updatePosition(self,position,velocity,dt,c=1.):
    return position + c * velocity * dt
  
  def updateVelocity(self,velocity,acceleration,dt,c=1.):
    return velocity + c * acceleration * dt
  
  def integrate(self,**kw):
    """
    nbody integration method
    Input:
      None
    Keywords:
      dt    : Time step. If not given the courant time step will be calculated and use.
              If scalar a constant time step will be used. If numpy array each iteration
              will step forward in the array and use that time step
      niter : Number of iterations. Will stop after this has been reached (Default: 1000)
      time  : integration will stop after this time has been reached. Takes presecense
              over niter (Default: None)
      rsoft : Softening length for gravitational potential (Default: 0)
    rerturns:
      position : particle positions at each iteration
      velocity : particle velocities at each iteration
      t        : time at each iteration
    """
    
    # ------- Keywords
    dt     = kw.get('dt',5e-8)
    niter  = kw.get('niter',1000)
    tmax   = kw.get('tmax',None) # yr
    nperiod = kw.get('nperiod',None)
    acc1    = kw.get('acc1',0.) # Msun/yr
    acc2    = kw.get('acc2',0.) # Msun/yr
    
    
    cyclic = False
    period = 1
    
    t0    = self.currentTime
    n     = 0
    
    if tmax is not None:
      tmax = tmax*unit.yr.to(self.code_time)
    
    if nperiod is not None:
      orbitalPeriod = self.orbitalPeriod(self.currentPosition,self.currentVelocity,self.currentMass)
      tmax          = nperiod*orbitalPeriod
    
    # ------- Write initial values to lists
    if len(self.timeList) == 0:
      self.timeList.append(self.currentTime)
    if len(self.positionList) == 0:
      self.positionList.append(self.currentPosition)
    if len(self.velocityList) == 0:
      self.velocityList.append(self.currentVelocity)
    if len(self.massList) == 0:
      self.massList.append(self.currentMass.copy())
    
    while True:
      # ------- increase mass
      if acc1 > 0:
        self.currentMass[0] += (acc1 * unit.Msun/unit.yr).to(self.code_mass/self.code_time).value*dt
      if acc2 > 0:
        self.currentMass[1] += (acc2 * unit.Msun/unit.yr).to(self.code_mass/self.code_time).value*dt

      ga       = self.gravitationalAcceleration(self.currentPosition,self.currentMass,rsoft=self.rsoft,cyclic=cyclic,period=period)
      
      # ------- set time step
#       courantdt = self.courantCondition(distance,ga,c=courantCoefficient)
#       ddt = min(courantdt,np.interp(self.currentTime,dttime,dt))
      
      # ------- Kick Drift Kick
      self.currentVelocity = self.updateVelocity(self.currentVelocity,np.nansum(ga,axis=1),dt,c=0.5) # kick
      self.currentPosition = self.updatePosition(self.currentPosition,self.currentVelocity,dt,c=1.)  # drift
      ga                   = self.gravitationalAcceleration(self.currentPosition,self.currentMass,rsoft=self.rsoft,cyclic=cyclic,period=period)
      self.currentVelocity = self.updateVelocity(self.currentVelocity,np.nansum(ga,axis=1),dt,c=0.5) # kick
      
      n += 1
      self.currentTime += dt
      
      # ------- yield results and add to lists
      #yield self.currentPosition, self.currentVelocity, self.currentTime, courantdt
      yield self.currentPosition, self.currentVelocity, self.currentMass, self.currentTime, dt
      self.timeList.append(self.currentTime)
      self.positionList.append(self.currentPosition)
      self.velocityList.append(self.currentVelocity)
      self.massList.append(self.currentMass.copy())
      
      # ------- Termination conditions
      if tmax is None:
        if n >= niter:
          break
      else:
        if self.currentTime - t0 > tmax:
          break
  
  def run(self,**kw):
    # ------- keywords
    countwrite = kw.get('countwrite',100)
    
    n = 0
    for p,v,m,t,dt in self.integrate(**kw):
      n += 1
      if n % countwrite == 0:
        print(n,t*self.code_time.to('yr'),m*self.code_mass.to('Msun'))
      continue

# class Nbody(object):
#   
#   def __init__(self,*args,**kw):
#     
#     # ------- Keywords
#     self.G      = kw.get('G',100.)
#     self.e      = kw.get('e',0.)
#     self.rmin   = kw.get('rmin',10.)
#     
#     # ------- Arguments
#     if len(args) == 1: # set initial values according to e and r0
#       self.currentMass = args[0]
# #       self.currentMass = args[0]*0.00033811191388411984
# #       self.G    = 39.48790616483102
#       
#       reducedMass      = self.currentMass[0]*self.currentMass[1]/(self.currentMass[0]+self.currentMass[1])
#       semiLactusRectum = self.rmin*(1+self.e)
#       semiMajorAxis    = semiLactusRectum/(1-self.e**2)
#       Vt               = np.sqrt(self.G*(self.currentMass[0] + self.currentMass[1])/semiLactusRectum)*(1-self.e)
#       
# #       l = a*(1-e**2)
# #       
# #       rmin = p/(1+e)
# #       p = 
# #       rmax = p/(1-e)
# #       
# #       P = 2*np.pi*np.sqrt(a**3/(G(m1+m2)))
# #       v = a/np.sqrt(a**3/(G(m1+m2)))
# #       a = G(m1+m2)/v**2
# #       
# #       a1 = G(m1+m2)/v**2 * m2/(m1+m2)
# #          = G*m2/v**2 
#       
#       # initialise arays
#       self.currentPosition = np.zeros((2,3),dtype=np.float64)
#       self.currentVelocity = np.zeros((2,3),dtype=np.float64)
#       
#       # set initial y positions
#       self.currentPosition[0,1] =  semiLactusRectum*self.currentMass[1]/(self.currentMass[0]+self.currentMass[1])
#       self.currentPosition[1,1] = -semiLactusRectum*self.currentMass[0]/(self.currentMass[0]+self.currentMass[1])
#       
#       self.currentVelocity[0,0] = -Vt*self.currentMass[1]/(self.currentMass[0]+self.currentMass[1])
#       self.currentVelocity[1,0] =  Vt*self.currentMass[0]/(self.currentMass[0]+self.currentMass[1])
#       
#     elif len(args) == 3:
#       self.currentPosition, self.currentVelocity, self.currentMass = args
#     else:
#       raise ValueError("So far only three inputs are supported")
#     
#     # ------- Constants
#     self.currentTime = 0.
#     
#     self.time     = []
#     self.position = []
#     self.velocity = []
#     self.mass     = []
#   
#   ######### Helper functions
#   def relativePosition(self,x,cyclic=False,period=1.):
#     """
#     Relative position between two points
#     
#     Input:
#       x1, x2 : The two points for which to calculate the relative position.
#     Keywords:
#       cyclic : Use cyclic position cordinates (Default: False)
#       period : period of cyclic coordinates. Only relevant if cyclic is True (Default: 1)
#     """
#   
#     # ------- Sanity checks
#     x = np.array(x)
#     
#     x1, x2 = x[:,np.newaxis], x[np.newaxis,:]
#   
#     # ------- Relative Euclidian distances
#     if cyclic:
#       return (((x2 - x1 + period/2.) % period) - period/2.)
#     else:
#       return x2 - x1
#   
#   def dist(self,x,**kw):
#     
#     return np.linalg.norm(self.relativePosition(x,**kw),axis=-1)
#   
#   ########## Force terms
#   def newtonianAcceleration(self,r):
#     """
#     Newtonian acceleration given distance, r
#     """
#     return r**(-2)
#   
#   def softenedAcceleration(self,r,rsoft):
#     """
#     Softened accelration given distance r and softening length rsoft
#     See Federath et al. (2010)
#     """
#     eps = r/rsoft
#   
#     f1 = 4./rsoft**2*(8/3*eps - 48/5*eps**3 + 8*eps**4)
#     f2 = 4./rsoft**2*(16/3*eps - 12*eps**2 + 48/5*eps**3 - 8/3*eps**4 - 1/60/eps**2)
#     f3 = 1/r**2
#   
#     return np.select([eps < 0.5, eps < 1., eps >= 1],[f1,f2,f3])
#     
#   def gravitationalAcceleration(self,position,mass,rsoft=0.,cyclic=False,period=1.):
#     """
#     gravitational acceleration on particles
#     Input:
#       position : Cartesian particle positions. numpy array of shape (nparticle,3)
#       position : Particle masses. numpy array of shape (nparticle,)
#     keywords:
#       rsoft    : softening length. If 0 newtonian acceleration is used
#     """
#     
#     relPos = self.relativePosition(position,cyclic=cyclic,period=period)
#     distance = np.linalg.norm(relPos,axis=-1)
#     
#     if rsoft == 0.:
#       f = self.newtonianAcceleration(distance)
#     else:
#       f = self.softenedAcceleration(distance,rsoft)
#     
#     return self.G * f[:,:,np.newaxis] * mass[np.newaxis,:,np.newaxis] * relPos / distance[:,:,np.newaxis]
#   
#   def courantCondition(self,distance,acceleration,c=0.5,dx=2**(-14)):
#     """
#     Courant condition for time step. See Federath et al. (2010)
#     """
#     
#     index = distance > 0
#     
#     #dt1 = c*np.min(dx/(2*np.linalg.norm(velocity/,axis=-1)))
#     
#     distance = np.where(distance > dx,dx,distance)
#     
#     dt2 = c*np.min(np.sqrt((distance/np.linalg.norm(acceleration,axis=-1))[index]))
#     #return min(dt1,dt2)
#     #return min(dt2,6e-8)
#     return dt2
#   
#   def updatePosition(self,position,velocity,dt,c=1.):
#     return position + c * velocity * dt
#   
#   def updateVelocity(self,velocity,acceleration,dt,c=1.):
#     return velocity + c * acceleration * dt
#   
#   def centreOfMass(self,x,mass,indices=None):
#     """
#     translation to centre of mass system
#     """
#     assert mass.ndim == 1, "mass array can only have one dimension"
#     assert x.shape[0] == mass.size, "m and x must have same first dimension"
#     
#     if indices is not None: # only cm over indicies
#       xcm = (mass[indices]*x[indices,:]).sum(axis=0)/mass[indices].sum()
#     else: # cm over everything
#       xcm = (mass[:,np.newaxis]*x).sum(axis=0)/mass.sum()
#     
#     return x - xcm
#   
#   def integrate(self,**kw):
#     """
#     nbody integration method
#     Input:
#       None
#     Keywords:
#       dt    : Time step. If not given the courant time step will be calculated and use.
#               If scalar a constant time step will be used. If numpy array each iteration
#               will step forward in the array and use that time step
#       niter : Number of iterations. Will stop after this has been reached (Default: 1000)
#       time  : integration will stop after this time has been reached. Takes presecense
#               over niter (Default: None)
#       rsoft : Softening length for gravitational potential (Default: 0)
#     rerturns:
#       position : particle positions at each iteration
#       velocity : particle velocities at each iteration
#       t        : time at each iteration
#     """
#     
#     # ------- Keywords
#     dt    = kw.get('dt',None)
#     dttime = kw.get('dttime',None)
#     niter = kw.get('niter',1000)
#     tmax  = kw.get('tmax',None)
#     rsoft = kw.get('rsoft',0.)
#     cyclic = kw.get('cyclic',False)
#     period = kw.get('period',1.)
#     courantCoefficient = kw.get('courantCoefficient',0.5)
#     ramses = kw.get('ramses',False)
#     
#     t0    = self.currentTime
#     n     = 0
#     
#     # ------- Write initial values to lists
#     if len(self.time) == 0:
#       self.time.append(self.currentTime)
#     if len(self.position) == 0:
#       self.position.append(self.currentPosition)
#     if len(self.velocity) == 0:
#       self.velocity.append(self.currentVelocity)
#     if len(self.mass) == 0:
#       self.mass.append(self.currentMass)
#     
#     if ramses:
#       ga = self.gravitationalAcceleration(self.currentPosition,self.currentMass,rsoft=rsoft,cyclic=cyclic,period=period)
#       self.currentVelocity = self.updateVelocity(self.currentVelocity,np.nansum(ga,axis=1),dt[0],c=0.5) # kick
#       dt = dt[1:]
#     
#     while True:
#       distance = self.dist(self.currentPosition,cyclic=cyclic,period=period)
#       ga       = self.gravitationalAcceleration(self.currentPosition,self.currentMass,rsoft=rsoft,cyclic=cyclic,period=period)
#       
#       # ------- set time step
#       courantdt = self.courantCondition(distance,ga,c=courantCoefficient)
#       ddt = min(courantdt,np.interp(self.currentTime,dttime,dt))
#       #if dt is None:
#         #ddt = self.courantCondition(distance,ga,c=courantCoefficient)
#       #elif isinstance(dt,np.ndarray):
#         #ddt = dt[n]
#       #else:
#         #ddt = dt
#       
#       # ------- Kick Drift Kick
#       self.currentVelocity = self.updateVelocity(self.currentVelocity,np.nansum(ga,axis=1),ddt,c=0.5) # kick
#       self.currentPosition = self.updatePosition(self.currentPosition,self.currentVelocity,ddt,c=1.)  # drift
#       ga                   = self.gravitationalAcceleration(self.currentPosition,self.currentMass,rsoft=rsoft,cyclic=cyclic,period=period)
#       self.currentVelocity = self.updateVelocity(self.currentVelocity,np.nansum(ga,axis=1),ddt,c=0.5) # kick
#       
#       n += 1
#       self.currentTime += ddt
#       
#       if cyclic:
#         self.currentPosition = self.currentPosition % period
#       
#       # ------- yield results and add to lists
#       #yield self.currentPosition, self.currentVelocity, self.currentTime, courantdt
#       yield self.currentPosition, self.currentVelocity, self.currentTime, ddt
#       self.time.append(self.currentTime)
#       self.position.append(self.currentPosition)
#       self.velocity.append(self.currentVelocity)
#       self.mass.append(self.currentMass)
#       
#       # ------- Termination conditions
#       if tmax is None:
#         if n >= niter:
#           break
#       else:
#         if self.currentTime - t0 > tmax:
#           break
#   
#   def run(self,**kw):
#     
#     # ------- Keywords
#     centreOfMass = kw.get('centreOfMass',False)
#     
#     pos, vel, time = [], [], []
#     for p,v,t in self.integrate(**kw):
#       continue
#     
#     return self.position, self.velocity, self.time
