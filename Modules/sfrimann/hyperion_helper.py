# -*- coding: utf-8 -*-
import multiprocessing
import numpy as np
import os

from hyperion.model.model import Model
from hyperion.grid import SphericalPolarGrid
from hyperion.util.interpolate import interp1d_fast
from hyperion.model.analytical_yso_model import AnalyticalYSOModel

def checkConvergence(x1,x2,percentile=99.,old_absolute=None):
  
  ratio = np.vstack((x1/x2,x2/x1)).max(axis=0)
  
  absolute = np.percentile(ratio,percentile)
  relative = max(absolute/old_absolute,old_absolute/absolute) if old_absolute is not None else None
  
  return absolute, relative

def find_last_iteration(file_handle):
    max_iteration = 0
    for group_name in file_handle:
        if "iteration" in group_name:
            iteration = int(group_name.split('_')[1])
            max_iteration = max(iteration, max_iteration)
    return max_iteration

def calculate_internal_heatsource(w,z,temperature,density,volumes,mstar,mu=2.279,alpha=0.01):
  """
  Function to calculate viscous heating
  Assumes viscous heating term from Eq. 4 of Min et al. (2011), http://adsabs.harvard.edu/abs/2011Icar..212..416M
  
  .. math::
    heating = 9/4 * \alpha * Pressure * Omega_kep
  
  
  The units in the file are erg/cm^3/s
  Input:
    rho        : density in g/cm^3
    T          : Temperature in K
    grid       : grid object
    parameters :  parameters object
  Keywords:
    modeldir : Directory in which to put the file (Default: current directory)
    binary   : if true writes a binary file instead of an ascii file (Default: False)
  Returns:
    None
  """
  
  from hyperion.util.constants import k, G, m_h
  
  p   = k*density*temperature/(mu*m_h) # ideal gas law. Pressure in erg/cm^3
  
  ww, zz = np.meshgrid(w,z)
  ww = ww[np.newaxis,:]
  
  Omega = np.sqrt(G*mstar/ww**3)
  
  heating = 9./4. * alpha * p * Omega # viscous heating in erg/cm^3/s
  
  heating = heating*volumes # viscous heating in erg/s
  lvisc   = heating.sum()

  return lvisc, heating

def hseq_profile(w, z, temperature, mstar, mu=2.279):
    """
    Compute the new (normalized) density profile
    corresponding to a given temperature profile

    Parameters
    ----------
    w : float
        The cylindrical radius at which to compute the profile (in cm)
    z : np.ndarray
        The z coordinates of the cells at radius w (in cm)
    temperature : np.ndarray
        The temperatures in the cells (in K)
    mstar : float
        The mass of the star (in g)
    """

    from hyperion.util.constants import G, m_h, k
    from hyperion.util.integrate import integrate, integrate_subset

    # Compute the integrand
    integrand = z / temperature / (w ** 2 + z ** 2) ** 1.5

    # Compute the integral for all cells
    # TODO - inefficient to compute integral from scratch - optimize
    i = np.array([integrate_subset(z, integrand, 0., zmax) for zmax in z])
    i[z < 0] = -i[z < 0]

    # Compute the factor for the integrand
    factor = G * mstar * mu * m_h / k

    # Compute the profile
    density = np.exp(-i * factor) / temperature

    # Normalize the density profile
    density = density / integrate(z, density)

    return density


# The mean molecular weight of H2 + He is given by:
#
# mu = 4 * (X + 1) / (X + 2)
#
# where X is the mass fraction of Helium to Hydrogren. Assuming
#
# X = 0.32476319350473615
#
# gives:
#
# mu = 2.279


def run_with_vertical_hseq(prefix, model, n_iter=10, mpi=False,
                           n_processes=multiprocessing.cpu_count(),
                           overwrite=False):
    """
    Run a model with vertical hydrostatic equilibrium.

    .. note:: this is an experimental function that is currently in
              development. Please use with care!

    The hydrostatic equilibrium condition is only applied to the disk
    components. The following requirements must be met:

    - The model should be an AnalyticalYSOModel
    - The model should be defined on a cylindrical polar grid
    - The stellar mass should be set
    - The model should include at least one disk

    The dust properties for the model can be specified as dust or dust+gas
    densities as this does not have an impact on this calculation - however,
    the hydrostatic equilibrium is computed assuming an H2 + He mix of gas
    (i.e. mu=2.279). Note that this calculation also ignores the effects of
    self-gravity in the disk, which might be important for more massive disks.

    Parameters
    ----------
    prefix : str
        The prefix for the output
    model : `~hyperion.model.analytical_yso_model.AnalyticalYSOModel`
        The model to run
    n_iter : int, optional
        The number of iterations to run the model for
    mpi : bool, optional
        Whether to run the model in parallel
    n_processes : int, optional
        The number of processes to use if ``mpi`` is ``True``
    overwrite : bool, optional
        Whether to overwrite previous files
    """

    from hyperion.grid import CylindricalPolarGrid
    from hyperion.model.model_output import ModelOutput
    from hyperion.util.integrate import integrate
    from hyperion.util.constants import lsun

    if not isinstance(model, AnalyticalYSOModel):
        raise TypeError("Can only run hydrostatic equilibrium for AnalyticalYSOModel instances")

    if model.grid['grid_type'] != 'cylindrical':
        raise TypeError("Can only run hydrostatic equilibrium for models with cylindrical polar grids")

    if model.star.mass is None:
        raise ValueError("Stellar mass needs to be defined for calculation of hydrostatic equilibrium")

    if len(model.disks) == 0:
        raise ValueError("Can only run hydrostatic equilibrium for models with disks")
    else:
        n_disks = len(model.disks)

    # Write out initial model
    model.write(prefix + '_00000.rtin', overwrite=overwrite, merge_if_possible=False)

    # Run the initial model
    mo = model.run(prefix + '_00000.rtout', overwrite=overwrite,
                   mpi=mpi, n_processes=n_processes)

    previous = prefix + '_00000.rtout'

    for iteration in range(1, n_iter + 1):

        # Read in output
        mo = ModelOutput(previous)

        # Extract the quantities
        g = mo.get_quantities()

        # Get wall positions
        rw, zw = g.w_wall, g.z_wall

        # Make a 2-d grid of wall positions
        R, Z = np.meshgrid(rw, zw)

        # Extract density and temperature
        density = g['density']
        temperature = g['temperature']

        # TODO: need to find a better way than just assuming the first n
        # density grids are disks

        for idisk in range(n_disks):

            # Vertically extrapolate temperatures
            for i in range(len(g.w)):
                for j in range(len(g.p)):
                    reset = temperature[idisk].array[j, :, i] < 1.
                    temperature[idisk].array[j, reset, i] = np.max(temperature[idisk].array[j, :, i])  # shouldn't be max, but will do for now

            # Compute new density
            for i in range(len(g.w)):
                for j in range(len(g.p)):
                    density[idisk].array[j, :, i] = hseq_profile(g.w[i], g.z, temperature[idisk].array[j, :, i], model.star.mass) * integrate(g.z, density[idisk].array[j, :, i])

        # Instantiate new model based on previous
        m = Model.read(previous)

        # Override the density
        m.grid['density'] = density

        # add internal heatsource
        lvisc, mp = calculate_internal_heatsource(m.grid.w,m.grid.z,temperature[0].array,density[0].array,m.grid.volumes,model.star.mass)
        if len(m.sources) == 1:
          m.add_map_source(luminosity=lvisc,map=mp,name='viscous_heating')
          with open('luminosity.txt','w') as f:
            print >> f, lvisc/lsun
        else:
          m.sources[-1].luminosity = lvisc
          m.sources[-1].map = mp

        # Write and run
        m.write('{0:s}_{1:05d}.rtin'.format(prefix, iteration), overwrite=overwrite)
        m.run('{0:s}_{1:05d}.rtout'.format(prefix, iteration),
              overwrite=overwrite, mpi=mpi, n_processes=n_processes)

        previous = '{0:s}_{1:05d}.rtout'.format(prefix, iteration)

def run_with_viscous_heating(prefix, model, n_iter=10, mpi=False,
                           n_processes=64,
                           overwrite=False,alpha=0.01):
    """
    Run a model with viscous heating.

    - The model should be an AnalyticalYSOModel
    - The stellar mass should be set
    - The model should include at least one disk

    The dust properties for the model can be specified as dust or dust+gas
    densities as this does not have an impact on this calculation - however,
    the hydrostatic equilibrium is computed assuming an H2 + He mix of gas
    (i.e. mu=2.279). Note that this calculation also ignores the effects of
    self-gravity in the disk, which might be important for more massive disks.

    Parameters
    ----------
    prefix : str
        The prefix for the output
    model : `~hyperion.model.analytical_yso_model.AnalyticalYSOModel`
        The model to run
    n_iter : int, optional
        The number of iterations to run the model for
    mpi : bool, optional
        Whether to run the model in parallel
    n_processes : int, optional
        The number of processes to use if ``mpi`` is ``True``
    overwrite : bool, optional
        Whether to overwrite previous files
    """

    from hyperion.grid import CylindricalPolarGrid
    from hyperion.model.model_output import ModelOutput
    from hyperion.util.integrate import integrate
    from hyperion.util.constants import lsun

    if not isinstance(model, AnalyticalYSOModel):
        raise TypeError("Can only run hydrostatic equilibrium for AnalyticalYSOModel instances")

    if model.grid['grid_type'] != 'cylindrical':
        raise TypeError("Can only run hydrostatic equilibrium for models with cylindrical polar grids")

    if model.star.mass is None:
        raise ValueError("Stellar mass needs to be defined for calculation of hydrostatic equilibrium")

    if len(model.disks) == 0:
        raise ValueError("Can only run hydrostatic equilibrium for models with disks")
    else:
        n_disks = len(model.disks)

    # Write out initial model
    model.write(prefix + '_00000.rtin', overwrite=overwrite, merge_if_possible=False)

    # Run the initial model
    print(n_processes)
    mo = model.run(prefix + '_00000.rtout', overwrite=overwrite,
                   mpi=mpi, n_processes=n_processes)

    previous = prefix + '_00000.rtout'

    for iteration in range(1, n_iter + 1):

        # Read in output
        mo = ModelOutput(previous)

        # Extract the quantities
        g = mo.get_quantities()

        # Extract density and temperature
        density = g['density']
        temperature = g['temperature']

        # Instantiate new model based on previous
        m = Model.read(previous)

        # add internal heatsource
        lvisc, mp = calculate_internal_heatsource(m.grid.w,m.grid.z,temperature[0].array,density[0].array,m.grid.volumes,model.star.mass,alpha=alpha)
        if len(m.sources) == 1:
          m.add_map_source(luminosity=lvisc,map=mp,name='viscous_heating')
          with open('luminosity.txt','w') as f:
            print(lvisc/lsun,file=f)
        else:
          m.sources[-1].luminosity = lvisc
          m.sources[-1].map = mp

        # Write and run
        m.write('{0:s}_{1:05d}.rtin'.format(prefix, iteration), overwrite=overwrite)
        m.run('{0:s}_{1:05d}.rtout'.format(prefix, iteration),
              overwrite=overwrite, mpi=mpi, n_processes=n_processes)

        previous = '{0:s}_{1:05d}.rtout'.format(prefix, iteration)
