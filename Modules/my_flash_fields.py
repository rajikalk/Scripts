#!/usr/bin/env python
#Created by Rajika Kuruwita, 2019
import yt
yt.enable_parallelism()

def _Neg_z(field, data):
    """
    returns the negative of the z-positions
    """
    return -1*data['z']

yt.add_field("Neg_z", function=_Neg_z, units=r"cm", sampling_type="local")
    
def _Neg_dz(field, data):
    """
    returns the negative of the dz
    """
    return -1*data['dz']
    
yt.add_field("Neg_dz", function=_Neg_z, units=r"cm", sampling_type="local")
