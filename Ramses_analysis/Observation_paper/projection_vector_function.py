import numpy as np

def xy_rotation_matrix(theta):
    """
    creates rotation matrix for given angle theta along the xy plane
    [cos(theta), -sin(theta), 0]
    [sin(theta), cos(theta) , 0]
    [0         , 0          , 1]
    """
    rot = np.array([[np.cos(theta), -1*np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0,0,1]])
    return rot
    
def z_rotation_matrix(phi):
    """
    creates rotation matrix for given angle phi along the yz plane
    [1, 0       , 0        ]
    [0, cos(phi), -sin(phi)]
    [0, sin(phi),  cos(phi)]
    """
    rot = np.array([[1,0,0], [0, np.cos(phi), -1*np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
    return rot
    
def projected_vector(vector, proj_vector):
    """
    Calculates the position of vector projected onto proj_vector
    """
    proj_v_x = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[0]
    proj_v_y = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[1]
    proj_v_z = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[2]
    proj_v = np.array([proj_v_x,proj_v_y,proj_v_z]).T
    return proj_v

def Get_projection_vectors(star_pos_1, star_pos_2, projected_separation=200):
    """
    star_pos_1: Position of the primary star in au. A list, or an array, e.g. [x, y, z]
    star_pos_2: Position of the secondary star in au. A list, or an array, e.g. [x, y, z]
    projected_separation: Desired projected separation in au. Default: 200au
    """
    try:
        separation = star_pos_1 - star_pos_2
    except:
        separation = np.array(star_pos_1) - np.array(star_pos_2)
    separation_magnitude = np.sqrt(separation[0]**2 + separation[1]**2 + separation[2]**2)
    sepration_unit = separation/separation_magnitude
    
    #Angle around xy plane:
    sep_xy_mag = np.sqrt(separation[0]**2 + separation[1]**2)
    sep_xy_unit = np.array([separation[0]/sep_xy_mag, separation[1]/sep_xy_mag, 0])
    theta = np.arccos(np.dot(sep_xy_unit, np.array([0,1,0])))
    #Angle from z axis:
    phi = np.arccos(np.dot(sepration_unit, [0,0,1]))
    
    #Calculate angle alpha between projection and separation vector such that the projected sepation is what you selected (default 200AU)
    alpha = np.arcsin(projected_separation/separation_magnitude)
    proj_length = np.sqrt(separation_magnitude**2 - projected_separation**2)
    radius = proj_length*(projected_separation/separation_magnitude)
    sep_z = np.sqrt(proj_length**2 - radius**2)
    
    vectors_along_cone = np.array([[radius, 0, sep_z],\
                                   [radius/np.sqrt(2), radius/np.sqrt(2), sep_z],\
                                   [0, radius, sep_z],\
                                   [-1*radius/np.sqrt(2), radius/np.sqrt(2), sep_z],\
                                   [-radius, 0, sep_z],\
                                   [-1*radius/np.sqrt(2), -1*radius/np.sqrt(2), sep_z],\
                                   [0, -radius, sep_z],\
                                   [radius/np.sqrt(2), -1*radius/np.sqrt(2), sep_z]])
    
    
    #Figure out how to rotate the projection vectors to be the same orientation as the separation vector.
    #Rotate around XY plane
    z_rot = z_rotation_matrix(-1*phi)
    #z_rot_rev = z_rotation_matrix(phi)
    #xy_rot_rev = xy_rotation_matrix(-1*theta)
    theta_pos = np.dot(xy_rotation_matrix(theta), np.dot(z_rot, [0,0,1])) - sepration_unit
    theta_neg = np.dot(xy_rotation_matrix(-1*theta), np.dot(z_rot, [0,0,1])) - sepration_unit
    theta_pos_len = np.sqrt(np.sum(theta_pos**2))
    theta_neg_len = np.sqrt(np.sum(theta_neg**2))
    if theta_pos_len<theta_neg_len:
        xy_rot = xy_rotation_matrix(theta)
    elif theta_pos_len>theta_neg_len:
        xy_rot = xy_rotation_matrix(-1*theta)
    else:
        print("PROBLEM WITH FINDING CORRECT ROTATION")
        import pdb
        pdb.set_trace()
    projection_vectors = []
    north_vectors = []
    for vector in vectors_along_cone:
        if separation_magnitude > projected_separation:
            proj_vector = np.dot(xy_rot, np.dot(z_rot, vector))
        
            Proj_sep_proj = projected_vector(separation,proj_vector)
            sep_on_proj_length = np.sqrt(Proj_sep_proj[0]**2 + Proj_sep_proj[1]**2 + Proj_sep_proj[2]**2)
            calculated_proj_separation = np.sqrt(separation_magnitude**2 - sep_on_proj_length**2)
            if np.round(calculated_proj_separation) != projected_separation:
                print("CALCULATED PROJECTED SEPARATION IS", np.round(calculated_proj_separation), "FOR FILE", fn, "AT TIME", time_val)
                import pdb
                pdb.set_trace()
            projection_vectors.append(proj_vector)
            north_vector = separation - Proj_sep_proj
            north_vectors.append(north_vector)
        else:
            projection_vectors.append(np.array([np.nan, np.nan, np.nan]))
            north_vectors.append(np.array([np.nan, np.nan, np.nan]))
            
    unit_projction_vectors = (np.array(projection_vectors).T/np.sqrt(np.sum(np.square(projection_vectors), axis=1))).T
    unit_north_vectors = (np.array(north_vectors).T/np.sqrt(np.sum(np.square(north_vectors), axis=1))).T
    return unit_projction_vectors, unit_north_vectors
            
def convert_to_spherical(xyz):
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    r = np.sqrt(xy + xyz[:,2]**2)
    #theta = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    theta = np.arctan2(np.sqrt(xy), xyz[:,2])/np.pi
    #phi = np.arctan2(xyz[:,1], xyz[:,0]) # for elevation angle defined from XY-plane up
    phi = np.arctan2(xyz[:,1], xyz[:,0])/np.pi
    sph_vectors = np.array([r, theta, phi]).T
    return sph_vectors
    
