#!/usr/bin/env python
import binary_orbit as bo
import numpy as np
import csv
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import sys

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--input_file", help="what file do you want to read in?", type=str)
    parser.add_argument("-bayes_f", "--bayes_file", help="where do you want to read out the calculated bayes factors?", type=str, default='bayes_factors.csv')
    parser.add_argument("-restart", "--restart_calc", help="Do you want to restart calculations?", type=str, default='False')
    parser.add_argument("-v_group", "--group_velocity", help="Do you want to use the group velocity or the mean velocity?", type=str, default='True')
    parser.add_argument("-v_sigma", "--group_velocity_sigma", help="To how many sigma do you want to use", type=float, default=3.0)
    parser.add_argument("-n_orbs", "--no_orbits", default=1e5, type=float)
    parser.add_argument("-n_sys", "--no_systems", default=1e2, type=float)
    parser.add_argument("-astrostd", "--astrophysical_std", default=1.0, type=float)
    args = parser.parse_args()
    return args

def __unicode__(self):
    return unicode(self.some_field) or u''

#=======MAIN=======
def main():
    
    rank = CW.Get_rank()
    size = CW.Get_size()

    args = parse_inputs()

    n_orb=int(args.no_orbits)
    n_systems = int(args.no_systems)
    q_min = 0.05
    my_orb = bo.random_orbits(n_orb=n_orb)
    US_group_vel = 10.
    UCL_group_vel = 4.
    #Madsen, 2002 gives the STANDARD ERROR of the US and UCL velcs to be 1.3 and 1.9km/s
    US_group_std = 1.3 * args.group_velocity_sigma #From Preibisch et al., 2008
    UCL_group_std = 1.3 * args.group_velocity_sigma
    standard_std = {'F':1.08, 'G':0.63, 'K':1.43, 'M':2.27}# 2.0
    astrophysical_std = args.astrophysical_std #Astrophysical radial velocity uncertainty

    Object = []
    Region = []
    IR_excess = []
    Temp_sptype = []
    Pref_template = []
    Obs_info = []
    all_bayes = [[],[]]

    RV_standard_info = {}
    
    sys.stdout.flush()
    CW.Barrier()
    
    #Read in RV standard list
    header = 0
    with open('/home/100/rlk100/RV_standard_list.csv', 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            if header != 0:
                RV_standard_info[row[0]] = (float(row[5]), float(row[6]), float(row[7]))
            else:
                header = 1
        f.close()

    sys.stdout.flush()
    CW.Barrier()

    print "Reading in current spreadsheet", args.input_file
    header = 0
    with open(args.input_file, 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            if header != 0:
                if 'U4' in row[0]:
                    row[0] = 'UCAC4' + row[0].split('U4')[-1]
                Object.append(row[0])
                Region.append(row[1])
                IR_excess.append(row[5])
                Pref_template.append(row[14])
                Temp_sptype.append(row[15])
                if len(row) > 16:
                    Obs = np.array(row[16:])
                    Obs = np.delete(Obs, np.where(Obs==''))
                    #if len(Obs) > 5:
                    #    Obs = np.reshape(Obs, (len(Obs)/5, 5))
                    Obs = np.reshape(Obs, (len(Obs)/6, 6))
                    for ind_obs in Obs:
                        if '/' in ind_obs[0]:
                            new_format = '20' + ind_obs[0].split('/')[-1] + '-' + ind_obs[0].split('/')[-2] + '-' + ("%02d" % int(ind_obs[0].split('/')[-3]))
                            ind_obs[0] = new_format
                else:
                    Obs = np.array([])
                Obs_info.append(Obs)
            if header == 0:
                header = 1
        f.close()
    del header

    sys.stdout.flush()
    CW.Barrier()

    Obj_bayes = np.nan*np.zeros(len(Object))

    #Read in currently calculated Bayes Factors:
    if args.restart_calc != 'False':
        print "Reading in calulated Bayes factors"
        header = 0
        with open(args.bayes_file, 'rU') as f:
            reader = csv.reader(f)
            for row in reader:
                if header != 0:
                    ind = Object.index(row[0])
                    Obj_bayes[ind] = float(row[2])
                    if row[1] == 'US':
                        all_bayes[0].append(float(row[2]))
                    else:
                        all_bayes[1].append(float(row[2]))
                    del ind
                else:
                    header = 1
            f.close()
        del header

    sys.stdout.flush()
    CW.Barrier()

    if args.restart_calc != 'False' and rank == 0:
        print "Creating new bayes file"
        f = open(args.bayes_file, 'w')
        f.write('Object,Region,Bayes_factor\n')
        f.close()

    sys.stdout.flush()
    CW.Barrier()

    inds = range(len(Object))
    skip_inds = np.where(np.array(IR_excess)=='NN')[0]
    for skit in skip_inds:
        inds.remove(skit)
    skip_inds = np.where(np.array(Pref_template)=='')[0]
    for skit in skip_inds:
        inds.remove(skit)
    del skip_inds
    del IR_excess

    rit = 0
    sys.stdout.flush()
    CW.Barrier()
    for obj in inds:

        Pref_template_name = Pref_template[obj].split('_')[0]
        if np.isnan(Obj_bayes[obj]) and rank == rit:
            print "Doing object:", Object[obj], "on rank:", rank            
            likelihoods = []
            single_likelihoods = []
        
            #Produces masses within +/- 10% of the mass of the template.
            #!!! Mike suggests a single mass.
            M_1 = (np.random.random(n_systems)*(RV_standard_info[Pref_template_name][1]-RV_standard_info[Pref_template_name][0])) + RV_standard_info[Pref_template_name][0]
            
            #Generates mass ratios with minium mass ratio of q_min (default 0.01?, should this be dependant on the primary mass? Because sometimes low mass ratios could give very low mass companions i.e. BD mass...)
            #!!! Mike suggests 0.05 due to brown dwarf desert.
            q = (np.random.random(n_systems)*(1-q_min)) + q_min
            
            #from Primary masses and mass ratios, secondary masses can get calculated
            M_2 = M_1 * q
            
            #Get dates of the observations of the object
            jds = Obs_info[obj][:,1].astype(np.float)
            
            #get observed data, and add in the error in the standards in quadrature.
            #This relates to the spectrograph stability
            #There is also an astrophysical error due to these objects being rapid rotators etc.
            RV_standard_err = standard_std[Temp_sptype[obj][0]]
            err = np.sqrt(Obs_info[obj][:,3].astype(float)**2 + RV_standard_err**2 + astrophysical_std**2)
            observed_rv = Obs_info[obj][:,2].astype(float)
            
            #IN A LOOP iterate over random orbits:
            for orb in range(n_orb):
                #FIXME: Figure out which velocity to use!
                if Region[obj] == 'US':
                    if args.group_velocity == 'True':
                        v_group = np.random.normal(US_group_vel, np.sqrt(US_group_std**2 + RV_standard_err**2), n_systems)
                    else:
                        v_group = np.random.normal(np.mean(observed_rv), np.sqrt(US_group_std**2 + RV_standard_err**2), n_systems)
                else:
                    if args.group_velocity == 'True':
                        v_group = np.random.normal(UCL_group_vel, np.sqrt(UCL_group_std**2 + RV_standard_err**2), n_systems)
                    else:
                        v_group = np.random.normal(np.mean(observed_rv), np.sqrt(UCL_group_std**2 + RV_standard_err**2), n_systems)
                

                #generate orbit?
                #!!! Find just one set of orbital parameters at at a time, and
                #scale the RVS. OR if you really want you can compute a, i etc
                #yourself and plug these into my_orb, but some RV scalign is still needed.
                rho, theta, normalised_vr = bo.binary_orbit(my_orb, jds, plot_orbit_no=orb)
                for system in range(n_systems):
                    actual_vr = bo.scale_rv(normalised_vr, my_orb['P'][orb], M_1[system], M_2[system], my_orb['i'][orb], group_velocity=v_group[system])
                    
                    this_likelihood = bo.calc_likelihood(actual_vr, observed_rv, err)
                    likelihoods.append(this_likelihood)
                    #THEN CALCULATE PROBABILITY OF BEING A SINGLE STAR
                    single_likelihoods.append(bo.calc_likelihood(v_group[system], observed_rv, err))
                    del actual_vr
                    del this_likelihood
                del v_group
            del M_1
            del q
            del M_2
            del jds
            del RV_standard_err
            del err
            del observed_rv
            
            #THEN CALCULATE BAYES FACTOR
            bayes_factor = np.mean(likelihoods)/np.mean(single_likelihoods)
            print("Bayes Factor: {0:5.2f} for ".format(bayes_factor) + Object[obj]), "on rank", rank, "with SpT", Temp_sptype[obj]
            del likelihoods
            del single_likelihoods
            if Region[obj] == 'US':
                send_data = [0.0,float(obj),bayes_factor, Temp_sptype[obj]]
                #print "Sending data:", send_data, "from rank:", rank
                if rank == 0:
                    bayes_update = send_data
                else:
                    CW.send(send_data, dest=0, tag=rank)
            else:
                send_data = [1.0,float(obj),bayes_factor, Temp_sptype[obj]]
                #print "Sending data:", send_data, "from rank:", rank
                if rank == 0:
                    bayes_update = send_data
                else:
                    CW.send(send_data, dest=0, tag=rank)
            del send_data
            if rank == 0:
                all_bayes[int(bayes_update[0])].append(bayes_update[2])
                Obj_bayes[int(bayes_update[1])] = bayes_update[2]
                print "Updated Bayes factors retrieved from rank 0 for object", Object[int(bayes_update[1])]
                f = open(args.bayes_file, 'a')
                write_string = Object[int(bayes_update[1])] + ',' + Region[int(bayes_update[1])] + ',' + str(bayes_update[2]) + ',' + str(bayes_update[3]) + '\n'
                f.write(write_string)
                f.close()
                del bayes_update
                del write_string
                    
    
        rit = rit +1
        if rit == size:
            sys.stdout.flush()
            CW.Barrier()
            rit = 0
            if rank == 0:
                
                print "UPDATING CALCULATED BAYES VALUES"
                for orit in range(1,size):
                    bayes_update = CW.recv(source=orit, tag=orit)
                    all_bayes[int(bayes_update[0])].append(bayes_update[2])
                    Obj_bayes[int(bayes_update[1])] = bayes_update[2]
                    print "Updated Bayes factors retrieved from rank", orit, "for object", Object[int(bayes_update[1])]
                    f = open(args.bayes_file, 'a')
                    write_string = Object[int(bayes_update[1])] + ',' + Region[int(bayes_update[1])] + ',' + str(bayes_update[2]) + ',' + str(bayes_update[3]) + '\n'
                    f.write(write_string)
                    f.close()
                    del bayes_update
                    del write_string
            sys.stdout.flush()
            CW.Barrier()

    sys.stdout.flush()
    CW.Barrier()
    if rank == 0:
        
        print "UPDATING CALCULATED BAYES VALUES"
        for orit in range(1,size):
            bayes_update = CW.recv(source=orit, tag=orit)
            all_bayes[int(bayes_update[0])].append(bayes_update[2])
            Obj_bayes[int(bayes_update[1])] = bayes_update[2]
            print "Updated Bayes factors retrieved from rank", orit, "for object", Object[int(bayes_update[1])]
            f = open(args.bayes_file, 'a')
            write_string = Object[int(bayes_update[1])] + ',' + Region[int(bayes_update[1])] + ',' + str(bayes_update[2]) + ',' + str(bayes_update[3]) + '\n'
            f.write(write_string)
            f.close()
            del bayes_update
            del write_string
        sys.stdout.flush()
        CW.Barrier()
    print "Finished Calculating bayes factors!"

if __name__ == '__main__': main()
