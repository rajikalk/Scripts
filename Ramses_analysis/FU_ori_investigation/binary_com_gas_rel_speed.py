#!/usr/bin/env python
"""
For every RAMSES output of the Sink 45 Level 20 zoom-in:
  1. find the sink nearest to sink 45 and the centre of mass (CoM) of the pair,
  2. compute the CoM velocity from the two sinks' masses and velocities,
  3. take the gas in a sphere (default 500 au) centred on the CoM and compute its bulk
     velocity with the mean, median and mass-weighted estimators,
  4. save the speed of that bulk gas velocity relative to the CoM velocity.

Run with MPI, e.g.  mpirun -np 48 python binary_com_gas_rel_speed.py

Pickle format (CoM_gas_rel_speed.pkl, and per-rank CoM_gas_rel_speed_<rank>.pkl):
    "format_version"     : 1
    "Estimators"         : ["mean", "median", "mass_weighted"]  (column order of Rel_speed)
    "Measuring_radius_au": float
    "File"               : list of N output paths, relative to the data directory
    "Time"               : (N,)   time since sink 45 formed [yr]
    "Rel_speed"          : (N, 3) |v_gas,bulk - v_CoM| [km/s]
    "Companion_id"       : (N,)   id of the sink nearest to sink 45
Snapshots are sorted by Time and unique by File.
"""
import glob
import sys
import argparse
import numpy as np
import pickle
import my_ramses_module as mym
import os
import yt
yt.enable_parallelism()
import my_ramses_fields_short as myf
import gc
from mpi4py.MPI import COMM_WORLD as CW

#-----------------------------------------------------
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

ESTIMATORS = ["mean", "median", "mass_weighted"]
SAVE_NAME = "CoM_gas_rel_speed"

parser = argparse.ArgumentParser()
parser.add_argument("-sph_rad", "--measuring_sphere_radius", default=500, type=float, help="Measuring sphere radius [au]")
parser.add_argument("-data_dir", "--data_directory", default='/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/', type=str)
parser.add_argument('files', nargs='*')
args = parser.parse_args()

#------------------------------------------------------
# Pickle helpers
def empty_dict():
    return {"format_version": 1, "Estimators": list(ESTIMATORS),
            "Measuring_radius_au": float(args.measuring_sphere_radius),
            "File": [], "Time": [], "Rel_speed": [], "Companion_id": []}

def load_pickle(path):
    with open(path, 'rb') as file_open:
        d = pickle.load(file_open)
    if d["Measuring_radius_au"] != float(args.measuring_sphere_radius):
        raise SystemExit("ERROR: " + path + " was made with a " + str(d["Measuring_radius_au"]) + " au sphere, but this run uses " +
                         str(args.measuring_sphere_radius) + " au. Run in a separate directory or remove the old pickles.")
    return d

def merge_dicts(dicts):
    """Combine dicts snapshot-by-snapshot, unique by File (later dicts win), sorted by Time."""
    merged = empty_dict()
    entries = {}
    for d in dicts:
        for i, fname in enumerate(d["File"]):
            entries[fname] = (d, i)
    order = sorted(entries, key=lambda f: entries[f][0]["Time"][entries[f][1]])
    merged["File"] = order
    merged["Time"] = np.array([entries[f][0]["Time"][entries[f][1]] for f in order], dtype=float)
    merged["Rel_speed"] = np.array([entries[f][0]["Rel_speed"][entries[f][1]] for f in order], dtype=float).reshape(-1, len(ESTIMATORS))
    merged["Companion_id"] = np.array([entries[f][0]["Companion_id"][entries[f][1]] for f in order], dtype=int)
    return merged

def save_pickle(d, path):
    tmp_path = path + ".tmp"
    with open(tmp_path, 'wb') as file_open:
        pickle.dump(d, file_open)
    os.replace(tmp_path, path)   # atomic, so a crash mid-write can't corrupt the pickle

def existing_pickles():
    return ([SAVE_NAME+'.pkl'] if os.path.exists(SAVE_NAME+'.pkl') else []) + sorted(glob.glob(SAVE_NAME+"_*.pkl"))

#------------------------------------------------------
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)

sink_id = 45
sink_form_time_yr = yt.YTQuantity(22926444.19370405, 'yr').value

# Outputs either directly under data_dir/data/ or in event subdirectories (data_dir/Event_*/data/)
data_dir = os.path.join(args.data_directory, '')
all_files = sorted(set(glob.glob(data_dir+"data/*/info*.txt") + glob.glob(data_dir+"*/data/*/info*.txt")))
if rank == 0:
    print("Found", len(all_files), "outputs under", data_dir)
    if len(all_files) == 0:
        print("WARNING: no info*.txt files found - check -data_dir")

# Skip outputs that are already done
done_files = set()
for pickle_file in existing_pickles():
    done_files.update(load_pickle(pickle_file)["File"])

rank_pickle = SAVE_NAME+'_'+str(rank)+'.pkl'
my_dict = empty_dict()
if os.path.exists(rank_pickle):
    d = load_pickle(rank_pickle)
    for key in ["File", "Time", "Rel_speed", "Companion_id"]:
        my_dict[key] = list(d[key])

files = [fn for fn in all_files if os.path.relpath(fn, data_dir) not in done_files]
if rank == 0:
    print(len(all_files) - len(files), "of", len(all_files), "outputs already done,", len(files), "to do")

sys.stdout.flush()
CW.Barrier()

#------------------------------------------------------
def calculate_rel_speed(ds, fn):
    au_cm = ds.quan(1, 'au').in_cgs().value
    time_val = ds.current_time.in_units('yr').value - sink_form_time_yr

    # All sink properties, shape (3, N_sinks) for positions and velocities
    sink_pos = np.array([ds.r["gas", "sink_particle_pos"+ax].in_cgs().value for ax in 'xyz'])
    sink_vel = np.array([ds.r["gas", "sink_particle_vel"+ax].in_cgs().value for ax in 'xyz'])
    sink_mass = ds.r["gas", "sink_particle_mass"].in_cgs().value

    # Nearest sink to sink 45
    separations = np.sqrt(np.sum((sink_pos - sink_pos[:, [sink_id]])**2, axis=0))
    separations[sink_id] = np.inf
    companion_id = int(np.argmin(separations))

    # Centre of mass of the pair, from the sink properties
    pair = [sink_id, companion_id]
    total_mass = np.sum(sink_mass[pair])
    com_pos = np.sum(sink_pos[:, pair]*sink_mass[pair], axis=1)/total_mass
    com_vel = np.sum(sink_vel[:, pair]*sink_mass[pair], axis=1)/total_mass

    # Gas in the measuring sphere around the CoM (cell centre within the radius)
    radius_cm = args.measuring_sphere_radius*au_cm
    sph = ds.sphere(ds.arr(com_pos, 'cm'), ds.quan(radius_cm, 'cm'))
    sep_vector = np.array([sph["index", ax].in_cgs().value - com_pos[i] for i, ax in enumerate('xyz')])
    in_sph = np.sqrt(np.sum(sep_vector**2, axis=0)) <= radius_cm
    gas_mass = sph["gas", "mass"].in_cgs().value[in_sph]
    gas_vel = np.array([sph["ramses", ax+"-velocity"].in_cgs().value[in_sph] for ax in 'xyz'])
    del sph, sep_vector

    # Bulk gas velocity [mean, median, mass-weighted], relative to the CoM, as a speed
    bulk_vels = [np.mean(gas_vel, axis=1), np.median(gas_vel, axis=1), np.sum(gas_vel*gas_mass, axis=1)/np.sum(gas_mass)]
    rel_speed = np.array([np.sqrt(np.sum((v - com_vel)**2)) for v in bulk_vels])/1.e5   # km/s

    return {"File": os.path.relpath(fn, data_dir), "Time": time_val, "Rel_speed": rel_speed, "Companion_id": companion_id}

if len(files) > 0:
    for fn in yt.parallel_objects(files, njobs=size):   # one output per rank at a time
        print('Reading file', fn, 'on rank', rank)
        sys.stdout.flush()
        ds = yt.load(fn, units_override=units_override)
        result = calculate_rel_speed(ds, fn)
        for key in result.keys():
            my_dict[key].append(result[key])
        del ds, result
        gc.collect()

        save_pickle(merge_dicts([my_dict]), rank_pickle)
        print("RANK "+str(rank)+": DONE FILE", files.index(fn), "OF", len(files))
        sys.stdout.flush()

print('Finished on rank', rank)
sys.stdout.flush()
CW.Barrier()

#------------------------------------------------------
if rank == 0:
    save_dict = merge_dicts([load_pickle(p) for p in existing_pickles()])
    save_pickle(save_dict, SAVE_NAME+'.pkl')
    print("Saved", SAVE_NAME+'.pkl', "with", len(save_dict["Time"]), "outputs")
    companions = np.unique(save_dict["Companion_id"])
    print("Nearest sink(s) to sink", sink_id, "over the run:", companions)
    sys.stdout.flush()
