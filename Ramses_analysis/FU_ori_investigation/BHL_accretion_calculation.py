#!/usr/bin/env python
"""
Estimate the Bondi-Hoyle-Lyttleton (BHL) accretion rate onto a sink from the gas
inside a measuring sphere, for every RAMSES output of one burst event.

Pickle format (format_version 2), used for both the per-rank pickles
(BHL_accretion_<rank>.pkl) and the merged pickle (BHL_accretion.pkl):
    "format_version"     : 2
    "Estimators"         : ["mean", "median", "mass_weighted"]  (column order below)
    "Alpha"              : [1, 2]  (BHL_Acc_acc_low uses alpha=1, _high uses alpha=2)
    "Measuring_radius_au": float
    "Rel_vel_method"     : "cell" or "bulk"
    "File"               : list of N output paths, relative to the event's data/ dir
    "Time"               : (N,)   time since sink formation [yr]
    "BHL_Acc_acc_low"    : (N, 3) alpha * 2 pi (G M)^2 rho / (v^2 + cs^2)^1.5 with alpha=1 [Msun/yr]
    "BHL_Acc_acc_high"   : (N, 3) same with alpha=2 [Msun/yr]
    "Bondi_radius"       : (N, 3) 2 G M / (cs^2 + v^2) [au]
    "Rel_speed"          : (N, 3) speed of the gas in the measuring sphere relative to the sink, v in the
                           formulas above (definition set by Rel_vel_method) [km/s]
    "Density"            : list of N arrays, per-cell density in the sphere [g/cm^3]
    "Rel_kep"            : list of N arrays, per-cell v_tan / v_kep in the sphere
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
ARRAY_KEYS = ["BHL_Acc_acc_low", "BHL_Acc_acc_high", "Bondi_radius", "Rel_speed"]   # (N, 3)
LIST_KEYS = ["Density", "Rel_kep"]                                      # ragged, one array per snapshot
META_KEYS = ["Measuring_radius_au", "Rel_vel_method"]

parser = argparse.ArgumentParser()
parser.add_argument("-event_id", "--event_identifier", default=2, type=int)
parser.add_argument("-ax", "--axis", default='xy', type=str)
parser.add_argument("-sph_rad", "--measuring_sphere_radius", default=5, type=float)
parser.add_argument("-vrel", "--rel_vel_method", default='cell', choices=['cell', 'bulk'],
                    help="'cell': estimators of the per-cell speed |v_gas - v_sink| (original behaviour, includes disc rotation). "
                         "'bulk': speed of the sphere's bulk (mean/median/mass-weighted) velocity relative to the sink.")
parser.add_argument("-est", "--plot_estimator", default='mass_weighted', choices=ESTIMATORS,
                    help="Which estimator to show in the plot")
parser.add_argument('files', nargs='*')
args = parser.parse_args()

#------------------------------------------------------
# Pickle helpers
def empty_bhl_dict():
    d = {"format_version": 2, "Estimators": list(ESTIMATORS), "Alpha": [1, 2],
         "Measuring_radius_au": float(args.measuring_sphere_radius),
         "Rel_vel_method": args.rel_vel_method,
         "File": [], "Time": []}
    for key in ARRAY_KEYS + LIST_KEYS:
        d[key] = []
    return d

def load_bhl_pickle(path):
    """Return the dict in path if it is in the current format, else None."""
    with open(path, 'rb') as file_open:
        d = pickle.load(file_open)
    if isinstance(d, dict) and d.get("format_version") == 2:
        return d
    if rank == 0:
        print("WARNING: ignoring", path, "- old pickle format (its BHL arrays were mis-sorted), it will be recomputed/overwritten")
    return None

def check_meta(d, path):
    for key in META_KEYS:
        if d[key] != empty_bhl_dict()[key]:
            raise SystemExit("ERROR: " + path + " was made with " + key + "=" + str(d[key]) +
                             " but this run uses " + str(empty_bhl_dict()[key]) +
                             ". Run in a separate directory or remove the old pickles.")

def merge_bhl_dicts(dicts):
    """Combine dicts snapshot-by-snapshot, unique by File (later dicts win), sorted by Time."""
    merged = empty_bhl_dict()
    for d in dicts:
        if d is not None:
            for key in META_KEYS:
                merged[key] = d[key]
    entries = {}
    for d in dicts:
        if d is None:
            continue
        for i, fname in enumerate(d["File"]):
            entries[fname] = (d, i)
    order = sorted(entries, key=lambda f: entries[f][0]["Time"][entries[f][1]])
    merged["File"] = order
    merged["Time"] = np.array([entries[f][0]["Time"][entries[f][1]] for f in order], dtype=float)
    for key in ARRAY_KEYS:
        # keys missing from older pickles (e.g. Rel_speed) are filled with NaN
        merged[key] = np.array([entries[f][0][key][entries[f][1]] if key in entries[f][0] else [np.nan]*len(ESTIMATORS)
                                for f in order], dtype=float).reshape(-1, len(ESTIMATORS))
    for key in LIST_KEYS:
        merged[key] = [np.asarray(entries[f][0][key][entries[f][1]]) for f in order]
    return merged

def save_pickle(d, path):
    tmp_path = path + ".tmp"
    with open(tmp_path, 'wb') as file_open:
        pickle.dump(d, file_open)
    os.replace(tmp_path, path)   # atomic, so a crash mid-write can't corrupt the pickle

def estimators(values, weights):
    """[mean, median, mass-weighted mean]"""
    return np.array([np.mean(values), np.median(values), np.sum(values*weights)/np.sum(weights)])

#------------------------------------------------------
time_bounds = [[3800, 4900],[5575, 5700], [6580, 6720], [7295, 7365], [7850, 7900]]

event_it = args.event_identifier
start_time = time_bounds[event_it -1][0]
end_time = time_bounds[event_it -1][1]
if event_it == 4 and os.getcwd().split('/')[-1] == 'End_7340':
    end_time = 7340
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)

sim_data_dir = '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/Event_'+str(event_it)+'/data/'
all_files = sorted(glob.glob(sim_data_dir+"*/info*.txt"))#[::10]

sink_id = 45
sink_form_time_yr = yt.YTQuantity(22926444.19370405, 'yr').value

# Work out which outputs are already done, from every existing pickle
rank_pickle = 'BHL_accretion_'+str(rank)+'.pkl'
existing_pickles = (['BHL_accretion.pkl'] if os.path.exists('BHL_accretion.pkl') else []) + sorted(glob.glob("BHL_accretion_*.pkl"))
done_files = set()
for pickle_file in existing_pickles:
    d = load_bhl_pickle(pickle_file)
    if d is not None:
        check_meta(d, pickle_file)
        done_files.update(d["File"])

# This rank keeps appending to its own pickle (if it is in the current format)
my_dict = empty_bhl_dict()
if os.path.exists(rank_pickle):
    d = load_bhl_pickle(rank_pickle)
    if d is not None:
        for key in ["File", "Time"] + ARRAY_KEYS + LIST_KEYS:
            if key in d:
                my_dict[key] = list(d[key])
            else:   # key added after this pickle was made (e.g. Rel_speed)
                my_dict[key] = [[np.nan]*len(ESTIMATORS)]*len(d["File"])

files = [fn for fn in all_files if os.path.relpath(fn, sim_data_dir) not in done_files]
if rank == 0:
    print(len(all_files) - len(files), "of", len(all_files), "outputs already done,", len(files), "to do")

sys.stdout.flush()
CW.Barrier()
gc.collect()

#------------------------------------------------------
G = float(yt.units.gravitational_constant_cgs.in_cgs().value)

def calculate_bhl(ds, fn):
    au_cm = ds.quan(1, 'au').in_cgs().value
    msun_g = ds.quan(1, 'Msun').in_cgs().value
    yr_s = ds.quan(1, 'yr').in_cgs().value

    time_val = ds.current_time.in_units('yr').value - sink_form_time_yr

    # Sink properties
    sink_pos = np.array([ds.r["gas", "sink_particle_pos"+ax][sink_id].in_cgs().value for ax in 'xyz'])
    sink_vel = np.array([ds.r["gas", "sink_particle_vel"+ax][sink_id].in_cgs().value for ax in 'xyz'])
    sink_mass = ds.r["gas", "sink_particle_mass"][sink_id].in_cgs().value

    # Only read the cells around the sink, instead of the whole domain
    radius_cm = args.measuring_sphere_radius*au_cm
    sph = ds.sphere(ds.arr(sink_pos, 'cm'), ds.quan(radius_cm, 'cm'))
    sep_vector = np.array([sph["index", ax].in_cgs().value - sink_pos[i] for i, ax in enumerate('xyz')])
    radii = np.sqrt(np.sum(sep_vector**2, axis=0))
    in_sph = radii <= radius_cm   # same cell selection as before (cell centre within the radius)
    sep_vector = sep_vector[:, in_sph]
    radii = radii[in_sph]

    gas_mass = sph["gas", "mass"].in_cgs().value[in_sph]
    density = sph["gas", "Density"].in_cgs().value[in_sph]
    gas_vel = np.array([sph["ramses", ax+"-velocity"].in_cgs().value[in_sph] for ax in 'xyz'])
    gamma = np.asarray(sph["gas", "Gamma"])[in_sph]
    pressure = sph["gas", "Pressure"].in_cgs().value[in_sph]
    del sph
    sound_speed = np.sqrt(gamma*pressure/density)

    # Keplerian velocity from sink mass + gas mass enclosed (cells at the same radius count as enclosed)
    sort_inds = np.argsort(radii)
    cumulative_mass = np.cumsum(gas_mass[sort_inds])
    enclosed_mass = cumulative_mass[np.searchsorted(radii[sort_inds], radii, side='right') - 1] + sink_mass
    with np.errstate(divide='ignore', invalid='ignore'):
        keplerian_velocity = np.where(radii > 0, np.sqrt(G*enclosed_mass/np.where(radii > 0, radii, 1.0)), np.nan)

    # Velocity relative to the sink, split into radial and tangential parts
    rel_vel = gas_vel - sink_vel[:, None]
    rel_speed_cells = np.sqrt(np.sum(rel_vel**2, axis=0))
    with np.errstate(divide='ignore', invalid='ignore'):
        rad_speed = np.where(radii > 0, np.sum(rel_vel*sep_vector, axis=0)/np.where(radii > 0, radii, 1.0), 0.0)
        tang_vel = np.sqrt(np.clip(rel_speed_cells**2 - rad_speed**2, 0, None))
        rel_kep = tang_vel/keplerian_velocity

    # Estimators [mean, median, mass-weighted]
    density_est = estimators(density, gas_mass)
    sound_speed_est = estimators(sound_speed, gas_mass)
    if args.rel_vel_method == 'cell':
        rel_speed_est = estimators(rel_speed_cells, gas_mass)
    else:
        bulk_vels = [np.mean(gas_vel, axis=1), np.median(gas_vel, axis=1), np.sum(gas_vel*gas_mass, axis=1)/np.sum(gas_mass)]
        rel_speed_est = np.array([np.sqrt(np.sum((v - sink_vel)**2)) for v in bulk_vels])

    # BHL rate and radius, each using matching estimators
    v2_plus_cs2 = rel_speed_est**2 + sound_speed_est**2
    bhl_alpha_1 = 2*np.pi*(G*sink_mass)**2*density_est/v2_plus_cs2**1.5   # g/s
    bhl_alpha_1 = bhl_alpha_1*yr_s/msun_g                                  # Msun/yr
    bondi_radius = 2*G*sink_mass/v2_plus_cs2/au_cm                         # au

    return {"File": os.path.relpath(fn, sim_data_dir), "Time": time_val,
            "BHL_Acc_acc_low": bhl_alpha_1, "BHL_Acc_acc_high": 2*bhl_alpha_1,
            "Bondi_radius": bondi_radius, "Rel_speed": rel_speed_est/1.e5,   # cm/s -> km/s
            "Density": density, "Rel_kep": rel_kep}

if len(files) > 0:
    for fn in yt.parallel_objects(files, njobs=size):   # one output per rank at a time
        print('Reading file', fn, 'on rank', rank)
        sys.stdout.flush()
        ds = yt.load(fn, units_override=units_override)
        result = calculate_bhl(ds, fn)
        for key in result.keys():
            my_dict[key].append(result[key])
        del ds, result
        gc.collect()

        save_pickle(merge_bhl_dicts([my_dict]), rank_pickle)
        print("RANK "+str(rank)+": CALCULATED BHL FOR FILE", files.index(fn), "OF", len(files))
        sys.stdout.flush()

print('Finished BHL Calculation on rank', rank)
sys.stdout.flush()
CW.Barrier()

#------------------------------------------------------
if rank == 0:
    pickle_files = (['BHL_accretion.pkl'] if os.path.exists('BHL_accretion.pkl') else []) + sorted(glob.glob("BHL_accretion_*.pkl"))
    save_dict = merge_bhl_dicts([load_bhl_pickle(p) for p in pickle_files])
    save_pickle(save_dict, 'BHL_accretion.pkl')
    print("Saved BHL_accretion.pkl with", len(save_dict["Time"]), "outputs")

    radius = args.measuring_sphere_radius

    try:
        sink_pickle = "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L20.pkl"
        file_open = open(sink_pickle, 'rb')
        particle_data, counter, _, _ = pickle.load(file_open)
        file_open.close()
        print("finished reading in pickle")
        sys.stdout.flush()
    except:
        try:
            sink_pickle = "/scratch/ek9/rlk100/RAMSES/Analysis/Event_plots/particle_data_L20.pkl"
            print("read pickle", sink_pickle)
            file_open = open(sink_pickle, 'rb')
            particle_data, counter, _, _ = pickle.load(file_open)
            file_open.close()
            print("finished reading in pickle")
            sys.stdout.flush()
        except:
            sink_pickle = "/home/100/rlk100/gdata/RAMSES/Analysis/Event_plots/particle_data_L20.pkl"
            print("read pickle", sink_pickle)
            file_open = open(sink_pickle, 'rb')
            particle_data, counter, _, _ = pickle.load(file_open)
            file_open.close()
            print("finished reading in pickle")
            sys.stdout.flush()

    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "mathtext.fontset": "stixsans"  # Force math to use a sans-serif look
    })

    #Ploting parameters
    two_col_width = 7.20472 #inches
    single_col_width = 3.50394 #inches
    page_height = 10.62472 #inches
    font_size = 9
    mym.set_global_font_size(font_size)

    est_ind = ESTIMATORS.index(args.plot_estimator)
    bhl_low = save_dict["BHL_Acc_acc_low"][:, est_ind]
    bhl_high = save_dict["BHL_Acc_acc_high"][:, est_ind]

    plt.clf()
    fig, ax = plt.subplots(figsize=(two_col_width, 0.6*two_col_width))
    ax.set_title("Burst event "+str(event_it), y=0.8)
    start_ind = np.argmin(abs(particle_data['time']-start_time))
    end_ind = np.argmin(abs(particle_data['time']-end_time))
    mdot = particle_data['mdot'].T[1][start_ind:end_ind]
    lns1 = ax.semilogy(particle_data['time'][start_ind:end_ind], mdot, color='b', ls='-', label="Accretion rate")
    band = ax.fill_between(save_dict["Time"], bhl_low, bhl_high, facecolor=to_rgba('g', 0.3), edgecolor='g', linewidth=0.5,
                           label="BHL prediction ("+args.plot_estimator.replace('_', '-')+r", $\alpha$=1–2)")

    # y-limits that include both the measured and predicted rates in the window
    in_window = (save_dict["Time"] >= start_time) & (save_dict["Time"] <= end_time)
    y_vals = np.concatenate([mdot, bhl_low[in_window], bhl_high[in_window]])
    y_vals = y_vals[np.isfinite(y_vals) & (y_vals > 0)]
    if len(y_vals) > 0:
        ax.set_ylim([0.8*np.min(y_vals), 1.25*np.max(y_vals)])

    axes_1_twin = ax.twinx()
    lns2 = axes_1_twin.plot(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind], ls='--', color='k', alpha=0.5, label="Separation")

    ax.set_xlabel('Time (yr)', labelpad=-0.2, fontsize=font_size)
    ax.set_ylabel(r'Accretion rate (M$_\odot$/yr)', labelpad=-0.2, fontsize=font_size)
    axes_1_twin.set_ylabel('Separation (au)', fontsize=font_size)
    ax.tick_params(axis='x', which='major', direction='in', color='k', top=True)
    ax.tick_params(axis='y', which='major', direction='in', color='k', right=True)
    ax.tick_params(axis='both', labelsize=font_size, labelfontfamily='sans-serif')
    ax.set_xlim([start_time, end_time])
    handles = [lns1[0], band, lns2[0]]
    ax.legend(handles, [h.get_label() for h in handles], loc='upper left')
    plt.savefig("BHL_Event_"+str(event_it)+"_radius_"+str(radius)+"au.png", format='png', bbox_inches='tight', pad_inches=0.02, dpi=300)
    print('Saved figure with BHL Accretion')
    sys.stdout.flush()
