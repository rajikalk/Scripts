#!/usr/bin/env python
"""
Plot the measured accretion rate against the BHL predictions from
BHL_accretion_calculation.py for several measuring-sphere radii and burst events.
Expects the format_version 2 pickles described in BHL_accretion_calculation.py.
"""
import glob
import sys
import argparse
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import os

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "mathtext.fontset": "stixsans"  # Force math to use a sans-serif look
    })

ESTIMATORS = ["mean", "median", "mass_weighted"]
ARRAY_KEYS = ["BHL_Acc_acc_low", "BHL_Acc_acc_high", "Bondi_radius"]
LIST_KEYS = ["Density", "Rel_kep"]

#-----------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-event_id", "--event_identifier", default=2, type=int)
parser.add_argument("-est", "--plot_estimator", default='mass_weighted', choices=ESTIMATORS,
                    help="Which estimator to show")
parser.add_argument('files', nargs='*')
args = parser.parse_args()

#------------------------------------------------------
def load_bhl_pickle(path):
    with open(path, 'rb') as file_open:
        d = pickle.load(file_open)
    if isinstance(d, dict) and d.get("format_version") == 2:
        return d
    print("WARNING: skipping", path, "- old pickle format (BHL arrays were mis-sorted). Rerun BHL_accretion_calculation.py")
    return None

def merge_bhl_dicts(dicts):
    """Combine dicts snapshot-by-snapshot, unique by File (later dicts win), sorted by Time."""
    entries = {}
    for d in dicts:
        if d is None:
            continue
        for i, fname in enumerate(d["File"]):
            entries[fname] = (d, i)
    if len(entries) == 0:
        return None
    order = sorted(entries, key=lambda f: entries[f][0]["Time"][entries[f][1]])
    merged = {"File": order,
              "Time": np.array([entries[f][0]["Time"][entries[f][1]] for f in order], dtype=float)}
    for key in ARRAY_KEYS:
        merged[key] = np.array([entries[f][0][key][entries[f][1]] for f in order], dtype=float).reshape(-1, len(ESTIMATORS))
    for key in LIST_KEYS:
        merged[key] = [np.asarray(entries[f][0][key][entries[f][1]]) for f in order]
    return merged

def load_bhl_directory(directory):
    """Merged pickle plus any per-rank pickles (e.g. from a run still in progress)."""
    pickle_files = ([directory+'BHL_accretion.pkl'] if os.path.exists(directory+'BHL_accretion.pkl') else []) \
                   + sorted(glob.glob(directory+"BHL_accretion_*.pkl"))
    return merge_bhl_dicts([load_bhl_pickle(p) for p in pickle_files])

#------------------------------------------------------
#Ploting parameters
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 9

time_bounds = [[3800, 4900],[5575, 5700], [6580, 6720], [7295, 7365], [7850, 7900]]
#time_bounds = [[3800, 4900+500],[5575, 5700+500], [6580, 6720+500], [7295, 7365+500], [7850, 7900+500]]
colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']
event_IDS = [2, 3, 4, 5]
Radii = [1, 2, 3, 4, 5, 6]
est_ind = ESTIMATORS.index(args.plot_estimator)

try:
    sink_pickle = "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L20.pkl"
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_id, sink_form_time = pickle.load(file_open)
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


# Load the BHL results once, for both figures
bhl_data = {}
for event_it in event_IDS:
    for rad in Radii:
        directory = "/home/100/rlk100/rlk/RAMSES/Analysis/BHL_analytical_calc/Event_"+str(event_it)+"/Radius_"+str(rad) +"/"
        bhl_data[(event_it, rad)] = load_bhl_directory(directory)
        if bhl_data[(event_it, rad)] is None:
            directory = "/home/100/rlk100/gdata/RAMSES/Analysis/BHL_analytical_calc/Event_"+str(event_it)+"/Radius_"+str(rad) +"/"
            bhl_data[(event_it, rad)] = load_bhl_directory(directory)
            if bhl_data[(event_it, rad)] is None:
                print("No data right now for Event", event_it, "Radius", rad)
            

def get_window(event_it):
    start_time = time_bounds[event_it -1][0]
    end_time = time_bounds[event_it -1][1]
    if event_it == 4 and os.getcwd().split('/')[-1] == 'End_7340':
        end_time = 7340
    start_ind = np.argmin(abs(particle_data['time']-start_time))
    end_ind = np.argmin(abs(particle_data['time']-end_time))
    return start_time, end_time, start_ind, end_ind

def style_panel(ax, event_it, start_time, end_time):
    ax.set_xlabel('Time (yr)', labelpad=-0.2, fontsize=font_size)
    ax.tick_params(axis='x', which='major', direction='in', color='k', top=True)
    ax.tick_params(axis='y', which='major', direction='in', color='k', right=True)
    ax.tick_params(axis='both', labelsize=font_size, labelfontfamily='sans-serif')
    ax.set_xlim([start_time, end_time])
    ax.text(0.05, 0.95, "Event "+str(event_it), transform=ax.transAxes, ha='left', va='top', fontsize=font_size)
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, pos.width, pos.height * 0.85])
    if event_it != 2:
        plt.setp(ax.get_yticklabels(), visible=False)

def add_legend(ax, handles, title):
    ax.legend(handles, [h.get_label() for h in handles],
        loc='upper center',
        bbox_to_anchor=(2.0, -0.25),
        ncol=4,
        title=title,
        title_fontsize=font_size,
    )

estimator_label = args.plot_estimator.replace('_', '-')

#------------------------------------------------------
# Figure 1: measured accretion rate vs BHL prediction
fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(two_col_width, 0.6*single_col_width), sharey=True)
plt.subplots_adjust(hspace=-0.61)
plt.subplots_adjust(wspace=0.0)

for event_it in event_IDS:
    ax = axs.flatten()[event_it-2]
    start_time, end_time, start_ind, end_ind = get_window(event_it)
    handles = []
    lns1 = ax.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[1][start_ind:end_ind], color='k', ls='-', label="Accretion rate")
    handles.append(lns1[0])

    for rad in Radii:
        save_dict = bhl_data[(event_it, rad)]
        if save_dict is None:
            continue
        colour = colours[Radii.index(rad)]
        band = ax.fill_between(save_dict["Time"], save_dict["BHL_Acc_acc_low"][:, est_ind], save_dict["BHL_Acc_acc_high"][:, est_ind],
                               facecolor=to_rgba(colour, 0.15), edgecolor=colour, linewidth=0.5, label="r = "+str(rad)+" au")
        handles.append(band)

    axes_1_twin = ax.twinx()
    lns2 = axes_1_twin.plot(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind], ls='--', color='k', alpha=0.5, label="Separation")
    handles.append(lns2[0])

    if event_it == 2:
        ax.set_ylabel(r'Accretion rate (M$_\odot$/yr)', labelpad=-0.2, fontsize=font_size)
        ax.set_ylim([5.e-9, 1.e-4])
    if event_it == 5:
        axes_1_twin.set_ylabel('Separation (au)', fontsize=font_size)
    axes_1_twin.set_ylim([5, 70])   # same separation axis in every panel
    axes_1_twin.tick_params(axis='y', labelsize=font_size)
    if event_it != 5:
        plt.setp(axes_1_twin.get_yticklabels(), visible=False)
    style_panel(ax, event_it, start_time, end_time)
    if event_it == 2:
        add_legend(ax, handles, "BHL prediction ("+estimator_label+r", $\alpha$=1–2)")

fig.savefig("BHL_Event_multi_rad_"+args.plot_estimator+".png", format='png', bbox_inches='tight', pad_inches=0.02, dpi=300)
print('Saved figure with BHL Accretion')
sys.stdout.flush()

#------------------------------------------------------
# Figure 2: BHL (Bondi) radius, with a reference line at 0.79 au
fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(two_col_width, 0.6*single_col_width), sharey=True)
plt.subplots_adjust(hspace=-0.61)
plt.subplots_adjust(wspace=0.0)

for event_it in event_IDS:
    ax = axs.flatten()[event_it-2]
    start_time, end_time, start_ind, end_ind = get_window(event_it)
    handles = []

    for rad in Radii:
        save_dict = bhl_data[(event_it, rad)]
        if save_dict is None:
            continue
        colour = colours[Radii.index(rad)]
        lns = ax.semilogy(save_dict["Time"], save_dict["Bondi_radius"][:, est_ind], color=colour, ls='-', lw=1, label="r = "+str(rad)+" au")
        handles.append(lns[0])

    hline = ax.axhline(0.79, ls='--', color='k', alpha=0.5, lw=1, label="0.79 au")
    handles.append(hline)

    if event_it == 2:
        ax.set_ylabel('Radius (au)', labelpad=-0.2, fontsize=font_size)
    style_panel(ax, event_it, start_time, end_time)
    if event_it == 2:
        add_legend(ax, handles, r"Bondi radius, $2GM/(c_s^2+v^2)$ ("+estimator_label+")")

fig.savefig("Bondi_radius_Event_multi_rad_"+args.plot_estimator+".png", format='png', bbox_inches='tight', pad_inches=0.02, dpi=300)
print('Saved figure with Bondi radius')
sys.stdout.flush()
