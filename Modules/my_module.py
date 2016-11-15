#!/usr/bin/env python
import yt
from yt.utilities.exceptions import YTOutputNotIdentified
import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob

def find_sink_formation_time(files):
    for source in range(len(files)):
        f = h5py.File(files[source], 'r')
        if 'particlepositions' in f.keys():
            particles = True
            sink_form = f['time'][0]/yt.units.yr.in_units('s').value
            break
        else:
            sink_form = f['time'][0]/yt.units.yr.in_units('s').value
    return sink_form

def generate_frame_times(files, dt, start_time=None, presink_frames=25):
    try:
        file = files[-1]
        part_file=file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        sink_form_time = np.min(dd['particle_creation_time'].value/yt.units.yr.in_units('s').value)
        max_time = ds.current_time.in_units('yr').value - sink_form_time
    except YTOutputNotIdentified:
        f = h5py.File(files[-1], 'r')
        sink_form_time = find_sink_formation_time(files)
        max_time = f['time'][0]/yt.units.yr.in_units('s').value - sink_form_time
        f.close()

    if start_time == None:
        m_times = np.logspace(0.0, np.log10(sink_form_time), presink_frames) - sink_form_time
        m_times = m_times.tolist()
    elif start_time < 0.0:
        m_times = np.logspace(np.log10(start_time), np.log10(sink_form_time), presink_frames) - sink_form_time
        m_times = m_times.tolist()
    else:
        m_times = []

    postsink = 0.0
    while postsink < (max_time):
        if start_time != None:
            if postsink >= start_time:
                m_times.append(postsink)
        else:
            m_times.append(postsink)
        postsink = postsink + dt
    return m_times

def find_files(m_times, files):
    try:
        file = files[-1]
        part_file=file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        sink_form_time = np.min(dd['particle_creation_time'].value/yt.units.yr.in_units('s').value)
        yt_file = True
    except YTOutputNotIdentified:
        f = h5py.File(files[-1], 'r')
        sink_form_time = find_sink_formation_time(files)
        f.close()
        yt_file = False
    usable_files = []
    mit = 0
    min = 0
    max = len(files)-1
    pit = 0
    while mit < len(m_times):
        it = int(np.round(min + ((max - min)/2.)))
        #print 'search iterator =', it
        if yt_file:
            file = files[it]
            part_file=file[:-12] + 'part' + file[-5:]
            ds = yt.load(file, particle_filename=part_file)
            time = ds.current_time.in_units('yr').value-sink_form_time
        else:
            f = h5py.File(files[it], 'r')
            time = f['time'][0]/yt.units.yr.in_units('s').value-sink_form_time
        if pit == it or time == m_times[mit]:
            if m_times[mit] == 0.0:
                if time < 0.0:
                    it = it + 1
                    while time < 0.0:
                        if yt_file:
                            file = files[it]
                            part_file=file[:-12] + 'part' + file[-5:]
                            ds = yt.load(file, particle_filename=part_file)
                            time = ds.current_time.in_units('yr').value-sink_form_time
                        else:
                            f = h5py.File(files[it], 'r')
                            time = f['time'][0]/yt.units.yr.in_units('s').value-sink_form_time
                        if time > 0.0:
                            usable_files.append(file)
                            print "found time", time, "for m_time", m_times[mit]
                        else:
                            it = it + 1
                else:
                    usable_files.append(file)
                    print "found time", time, "for m_time", m_times[mit]
            else:
                usable_files.append(file)
                print "found time", time, "for m_time", m_times[mit]
            mit = mit + 1
            min = it
            max = len(files)-1
            pit = it
        elif time > m_times[mit]:
            max = it
            pit = it
        elif time < m_times[mit]:
            min = it
            pit = it
    return usable_files


