import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

text_font = 14
idl_file = 'iras2a_mom1_1.idl'
obs_rv = scipy.io.readsav('iras2a_mom1_1.idl')
image = obs_rv['mom1']
contour_arr = obs_rv['mom0']
con_max_it = int(np.max(obs_rv['mom0'])/obs_rv['rms'])+5 + 5
contour_levels = np.arange(10, con_max_it, 5) * obs_rv['rms']
x = obs_rv['xaxis']
y = obs_rv['yaxis']
X, Y = np.meshgrid(x, y)
xlim = [x[0], x[-1]]
ylim = [y[0], y[-1]]
vla_pos = [obs_rv['vla1xpos'], obs_rv['vla1ypos']]
angle = 113.9
line_length = 100
line_pos_x = [-1*line_length*np.sin(np.radians(angle))+vla_pos[0], line_length*np.sin(np.radians(angle))+vla_pos[0]]
line_pos_y = [-1*line_length*np.cos(np.radians(angle))+vla_pos[1], line_length*np.cos(np.radians(angle))+vla_pos[1]]
title = "HDCO"


plt.clf()
fig, ax = plt.subplots()
ax.set_xlabel('$x$ (AU)', labelpad=-1, fontsize=text_font)
ax.set_ylabel('$y$ (AU)', fontsize=text_font) #, labelpad=-20
ax.set_xlim(xlim)
ax.set_ylim(ylim)

v_cbar_min = 6 #center_vel_rv.in_units('km/s').value - 5
v_cbar_max = 8#center_vel_rv.in_units('km/s').value + 5
#plot = ax.pcolormesh(X, Y, vel_rad/100000, cmap=plt.cm.seismic_r, rasterized=True, vmin=v_cbar_min, vmax=v_cbar_max)
plot = ax.pcolormesh(X, Y, image, cmap='RdYlBu_r', rasterized=True, vmin=v_cbar_min, vmax=v_cbar_max)
ax.contour(X,Y,contour_arr, locator=plt.LogLocator(), linewidths=0.5, colors='k', levels=contour_levels)

plt.gca().set_aspect('equal')
cbar = plt.colorbar(plot, pad=0.0)

plt.plot(line_pos_x, line_pos_y, 'm--', lw=2)

label_string = 'Radial Velocity ($km/s$)'
cbar.set_label(r"{}".format(label_string), rotation=270, labelpad=14, size=text_font)

if len(title) > 0:
    title_text = ax.text((np.mean(xlim)), (ylim[1]-0.05*(ylim[1]-ylim[0])), title, va="center", ha="center", color='k', fontsize=(text_font+4))

plt.tick_params(axis='both', which='major')# labelsize=16)
for line in ax.xaxis.get_ticklines():
    line.set_color('white')
for line in ax.yaxis.get_ticklines():
    line.set_color('white')

plt.savefig("rv_obs.jpg", format='jpg', bbox_inches='tight')
print('saved figure')
#plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
#print('Created frame of radial velocity', proj_number, 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
