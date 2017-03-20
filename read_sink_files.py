import csv
import numpy as np
import sys
import matplotlib.pyplot as plt
import yt

#============================================================================================
file = sys.argv[1]

particle_tags = []
time = []
position = []
velocity = []
acceleration = []
angular_momentum = []
mass = []
accretion = []

csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)

with open(file, 'r') as f:
    reader = csv.reader(f, dialect='dat')
    for row in reader:
        if row[0][0] != '[':
            part_tag = int(row[0])
            if part_tag not in particle_tags:
                particle_tags.append(part_tag)
                time.append([])
                position.append([])
                velocity.append([])
                acceleration.append([])
                angular_momentum.append([])
                mass.append([])
                accretion.append([])
            pit = particle_tags.index(part_tag)

            t = float(row[1]) - float(row[15])
            pos = [float(row[2]), float(row[3]), float(row[4])]
            vel = [float(row[5]), float(row[6]), float(row[7])]
            acc = [float(row[8]), float(row[9]), float(row[10])]
            ang = np.sqrt(float(row[11])**2. + float(row[12])**2. + float(row[13])**2.)
            m = float(row[14])
            mdot = float(row[15])

            time[pit].append(t)
            position[pit].append(pos)
            velocity[pit].append(vel)
            acceleration[pit].append(acc)
            angular_momentum[pit].append(ang)
            mass[pit].append(m)
            accretion.append(mdot)

particle_tags = np.array(particle_tags)
time = np.array(time)
position = np.array(position)
velocity = np.array(velocity)
acceleration = np.array(acceleration)
angular_momentum = np.array(angular_momentum)
mass = np.array(mass)
accretion = np.array(accretion)

for particle in range(len(particle_tag)):
    plt.plot(time[particle]/yt.units.year.in_units('s').value, angular_momentum[particle]/mass[particle])
plt.x



