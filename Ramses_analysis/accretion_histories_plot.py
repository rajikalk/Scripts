import numpy as np
import matplotlib.pyplot as plt
import yt
import matplotlib.pyplot as plt
import pickle

particle_pickle = "/g/data/ek9/rlk100/RAMSES/Global/stars_red_512.pkl"
file = open(particle_pickle, 'rb')
particle_data = pickle.load(file)
file.close()

#Plot sink histories for classification
for sink_id in range(len(particle_data['dm'][0])):
    plt.clf()
    plt.semilogy(particle_data['time'].T[sink_id], particle_data['dm'].T[sink_id])
    plt.xlabel("time")
    plt.ylabel("Accretion rate")
    plt.title("Sink ID = "+str(sink_id))
    plt.savefig("Accretion_history_"+str(sink_id)+".png")
    print("saved accretion history for sink", sink_id)
    
