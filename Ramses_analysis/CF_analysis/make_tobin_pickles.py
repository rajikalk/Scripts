import numpy as np
import matplotlib.pyplot as plt

def CF_err(CF, N_sys, z=1):
    Err_lower = (1/(1+((z**2)/N_sys)))*(CF + (z**2/(2*N_sys))) - (z/(1+((z**2)/N_sys)))*np.sqrt((((CF*(1-CF))/N_sys)+((z**2)/(4*N_sys))))
    Err_upper = (1/(1+((z**2)/N_sys)))*(CF + (z**2/(2*N_sys))) + (z/(1+((z**2)/N_sys)))*np.sqrt((((CF*(1-CF))/N_sys)+((z**2)/(4*N_sys))))
    return Err_lower, Err_upper

S_bins = np.logspace(0.75,4,14)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2

#Derive bimodal distribution fits
Perseus_singles = 33
Perseus_separations = [24.1, 25.5, 29.1, 75.3, 79.6, 83.4, 90.1, 93.3, 93.9, 103.7, 117.3, 185.4, 186, 225.4, 263.2, 297.4, 414.5, 549, 572.5, 885, 1560.3, 1820.1, 2176.5, 2431.5, 2782.9, 3197, 3281.8, 3974, 4012.3, 4188.1, 4819, 5937.4, 6101.2, 8305.3, 8665.4, 8868, 9585.9]

CF_per_bin_Tobin_per = []
CF_err_per = []
for bin_it in range(1, len(S_bins)):
    smaller_seps = len(np.argwhere(Perseus_separations < S_bins[bin_it-1]))
    N_b = len(np.argwhere((Perseus_separations < S_bins[bin_it])&(Perseus_separations > S_bins[bin_it-1])))
    larger_seps = len(np.argwhere(Perseus_separations > S_bins[bin_it]))
    N_s = Perseus_singles + smaller_seps + 2*larger_seps
    N_sys = N_s + N_b
    cf = N_b/N_sys
    cf_err = CF_err(cf, N_sys)
    CF_per_bin_Tobin_per.append(cf)
    CF_err_per.append(cf_err)

plt.clf()
plt.bar(bin_centers, CF_per_bin_Tobin_per, yerr=CF_err_per, width=0.25, fill=False, edgecolor='black')
#plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, fill=False, edgecolor='black')
plt.ylabel("Companion Frequency")
plt.xlabel("Log (AU)")
plt.xlim([1,4])
plt.ylim([0, 0.2])
plt.savefig("Tobin_2022_perseus.png")

import pdb
pdb.set_trace()

#2016
Tobin_objects = {
"Per-emb-1 (0)": [],
"Per-emb-3 (0)": [],
"Per-emb-9 (0)": [],
"Per-emb-14 (0)": [],
"Per-emb-15 (0)": [],
"Per-emb-19 (0/I)": [],
"Per-emb-20 (0/I)": [],
"Per-emb-23 (0)": [],
"Per-emb-24 (0/I)": [],
"Per-emb-25 (0/I)": [],
"Per-emb-29 (0)": [],
"Per-emb-30 (0/I)": [],
"Per-emb-31 (0/I)": [],
"L1451-MMS (0)": [],
"Per-emb-34 (I)": [],
"Per-emb-38 (I)": [],
"Per-emb-46 (I)": [],
"Per-emb-47 (I)": [],
"Per-emb-50 (I)": [],
"Per-emb-52 (I)": [],
"Per-emb-53 (I)": [],
"Per-emb-54 (I)": [],
"Per-emb-56 (I)": [],
"Per-emb-57 (I)": [],
"Per-emb-61 (I)": [],
"Per-emb-62 (I)": [],
"Per-emb-63 (I)": [],
"Per-emb-64 (I)": [],
"Per-emb-66 (I)": [],
"IRAS 03363+3207 (I?)": [],
"EDJ2009-263 (Flat)": [],
"IRAS 03295+3050 (II)": [],
"L1455IRS2 (Flat)": [],
"EDJ2009-385 (II)": [],
"EDJ2009-172 (II)": [],
"SVS3 (II)": [],
"EDJ2009-173 (II)": [],
"Per-emb-2": [18.4],
"Per-emb-5": [22.3],
"Per-emb-17": [63.9],
"Per-emb-22": [172.8],
"Per-emb-26+Per-emb-42": [1864.0],
"Per-emb-16+Per-emb-28": [3694.5],
"Per-emb-6+Per-emb-10": [7347.9],
"Per-emb-48": [79.5],
"Per-emb-40": [90.0],
"EDJ2009-183": [235.8],
"L1448IRS1": [327.4],
"Per-emb-35": [438.8],
"EDJ2009-156": [714.6],
"Per-emb-58+Per-emb-65": [6641.9],
"EDJ2009-269": [120.6],
"Per-emb-11": [678.8, 2177.8],
"Per-emb-32+EDJ2009+366": [1395.3, 8419.2],
"Per-emb-8+Per-emb-55": [142.1, 2198.2],
"Per-emb-37+EDJ2009+235+EDJ2009+233": [2427.8, 7752.0],
"B1-bS+Per-emb-41+B1-bN": [3210.1, 4000.8],
"Per-emb-36+Per-emb-27": [71.6, 142.6, 7226.6],
"Per-emb-12+Per-emb-13+IRAS4B'":[420.8, 2450.4, 6840.0],
"Per-emb-44+SVS13A2+SVS13B+SVS13C": [69.0, 1222.2, 3434.4, 7941.5],
"Per-emb-18+Per-emb-21+Per-emb-49": [19.6, 71.9, 3048.0, 6319.1],
"Per-emb-33+L1448IRS3A+L1448NW": [57.7, 60.7, 182.8, 1683.0, 4945.6]
}

#S:T:B:Q:5:6 == 37:15:5:2:2:1

CF_per_bin_Tobin = []
CF_errs = []
S_true = 37

for bin_it in range(1, len(S_bins)):
    N_comps_in_sys = []
    for key in Tobin_objects.keys():
        N_comps = len(Tobin_objects[key]) + 1
        smaller_seps = len(np.argwhere(np.array(Tobin_objects[key]) < S_bins[bin_it-1]))
        larger_seps = len(np.argwhere(np.array(Tobin_objects[key]) > S_bins[bin_it]))
        binaries = len(np.argwhere((np.array(Tobin_objects[key]) < S_bins[bin_it])&(np.array(Tobin_objects[key]) > S_bins[bin_it-1])))
        if N_comps == 2:
            if larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 3:
            if larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1]
            elif binaries == 2:
                N_comps_in_sys = N_comps_in_sys + [3]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 4:
            if larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 2 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 5:
            if larger_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1, 1]
            elif smaller_seps == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 3 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 3 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 6:
            if larger_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1, 1]
            elif binaries == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 3 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 3 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 4 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 4 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
    N_t = len(np.argwhere(np.array(N_comps_in_sys) == 3))
    N_b = len(np.argwhere(np.array(N_comps_in_sys) == 2))
    N_s = len(np.argwhere(np.array(N_comps_in_sys) == 1)) + S_true
    N_sys = N_s + N_b + N_t
    cf = (N_b+2*N_t)/N_sys
    N_comp = np.sum(np.array(N_comps_in_sys) - 1)
    CF_err = CF_err(cf, N_sys)
    #import pdb
    #pdb.set_trace()
    #CF_err = ((N_comp*(1-(N_comp/N_sys)))**0.5)*(1/N_sys)
    #CF_err = (N_comp*(1-(N_comp/N_sys))**0.5)*(1/N_sys)
    
    CF_per_bin_Tobin.append(cf)
    CF_errs.append(CF_err)
    
plt.clf()
plt.bar(bin_centers, CF_per_bin_Tobin, yerr=CF_errs, width=0.25, fill=False, edgecolor='black')
#plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, fill=False, edgecolor='black')
plt.ylabel("Companion Frequency")
plt.xlabel("Log (AU)")
plt.xlim([1,4])
plt.ylim([0, 0.2])
plt.savefig("Tobin_2016_full_sample.png")

#Just class 0s

Tobin_objects = {
"Per-emb-1 (0)": [],
"Per-emb-3 (0)": [],
"Per-emb-9 (0)": [],
"Per-emb-14 (0)": [],
"Per-emb-15 (0)": [],
"Per-emb-19 (0/I)": [],
"Per-emb-20 (0/I)": [],
"Per-emb-23 (0)": [],
"Per-emb-24 (0/I)": [],
"Per-emb-25 (0/I)": [],
"Per-emb-29 (0)": [],
"Per-emb-30 (0/I)": [],
"Per-emb-31 (0/I)": [],
"L1451-MMS (0)": [],
"Per-emb-2": [18.4],
"Per-emb-5": [22.3],
"Per-emb-17": [63.9],
"Per-emb-22": [172.8],
"Per-emb-26+Per-emb-42": [1864.0],
"Per-emb-8+Per-emb-55": [2198.2],
"Per-emb-16+Per-emb-28": [3694.5],
"Per-emb-6+Per-emb-10": [7347.9],
"Per-emb-27+Per-emb-36": [142.6, 7226.6],
"Per-emb-11": [678.8, 2177.8],
"Per-emb-32+EDJ2009+366": [1395.3],#"Per-emb-32+EDJ2009+366": [1395.3, 8419.2],
"Per-emb-37+EDJ2009+235+EDJ2009+233": [],#"Per-emb-37+EDJ2009+235+EDJ2009+233": [2427.8, 7752.0],
"B1-bS+Per-emb-41+B1-bN": [3210.1, 4000.8],
"Per-emb-18+Per-emb-21+Per-emb-49": [19.6, 3048.0, 6319.1],
"Per-emb-13+IRAS4B’+Per-emb-12": [420.8, 2450.4, 6840.0],
"Per-emb-44+SVS13A2+SVS13B+SVS13C": [69.0, 1222.2, 3434.4, 7941.5],
"Per-emb-33+L1448IRS3A+L1448NW": [57.7, 60.7, 182.8, 1683.0, 4945.6]
}

#S:T:B:Q:5:6 == 16:8:5:2:1:1

CF_per_bin_Tobin = []
CF_errs = []
#S_true = 17#14

for bin_it in range(1, len(S_bins)):
    N_comps_in_sys = []
    for key in Tobin_objects.keys():
        N_comps = len(Tobin_objects[key]) + 1
        smaller_seps = len(np.argwhere(np.array(Tobin_objects[key]) < S_bins[bin_it-1]))
        larger_seps = len(np.argwhere(np.array(Tobin_objects[key]) > S_bins[bin_it]))
        binaries = len(np.argwhere((np.array(Tobin_objects[key]) < S_bins[bin_it])&(np.array(Tobin_objects[key]) > S_bins[bin_it-1])))
        if N_comps == 1:
            N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 2:
            if larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 3:
            if larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1]
            elif binaries == 2:
                N_comps_in_sys = N_comps_in_sys + [3]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 4:
            if larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 2 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 5:
            if larger_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1, 1]
            elif smaller_seps == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 3 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 3 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif binaries == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1]
            elif smaller_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 6:
            if larger_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1, 1]
            elif binaries == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 3 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 3 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 4 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 4 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
    N_t = len(np.argwhere(np.array(N_comps_in_sys) == 3))
    N_b = len(np.argwhere(np.array(N_comps_in_sys) == 2))
    N_s = len(np.argwhere(np.array(N_comps_in_sys) == 1))
    N_sys = N_s + N_b + N_t
    cf = (N_b+2*N_t)/N_sys
    N_comp = np.sum(np.array(N_comps_in_sys) - 1)
    CF_err = CF_err(cf, N_sys)
    #import pdb
    #pdb.set_trace()
    CF_err = ((N_comp*(1-(N_comp/N_sys)))**0.5)*(1/N_sys)
    #CF_err = (N_comp*(1-(N_comp/N_sys))**0.5)*(1/N_sys)
    
    CF_per_bin_Tobin.append(cf)
    CF_errs.append(CF_err)
    
plt.clf()
plt.bar(bin_centers, CF_per_bin_Tobin, yerr=CF_errs, width=0.25, fill=False, edgecolor='black')
#plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, fill=False, edgecolor='black')
plt.ylabel("Companion Frequency")
plt.xlabel("Log (AU)")
plt.xlim([1,4])
plt.ylim([0, 0.2])
plt.savefig("Tobin_2016_class_0.png")

#2018
Tobin_objects = {
"Per-emb-1 (0)": [],
"Per-emb-3 (0)": [],
"Per-emb-9 (0)": [],
"Per-emb-14 (0)": [],
"Per-emb-15 (0)": [],
"Per-emb-19 (0/I)": [],
"Per-emb-20 (0/I)": [],
"Per-emb-23 (0)": [],
"Per-emb-24 (0/I)": [],
"Per-emb-25 (0/I)": [],
"Per-emb-29 (0)": [],
"Per-emb-30 (0/I)": [],
"Per-emb-31 (0/I)": [],
"L1451-MMS (0)": [],
"Per-emb-34 (I)": [],
"Per-emb-38 (I)": [],
"Per-emb-46 (I)": [],
"Per-emb-47 (I)": [],
"Per-emb-50 (I)": [],
"Per-emb-52 (I)": [],
"Per-emb-53 (I)": [],
"Per-emb-54 (I)": [],
"Per-emb-56 (I)": [],
"Per-emb-57 (I)": [],
"Per-emb-61 (I)": [],
"Per-emb-62 (I)": [],
"Per-emb-63 (I)": [],
"Per-emb-64 (I)": [],
"Per-emb-66 (I)": [],
"IRAS 03363+3207 (I?)": [],
"EDJ2009-263 (Flat)": [],
"IRAS 03295+3050 (II)": [],
"L1455IRS2 (Flat)": [],
"EDJ2009-385 (II)": [],
"EDJ2009-172 (II)": [],
"SVS3 (II)": [],
"EDJ2009-173 (II)": [],
"EDJ2009+233": [],
"SVS13C": [],
"EDJ2009+366": [],
"Per-emb-2": [24.0],
"Per-emb-5": [29.1],
"Per-emb-17": [83.3],
"Per-emb-22": [225.4],
"Per-emb-26+Per-emb-42": [2431.3],
"Per-emb-16+Per-emb-28": [4818.9],
"Per-emb-6+Per-emb-10": [9584.2],
"Per-emb-48": [103.7],
"Per-emb-40": [117.4],
"EDJ2009-183": [307.6],
"L1448IRS1": [427.0],
"Per-emb-35": [572.3],
"EDJ2009-156": [932.1],
"Per-emb-58+Per-emb-65": [8663.3],
"EDJ2009-269": [157.3],
"Per-emb-11": [885.4, 2840.6],
"Per-emb-32": [1820.0],
"Per-emb-8+Per-emb-55": [185.3, 2867.1],
"Per-emb-37+EDJ2009+235": [3166.7],
"B1-bS+Per-emb-41+B1-bN": [4187.1, 5218.4],
"Per-emb-36+Per-emb-27": [93.4, 186.0, 9426.0],
"Per-emb-12+Per-emb-13+IRAS4B'":[548.9, 3196.2, 8921.7],
"Per-emb-44+SVS13A2+SVS13B": [90.0, 1594.2, 4479.7],
"Per-emb-18+Per-emb-21+Per-emb-49": [25.6, 93.8, 3975.7, 8242.3],
"Per-emb-33+L1448IRS3A+L1448NW": [75.3, 79.2, 238.4, 2195.2, 6450.8]
}

CF_per_bin_Tobin = []
CF_errs = []

for bin_it in range(1, len(S_bins)):
    N_comps_in_sys = []
    for key in Tobin_objects.keys():
        N_comps = len(Tobin_objects[key]) + 1
        smaller_seps = len(np.argwhere(np.array(Tobin_objects[key]) < S_bins[bin_it-1]))
        larger_seps = len(np.argwhere(np.array(Tobin_objects[key]) > S_bins[bin_it]))
        binaries = len(np.argwhere((np.array(Tobin_objects[key]) < S_bins[bin_it])&(np.array(Tobin_objects[key]) > S_bins[bin_it-1])))
        if N_comps == 1:
            N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 2:
            if larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 3:
            if larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1]
            elif binaries == 2:
                N_comps_in_sys = N_comps_in_sys + [3]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 4:
            if larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 2 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 5:
            if larger_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1, 1]
            elif smaller_seps == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 3 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 3 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 6:
            if larger_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1, 1]
            elif binaries == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 3 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 3 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 4 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 4 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
    N_t = len(np.argwhere(np.array(N_comps_in_sys) == 3))
    N_b = len(np.argwhere(np.array(N_comps_in_sys) == 2))
    N_s = len(np.argwhere(np.array(N_comps_in_sys) == 1))
    N_sys = N_s + N_b + N_t
    cf = (N_b+2*N_t)/N_sys
    N_comp = np.sum(np.array(N_comps_in_sys) - 1)
    CF_err = CF_err(cf, N_sys)
    #import pdb
    #pdb.set_trace()
    CF_err = ((N_comp*(1-(N_comp/N_sys)))**0.5)*(1/N_sys)
    #CF_err = (N_comp*(1-(N_comp/N_sys))**0.5)*(1/N_sys)
    
    CF_per_bin_Tobin.append(cf)
    CF_errs.append(CF_err)
    
plt.clf()
plt.bar(bin_centers, CF_per_bin_Tobin, yerr=CF_errs, width=0.25, fill=False, edgecolor='black')
#plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, fill=False, edgecolor='black')
plt.ylabel("Companion Frequency")
plt.xlabel("Log (AU)")
plt.xlim([1,4])
plt.ylim([0, 0.2])
plt.savefig("Tobin_2018_full_sample.png")

Tobin_objects = {
"Per-emb-1 (0)": [],
"Per-emb-3 (0)": [],
"Per-emb-9 (0)": [],
"Per-emb-14 (0)": [],
"Per-emb-15 (0)": [],
"Per-emb-19 (0/I)": [],
"Per-emb-20 (0/I)": [],
"Per-emb-23 (0)": [],
"Per-emb-24 (0/I)": [],
"Per-emb-25 (0/I)": [],
"Per-emb-29 (0)": [],
"Per-emb-30 (0/I)": [],
"Per-emb-31 (0/I)": [],
"L1451-MMS (0)": [],
"SVS13C": [],
"Per-emb-2": [24.0],
"Per-emb-5": [29.1],
"Per-emb-17": [83.3],
"Per-emb-22": [225.4],
"Per-emb-26+Per-emb-42": [2431.3],
"Per-emb-8+Per-emb-55": [2867.2],
"Per-emb-16+Per-emb-28": [4818.9],
"Per-emb-6+Per-emb-10": [9584.2],
"Per-emb-27+Per-emb-36": [186.0, 9426.0],
"Per-emb-11": [885.4, 2840.6],
"Per-emb-32": [1820.0],
"Per-emb-37+EDJ2009+235": [],#"Per-emb-37+EDJ2009+235": [3166.7],
"B1-bS+Per-emb-41+B1-bN": [4187.1, 5218.4],
"Per-emb-18+Per-emb-21+Per-emb-49": [25.6, 3975.7, 8242.3],
"Per-emb-13+IRAS4B’+Per-emb-12": [548.9, 3196.2, 8921.7],
"Per-emb-44+SVS13A2+SVS13B": [90.0, 1594.2, 4479.7],
"Per-emb-33+L1448IRS3A+L1448NW": [75.3, 79.2, 238.4, 2195.2, 6450.8]
}

#S:T:B:Q:5:6 == 16:8:5:2:1:1

CF_per_bin_Tobin = []
CF_errs = []
#S_true = 17#14

for bin_it in range(1, len(S_bins)):
    N_comps_in_sys = []
    for key in Tobin_objects.keys():
        N_comps = len(Tobin_objects[key]) + 1
        smaller_seps = len(np.argwhere(np.array(Tobin_objects[key]) < S_bins[bin_it-1]))
        larger_seps = len(np.argwhere(np.array(Tobin_objects[key]) > S_bins[bin_it]))
        binaries = len(np.argwhere((np.array(Tobin_objects[key]) < S_bins[bin_it])&(np.array(Tobin_objects[key]) > S_bins[bin_it-1])))
        if N_comps == 1:
            N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 2:
            if larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 3:
            if larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1]
            elif binaries == 2:
                N_comps_in_sys = N_comps_in_sys + [3]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 4:
            if larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 2 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 5:
            if larger_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1, 1]
            elif smaller_seps == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 3 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 3 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif binaries == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1]
            elif smaller_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 6:
            if larger_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1, 1]
            elif binaries == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 3 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 3 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 4 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 4 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
    N_t = len(np.argwhere(np.array(N_comps_in_sys) == 3))
    N_b = len(np.argwhere(np.array(N_comps_in_sys) == 2))
    N_s = len(np.argwhere(np.array(N_comps_in_sys) == 1))
    N_sys = N_s + N_b + N_t
    cf = (N_b+2*N_t)/N_sys
    N_comp = np.sum(np.array(N_comps_in_sys) - 1)
    CF_err = CF_err(cf, N_sys)
    #import pdb
    #pdb.set_trace()
    CF_err = ((N_comp*(1-(N_comp/N_sys)))**0.5)*(1/N_sys)
    #CF_err = (N_comp*(1-(N_comp/N_sys))**0.5)*(1/N_sys)
    
    CF_per_bin_Tobin.append(cf)
    CF_errs.append(CF_err)
    
plt.clf()
plt.bar(bin_centers, CF_per_bin_Tobin, yerr=CF_errs, width=0.25, fill=False, edgecolor='black')
#plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, fill=False, edgecolor='black')
plt.ylabel("Companion Frequency")
plt.xlabel("Log (AU)")
plt.xlim([1,4])
plt.ylim([0, 0.2])
plt.savefig("Tobin_2018_class_0.png")


#2018 Class 0+I
Tobin_objects = {
"Per-emb-1 (0)": [],
"Per-emb-3 (0)": [],
"Per-emb-9 (0)": [],
"Per-emb-14 (0)": [],
"Per-emb-15 (0)": [],
"Per-emb-19 (0/I)": [],
"Per-emb-20 (0/I)": [],
"Per-emb-23 (0)": [],
"Per-emb-24 (0/I)": [],
"Per-emb-25 (0/I)": [],
"Per-emb-29 (0)": [],
"Per-emb-30 (0/I)": [],
"Per-emb-31 (0/I)": [],
"L1451-MMS (0)": [],
"Per-emb-34 (I)": [],
"Per-emb-38 (I)": [],
"Per-emb-46 (I)": [],
"Per-emb-47 (I)": [],
"Per-emb-50 (I)": [],
"Per-emb-52 (I)": [],
"Per-emb-53 (I)": [],
"Per-emb-54 (I)": [],
"Per-emb-56 (I)": [],
"Per-emb-57 (I)": [],
"Per-emb-61 (I)": [],
"Per-emb-62 (I)": [],
"Per-emb-63 (I)": [],
"Per-emb-64 (I)": [],
"Per-emb-66 (I)": [],
"IRAS 03363+3207 (I?)": [],
"SVS13C (0)": [],
"Per-emb-2": [24.0],
"Per-emb-5": [29.1],
"Per-emb-17": [83.3],
"Per-emb-22": [225.4],
"Per-emb-26+Per-emb-42": [2431.3],
"Per-emb-8+Per-emb-55": [185.3, 2867.2],
"Per-emb-16+Per-emb-28": [4818.9],
"Per-emb-6+Per-emb-10": [9584.2],
"Per-emb-27+Per-emb-36": [93.4, 186.0, 9426.0],
"Per-emb-11": [885.4, 2840.6],
"Per-emb-32": [1820.0],
"Per-emb-37+EDJ2009+235": [],#"Per-emb-37+EDJ2009+235": [3166.7],
"B1-bS+Per-emb-41+B1-bN": [4187.1, 5218.4],
"Per-emb-18+Per-emb-21+Per-emb-49": [25.6, 93.8, 3975.7, 8242.3],
"Per-emb-13+IRAS4B’+Per-emb-12": [548.9, 3196.2, 8921.7],
"Per-emb-44+SVS13A2+SVS13B": [90.0, 1594.2, 4479.7],
"Per-emb-33+L1448IRS3A+L1448NW": [75.3, 79.2, 238.4, 2195.2, 6450.8],
"Per-emb-48": [103.7],
"Per-emb-40": [117.4],
"EDJ2009-183": [307.6],
"L1448IRS1": [427.0],
"Per-emb-35": [572.3],
"Per-emb-58+Per-emb-65": [8663.3]
}

CF_per_bin_Tobin = []
CF_errs = []
#S_true = 17#14

for bin_it in range(1, len(S_bins)):
    N_comps_in_sys = []
    for key in Tobin_objects.keys():
        N_comps = len(Tobin_objects[key]) + 1
        smaller_seps = len(np.argwhere(np.array(Tobin_objects[key]) < S_bins[bin_it-1]))
        larger_seps = len(np.argwhere(np.array(Tobin_objects[key]) > S_bins[bin_it]))
        binaries = len(np.argwhere((np.array(Tobin_objects[key]) < S_bins[bin_it])&(np.array(Tobin_objects[key]) > S_bins[bin_it-1])))
        if N_comps == 1:
            N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 2:
            if larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1]
        elif N_comps == 3:
            if larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1]
            elif binaries == 2:
                N_comps_in_sys = N_comps_in_sys + [3]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 4:
            if larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 2 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 5:
            if larger_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1]
            elif binaries == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1, 1]
            elif smaller_seps == 1 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 1 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 3 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 3 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif binaries == 2 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1]
            elif smaller_seps == 4:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
        elif N_comps == 6:
            if larger_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1, 1, 1]
            elif binaries == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [2, 2, 1, 1]
            elif smaller_seps == 2 and larger_seps == 3:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1, 1]
            elif smaller_seps == 2 and binaries == 1 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [2, 1, 1]
            elif smaller_seps == 3 and larger_seps == 2:
                N_comps_in_sys = N_comps_in_sys + [1, 1, 1]
            elif smaller_seps == 3 and binaries == 1 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [2, 1]
            elif smaller_seps == 4 and larger_seps == 1:
                N_comps_in_sys = N_comps_in_sys + [1, 1]
            elif smaller_seps == 4 and binaries == 1:
                N_comps_in_sys = N_comps_in_sys + [2]
            elif smaller_seps == 5:
                N_comps_in_sys = N_comps_in_sys + [1]
            else:
                import pdb
                pdb.set_trace()
    N_t = len(np.argwhere(np.array(N_comps_in_sys) == 3))
    N_b = len(np.argwhere(np.array(N_comps_in_sys) == 2))
    N_s = len(np.argwhere(np.array(N_comps_in_sys) == 1))
    N_sys = N_s + N_b + N_t
    cf = (N_b+2*N_t)/N_sys
    N_comp = np.sum(np.array(N_comps_in_sys) - 1)
    CF_err = CF_err(cf, N_sys)
    #import pdb
    #pdb.set_trace()
    CF_err = ((N_comp*(1-(N_comp/N_sys)))**0.5)*(1/N_sys)
    #CF_err = (N_comp*(1-(N_comp/N_sys))**0.5)*(1/N_sys)
    
    CF_per_bin_Tobin.append(cf)
    CF_errs.append(CF_err)
    
plt.clf()
plt.bar(bin_centers, CF_per_bin_Tobin, yerr=CF_errs, width=0.25, fill=False, edgecolor='black')
#plt.bar(bin_centers, CF_per_bin_Tobin, width=0.25, fill=False, edgecolor='black')
plt.ylabel("Companion Frequency")
plt.xlabel("Log (AU)")
plt.xlim([1,4])
plt.ylim([0, 0.2])
plt.savefig("Tobin_2018_class_0_I.png")
