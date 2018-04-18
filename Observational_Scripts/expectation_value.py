%run /home/100/rlk100/Scripts/bayesian_analysis.py -file /home/100/rlk100/WiFeS_target_spreadsheet_cleaned.csv -suffix '_resubmission' -bayes_f bayes_factors_cleaned.csv
SpT_bin_1 = ['F6', 'F7', 'F8', 'F9', 'G0', 'G1', 'G2']
SpT_bin_2 = ['G3', 'G5', 'G6', 'G7', 'G9', 'K0', 'K1', 'K3']
SpT_bin_3 = ['K4', 'K6', 'M1', 'M2', 'M3']
bin_1_inds_US = np.array([])
bin_2_inds_US = np.array([])
bin_3_inds_US = np.array([])
bin_1_inds_UCL = np.array([])
bin_2_inds_UCL = np.array([])
bin_3_inds_UCL = np.array([])
for SpT in SpT_bin_1:
    inds = np.where(np.array(Obj_bayes_SpT[0])==SpT)[0]
    if len(inds) > 0:
        bin_1_inds_US = np.concatenate((bin_1_inds_US, inds))

for SpT in SpT_bin_2:
    inds = np.where(np.array(Obj_bayes_SpT[0])==SpT)[0]
    if len(inds) > 0:
        bin_2_inds_US = np.concatenate((bin_2_inds_US, inds))

for SpT in SpT_bin_3:
    inds = np.where(np.array(Obj_bayes_SpT[0])==SpT)[0]
    if len(inds) > 0:
        bin_3_inds_US = np.concatenate((bin_3_inds_US, inds))

for SpT in SpT_bin_1:
    inds = np.where(np.array(Obj_bayes_SpT[1])==SpT)[0]
    if len(inds) > 0:
        bin_1_inds_UCL = np.concatenate((bin_1_inds_UCL, inds))

for SpT in SpT_bin_2:
    inds = np.where(np.array(Obj_bayes_SpT[1])==SpT)[0]
    if len(inds) > 0:
        bin_2_inds_UCL = np.concatenate((bin_2_inds_UCL, inds))

for SpT in SpT_bin_3:
    inds = np.where(np.array(Obj_bayes_SpT[1])==SpT)[0]
    if len(inds) > 0:
        bin_3_inds_UCL = np.concatenate((bin_3_inds_UCL, inds))

bin_1_frac = np.array([0.46, 0.50, 0.54])
bin_2_frac = np.array([0.38, 0.41, 0.44])
bin_3_frac = np.array([0.34, 0.38, 0.42])
(len(bin_1_inds_US)*0.31*bin_1_frac + len(bin_2_inds_US)*0.31*bin_2_frac + len(bin_3_inds_US)*0.31*bin_3_frac)/len(all_bayes[0])
(len(bin_1_inds_UCL)*0.31*bin_1_frac + len(bin_2_inds_UCL)*0.31*bin_2_frac + len(bin_3_inds_UCL)*0.31*bin_3_frac)/len(all_bayes[1])
0.13229459 - 0.12043919, 0.12043919 - 0.10858378
0.13259298 - 0.12138947, 0.12138947 - 0.11018596
for x_gamma in range(len(gamma)-2)[1:]:
    area = np.trapz(p_gamma_US_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
    delta = np.abs(area - 0.16)
    if delta > prev_delta:
        break
prev_delta = 1.
for x_gamma in range(len(gamma)-2)[1:]:
    area = np.trapz(p_gamma_US_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
    delta = np.abs(area - 0.16)
    if delta > prev_delta:
        percentile = gamma[x_gamma-1]
        break
delta
for x_gamma in range(len(gamma)-2)[1:]:
    area = np.trapz(p_gamma_US_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
    delta = np.abs(area - 0.16)
    if delta > prev_delta:
        percentile = gamma[x_gamma-1]
        break
    prev_delta = delta
percentile
gamma[122]
prev_delta = 1.
for x_gamma in range(len(gamma)-2)[1:]:
    area = np.trapz(p_gamma_UCL_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
    delta = np.abs(area - 0.16)
    if delta > prev_delta:
        percentile = gamma[x_gamma-1]
        break
    prev_delta = delta
percentile
prev_delta = 1.
for x_gamma in range(len(gamma)-2)[1:]:
    area = np.trapz(p_gamma_US_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
    delta = np.abs(area - 0.84)
    if delta > prev_delta:
        percentile = gamma[x_gamma-1]
        break
    prev_delta = delta
percentile
prev_delta = 1.
for x_gamma in range(len(gamma)-2)[1:]:
    area = np.trapz(p_gamma_UCL_norm[:x_gamma], x=gamma[:x_gamma], dx=(gamma[1]-gamma[0]))
    delta = np.abs(area - 0.84)
    if delta > prev_delta:
        percentile = gamma[x_gamma-1]
        break
    prev_delta = delta
percentile
12.2 + (27-12.2)/2.
10.2 + (26.9-10.2)/2.
mean_US
mean_US - 12.2
mean_US*100. - 12.2
27 - mean_US*100.
mean_UCL*100. - 10.2
26.9 - mean_UCL*100.
prev_delta = 1.
