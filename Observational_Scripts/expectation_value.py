%run /home/100/rlk100/Scripts/bayesian_analysis.py -file Table_2_reviewer_response.csv -suffix '_reviewer_response' -bayes_f bayes_factors_reviewer_response.csv
bin_1_inds_US = []
bin_2_inds_US = []
bin_3_inds_US = []
bin_1_inds_UCL = []
bin_2_inds_UCL = []
bin_3_inds_UCL = []
for obj_ind in range(len(Obj_bayes_SpT[0])):
    if RV_standard_info[Obj_bayes_SpT[0][obj_ind]][1] > 1.02:
        bin_1_inds_US.append(obj_ind)
    elif RV_standard_info[Obj_bayes_SpT[0][obj_ind]][1] < 0.78:
        bin_3_inds_US.append(obj_ind)
    else:
        bin_2_inds_US.append(obj_ind)
for obj_ind in range(len(Obj_bayes_SpT[1])):
    if RV_standard_info[Obj_bayes_SpT[1][obj_ind]][1] > 1.02:
        bin_1_inds_UCL.append(obj_ind)
    elif RV_standard_info[Obj_bayes_SpT[1][obj_ind]][1] < 0.78:
        bin_3_inds_UCL.append(obj_ind)
    else:
        bin_2_inds_UCL.append(obj_ind)

bin_1_inds_US = np.array(bin_1_inds_US)
bin_2_inds_US = np.array(bin_2_inds_US)
bin_3_inds_US = np.array(bin_3_inds_US)
bin_1_inds_UCL = np.array(bin_1_inds_UCL)
bin_2_inds_UCL = np.array(bin_2_inds_UCL)
bin_3_inds_UCL = np.array(bin_3_inds_UCL)

bin_1_frac = np.array([0.50, 0.54, 0.58])
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
