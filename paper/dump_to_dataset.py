import pickle
import csv
from copy import copy

from scipy.spatial.distance import jaccard

OMNI = 'rec'
PARTIAL = 'partial'
NO_REC = 'no_rec'
rec_policy_keys = [OMNI, PARTIAL, NO_REC]
WORKING_DIR = ''
with open(WORKING_DIR + 'sim_results_t_25.p', 'rb') as fp:
    df = pickle.load(fp)

def d(i,j, N):
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))

def div_fun(CiT, T, N):
    div_score = 0.0
    for i in range(len(CiT)):
        for j in range(len(CiT)):
            if i == j: continue
            div_score = div_score + d(CiT[i],CiT[j], N)
    return div_score*(1.0/(T*(T-1)))

def homogeneity(C1, C2, type_sim="jaccard"):
    if type_sim == "jaccard":
        return 1 - jaccard(C1, C2)
    return None

def follow_rec(Ci, Rec, T, N):
    s = 0.0
    for i in range(len(Ci)):
        if Ci[i] == Rec[i]:
            s += 1
    return s / len(Ci)


INDIVIDUAL_FIELD_NAMES =['pop_idx', 'indiv_idx', 'regime', 'welfare', 'diversity_score', 'rho', 'beta', 'follow_recommendation', 'N', 'T', 'sigma']
with open(WORKING_DIR + 'rec_data_t_25.csv', 'w') as rec_csv:
    data_writer = csv.DictWriter(rec_csv, fieldnames=INDIVIDUAL_FIELD_NAMES)
    data_writer.writeheader()
    for key, value in df.items():
        for pop_idx in range(len(value)):
            N = key[0]
            T = key[1]
            dat = {
                'pop_idx': pop_idx,
                'N': N,
                'T': T,
                'rho': key[2],
                'beta': key[3],
                'sigma': key[4]
            }
            cur = value[pop_idx]
            welfare = cur['Welfare']
            consumption = cur['Consumption']
            for policy in rec_policy_keys:
            	dat['regime'] = policy
            	for indiv_idx in range(len(welfare[policy])):
                    dat['indiv_idx'] = indiv_idx
                    dat['welfare'] = welfare[policy][indiv_idx]
                    dat['diversity_score'] = div_fun(consumption[policy][indiv_idx], T, N)
                    dat['follow_recommendation'] = False
                    if policy == PARTIAL:
                        dat['follow_recommendation'] = follow_rec(consumption[policy][indiv_idx], cur['Rec'][policy][indiv_idx], T, N)
                    cur_dat = copy(dat)
                    data_writer.writerow(cur_dat)


INDIVIDUAL_FIELD_NAMES =['pop_idx', 'regime', 'rho', 'beta', 'N', 'T', 'sigma', 'jaccard']
with open(WORKING_DIR + 'homogeneity_data_t_25.csv', 'w') as rec_csv:
    data_writer = csv.DictWriter(rec_csv, fieldnames=INDIVIDUAL_FIELD_NAMES)
    data_writer.writeheader()
    for key, value in df.items():
        for pop_idx in range(len(value)):
            N = key[0]
            T = key[1]
            dat = {
                'pop_idx': pop_idx,
                'N': N,
                'T': T,
                'rho': key[2],
                'beta': key[3],
                'sigma': key[4]
            }
            cur = value[pop_idx]
            consumption = cur['Consumption']
            for policy in rec_policy_keys:
                dat['regime'] = policy
                iter_size = len(consumption[policy])
                tot = 0.0
                for indiv_idx1 in range(iter_size):
                    for indiv_idx2 in range(iter_size):
                        if indiv_idx1 == indiv_idx2: continue
                        c1 = consumption[policy][indiv_idx1]
                        c2 = consumption[policy][indiv_idx2]
                        dist = homogeneity(c1, c2, "jaccard")
                        tot += dist
                tot /= (iter_size * (iter_size - 1))
                dat['jaccard'] = tot
                cur_dat = copy(dat)
                data_writer.writerow(cur_dat)

INDIVIDUAL_FIELD_NAMES =['pop_idx', 'regime', 'rho', 'beta', 'N', 'T', 'sigma', 't', 'follow_recommendation', 'consumption_dist']
with open(WORKING_DIR + 'time_path_25.csv', 'w') as rec_csv:
    data_writer = csv.DictWriter(rec_csv, fieldnames=INDIVIDUAL_FIELD_NAMES)
    data_writer.writeheader()
    for key, value in df.items():
        for pop_idx in range(len(value)):
            N = key[0]
            T = key[1]
            dat = {
                'pop_idx': pop_idx,
                'N': N,
                'T': T,
                'rho': key[2],
                'beta': key[3],
                'sigma': key[4]
            }
            cur = value[pop_idx]
            consumption = cur['Consumption']
            for policy in rec_policy_keys:
                dat['regime'] = policy
                for indiv_idx in range(len(consumption)):
                    c1 = consumption[policy][indiv_idx]
                    prev_consumption = None
                    rec = None
                    if policy == PARTIAL:
                        rec = cur['Rec'][policy][indiv_idx]
                    for t in range(len(c1)):
                        dat['t'] = t
                        if prev_consumption != None:
                            dat['consumption_dist'] = d(prev_consumption, c1[t], N)
                        prev_consumption = c1[t]
                        dat['follow_recommendation'] = False
                        if policy == PARTIAL:
                            dat['follow_recommendation'] = c1[t] == rec[t]
                        cur_dat = copy(dat)
                        data_writer.writerow(cur_dat)
