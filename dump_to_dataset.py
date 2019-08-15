import pickle
import csv
import json
from copy import copy
import numpy as np


from scipy.spatial.distance import jaccard, euclidean 
from simulations import iota

OMNI = 'omni'
PARTIAL = 'partial'
NO_REC = 'no_rec'
rec_policy_keys = [OMNI, PARTIAL, NO_REC]
#WORKING_DIR = '/Users/guyaridor/Desktop/ExAnteFilterBubble/data/'
WORKING_DIR = '/home/guyaridor/ExAnteFilterBubble/data/'
with open(WORKING_DIR + 'new_sim.json', 'r') as fp:
    df = json.load(fp)

def d(i,j, N):
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))

def div_fun(CiT, T, N):
    div_score = 0.0
    for i in range(len(CiT)):
        for j in range(len(CiT)):
            if i == j: continue
            div_score = div_score + d(CiT[i],CiT[j], N)#euclidean( iota(CiT[i], N), iota(CiT[j], N) )
    return div_score*(1.0/(T*(T-1)))

def homogeneity(C1, C2, type_sim="jaccard"):
    if type_sim == "jaccard":
        return 1 - jaccard(C1, C2)
    return None

def follow_rec(Ci, Rec, N, T):
    s = 0.0
    for i in range(len(Ci)):
        if Ci[i] == Rec[i]:
            s += 1
    return s / T

def parse_pickle_key(key):
    key = eval(key)
    dat = {
        'N': float(key[0]),
        'T': float(key[1]),
        'rho': key[2],
        'beta': key[3],
        'sigma': key[4],
        'alpha': key[5],
        'epsilon': key[6]
        #'pop_idx': pop_idx
    }
    return dat

INDIVIDUAL_FIELD_NAMES =['pop_idx', 'indiv_idx', 'regime', 'welfare', 'diversity_score', 'rho', 'beta', 'epsilon', 'follow_recommendation', 'N', 'T', 'sigma', 'alpha', 'nr_pop', 'nr_ind']
with open(WORKING_DIR + 'rec_data.csv', 'w') as rec_csv:
    data_writer = csv.DictWriter(rec_csv, fieldnames=INDIVIDUAL_FIELD_NAMES)
    data_writer.writeheader()
    for key, value in df.items():
        for pop_idx in range(len(value)):
            dat = parse_pickle_key(key)
            T = dat['T']
            N = dat['N']
            dat['pop_idx'] = pop_idx
            cur = value[pop_idx]
            welfare = cur['Welfare']
            consumption = cur['Consumption']
            for policy in rec_policy_keys:
                dat['regime'] = policy
                for indiv_idx in range(len(welfare[policy])):
                    dat['indiv_idx'] = indiv_idx
                    dat['welfare'] = welfare[policy][indiv_idx]
                    consumption_arr = np.array(consumption[policy])
                    dat['diversity_score'] = div_fun(consumption_arr[:,indiv_idx], T, N)
                    dat['follow_recommendation'] = False
                    if policy == PARTIAL:
                        follow_rec_arr = np.array(cur['Rec'][policy])
                        dat['follow_recommendation'] = follow_rec(consumption_arr[:,indiv_idx], follow_rec_arr[:,indiv_idx], T, N)
                    cur_dat = copy(dat)
                    data_writer.writerow(cur_dat)

INDIVIDUAL_FIELD_NAMES =['pop_idx', 'regime', 'rho', 'beta', 'N', 'T', 'sigma', 'jaccard', 'beta', 'alpha', 'epsilon', 'nr_pop', 'nr_ind']
with open(WORKING_DIR + 'homogeneity_data.csv', 'w') as rec_csv:
    data_writer = csv.DictWriter(rec_csv, fieldnames=INDIVIDUAL_FIELD_NAMES)
    data_writer.writeheader()
    for key, value in df.items():
        for pop_idx in range(len(value)):
            dat = parse_pickle_key(key)
            T = dat['T']
            N = dat['N']
            cur = value[pop_idx]
            consumption = cur['Consumption']
            for policy in rec_policy_keys:
                dat['regime'] = policy
                consumption_arr = np.array(consumption[policy])
                iter_size = len(consumption_arr[0,:])
                tot = 0.0
                for indiv_idx1 in range(iter_size):
                    for indiv_idx2 in range(iter_size):
                        if indiv_idx1 == indiv_idx2: continue
                        c1 = consumption_arr[:,indiv_idx1]
                        c2 = consumption_arr[:,indiv_idx2]
                        dist = homogeneity(c1, c2, "jaccard")
                        tot += dist
                tot /= (iter_size * (iter_size - 1))
                dat['jaccard'] = tot
                cur_dat = copy(dat)
                data_writer.writerow(cur_dat)


INDIVIDUAL_FIELD_NAMES =['pop_idx', 'regime', 'rho', 'beta', 'epsilon', 'alpha', 'N', 'T', 'sigma', 't', 'follow_recommendation', 'consumption_dist']
with open(WORKING_DIR + 'time_path.csv', 'w') as rec_csv:
    data_writer = csv.DictWriter(rec_csv, fieldnames=INDIVIDUAL_FIELD_NAMES)
    data_writer.writeheader()
    for key, value in df.items():
        for pop_idx in range(len(value)):
            dat = parse_pickle_key(key)
            N = dat['N']
            dat['pop_idx'] = pop_idx
            cur = value[pop_idx]
            consumption = cur['Consumption']
            for policy in rec_policy_keys:
                consumption_arr = np.array(consumption[policy])
                dat['regime'] = policy
                if policy == PARTIAL:
                    follow_rec_arr = np.array(cur['Rec'][policy])
                for indiv_idx in range(len(consumption_arr[0,:])):
                    c1 = consumption_arr[:,indiv_idx]
                    prev_consumption = None
                    rec = None
                    if policy == PARTIAL:
                        rec = follow_rec_arr[:,indiv_idx]
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
