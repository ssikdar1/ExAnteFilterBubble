import pickle
import csv
import json
from copy import copy


from scipy.spatial.distance import jaccard, euclidean 
from simulations import iota

OMNI = 'omni'
PARTIAL = 'partial'
NO_REC = 'no_rec'
rec_policy_keys = [OMNI, PARTIAL, NO_REC]
WORKING_DIR = '/Users/guyaridor/Desktop/'
with open(WORKING_DIR + 'new_sim.json', 'r') as fp:
    df = json.load(fp)

def d(i,j, N):
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))

def div_fun(CiT, T, N):
    div_score = 0.0
    for i in range(len(CiT)):
        for j in range(len(CiT)):
            if i == j: continue
            div_score = div_score + d(i,j,N)#euclidean( iota(CiT[i], N), iota(CiT[j], N) )
    return div_score*(1.0/(T*(T-1)))

def homogeneity(C1, C2, type_sim="jaccard"):
    if type_sim == "jaccard":
        return jaccard(C1, C2)
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

INDIVIDUAL_FIELD_NAMES =['pop_idx', 'indiv_idx', 'regime', 'welfare', 'diversity_score', 'rho', 'beta', 'follow_recommendation', 'N', 'T', 'sigma', 'beta', 'alpha', 'epsilon', 'nr_pop', 'nr_ind']
with open(WORKING_DIR + 'rec_data.csv', 'w') as rec_csv:
    data_writer = csv.DictWriter(rec_csv, fieldnames=INDIVIDUAL_FIELD_NAMES)
    data_writer.writeheader()
    for key, value in df.items():
        for pop_idx in range(len(value)):
            dat = parse_pickle_key(key)
            T = dat['T']
            N = dat['N']
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

