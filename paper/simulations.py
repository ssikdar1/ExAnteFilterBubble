import pandas as pd
import numpy as np
import pickle
import cProfile
from copy import copy

from joblib import Parallel, delayed

import itertools
from scipy import stats
from scipy.stats import beta


import pandas as pd
import numpy as np
import pickle
import cProfile
import numba

from joblib import Parallel, delayed

import itertools
from scipy import stats
from scipy.stats import beta

WRITE_FILE = 'sim_results_t_25.p'
OMNI = 'rec'
PARTIAL = 'partial'
NO_REC = 'no_rec'


###SETUP FUNCTIONS
def d(i,j, N):
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))

def cov_mat_fun(sigm, rh, N):
    cov_mat = np.ones((N,N))
    for i in range(0,N):
        for j in range(0,N):
            cov_mat[i,j] = rh**d(i,j, N)
    cov_mat = sigm * cov_mat
    return cov_mat



### Welfare Functions - Statistic Calculation Functions


def w_fun(CiT,Ui):
    w_score = 0.0
    for i in range(len(CiT)):
        w_score = w_score + Ui[CiT[i]]
    return w_score*(T**(-1))


## Bayesian updating
@numba.jit
def inv_nla_jit(A):
  return np.linalg.inv(A)

def update_Ui(Cit,Ui,mu_Ui, Sigma_Ui, Nset):
    x1 = Cit
    x2 = [n for n in Nset if n not in Cit]
    Nit = [n for n in Nset if n not in Cit]
    mu1 = np.array([mu_Ui[n] for n in x1]).reshape((1,len(x1)))
    mu2 = np.array([mu_Ui[n] for n in x2]).reshape((1,len(x2)))
    Sigma11 = np.ones((len(x1),len(x1)))
    Sigma21 = np.ones((len(x2),len(x1)))
    for i in range(len(Cit)):
        for j in range(len(Cit)):
            Sigma11[i,j] = Sigma_Ui[Cit[i],Cit[j]]
    for i in range(len(Cit)):
        for j in range(len(Nit)):
            Sigma21[j,i] = Sigma_Ui[Cit[i], Nit[j]]
    a = np.array([Ui[n] for n in x1]).reshape((1,len(x1)))
    inv_mat = inv_nla_jit(Sigma11)
    inner = np.matmul(Sigma21, inv_mat)
    mubar = mu2 + (np.matmul(inner,(a-mu1).T)).T
    mu_new = mu_Ui
    for i in range(len(x1)):
        mu_new[x1[i]] = Ui[i]
    for i in range(len(x2)):
        mu_new[x2[i]] = mubar[0,i]
    return mu_new

## CHOICE FUNCTIONS
def choice_helper(Cit,mu, Nset):
    x2 = [n for n in Nset if n not in Cit]
    cit = x2[np.argmax([mu[i] for i in x2])]
    return cit

def choice_ind(U_i,mu_U_i, Sigma_U_i,T,N, Nset):
    C_iT = []
    for t in range(T):
        mu_Uit = copy(mu_U_i)
        if len(C_iT) > 0:
            mu_Uit = update_Ui(C_iT,U_i,mu_U_i, Sigma_U_i, Nset)
        c_it = choice_helper(C_iT,mu_Uit, Nset)
        C_iT = C_iT + [c_it]
    return C_iT


def choice_omni(U_i,T,N, Nset):
    C_iT = []
    for t in range(T):
        c_it = choice_helper(C_iT,U_i, Nset)
        C_iT = C_iT + [c_it]
    return C_iT


def choice_part(V_i, mu_V_i,Sigma_V_i,V,T,N, Nset):
    C_iT = []
    R_iT = []
    for t in range(T):
        mu_Vit = mu_V_i
        if len(C_iT) > 0:
            mu_Vit = update_Ui(C_iT,V_i,mu_V_i, Sigma_V_i, Nset)
        mu_Uit = beta*mu_Vit+(1-beta)*V
        c_it = choice_helper(C_iT,mu_Uit, Nset)
        Nit = [n for n in Nset if n not in C_iT]
        r_it = Nit[np.argmax([V[i] for i in Nit])]
        R_iT = R_iT + [r_it]
        C_iT = C_iT + [c_it]

    return C_iT, R_iT

def simulate(
    N,
    T,
    sigma,
    sigma_i,
    sigma_ibar,
    beta,
    nr_ind,
    Sigma_V_i,
    Sigma_V,
    Sigma_V_ibar,
    seed=1.0
):
    print("iteration")
    print(seed)
    np.random.seed(int(seed))
    Nset = range(0,N)
    mu_V = np.zeros(N)
    V = np.random.multivariate_normal(mu_V, Sigma_V)
    mu_V.reshape((1,N))
    C_pop = { NO_REC: [], OMNI: [], PARTIAL: []}
    W_pop = { NO_REC: [], OMNI: [], PARTIAL: []}
    R_pop = { NO_REC: [], OMNI: [], PARTIAL: []}
    for it_ind in range(nr_ind):
        mu_V_ibar = np.random.multivariate_normal(np.zeros(N), Sigma_V_ibar)
        mu_V_i = copy(mu_V_ibar)
        V_i = np.random.multivariate_normal(mu_V_i, Sigma_V_i)
        mu_V_i.reshape((1,N))
        U_i = beta * V_i + (1-beta) * V
        mu_U_i = beta * mu_V_i + (1-beta) * mu_V
        #print(mu_U_i)
        #print(U_i)

        ##No Rec Case
        Sigma_U_i = beta**2 * Sigma_V_i + (1-beta)**2 * Sigma_V
        C_iT = choice_ind(U_i,copy(mu_U_i), Sigma_U_i,T,N, Nset)
        C_pop[NO_REC] += [C_iT]
        w_val = w_fun(C_iT,U_i)
        W_pop[NO_REC] += [w_val]
        #print(C_iT)
        ## OMNI CASE
        C_iT = choice_omni(U_i,T,N, Nset)
        C_pop[OMNI] += [C_iT]
        w_val = w_fun(C_iT,U_i)
        W_pop[OMNI] += [w_val]
        #print(C_iT)
        ## PARTIAL REC Case

        mu_V_i = copy(mu_V_ibar)
        mu_V_i.reshape((1,N))
        C_iT, R_iT = choice_part(V_i,mu_V_i, Sigma_V_i,V,T,N, Nset)
        C_pop[PARTIAL] += [C_iT]
        w_val = w_fun(C_iT,U_i)
        W_pop[PARTIAL] += [w_val]
        R_pop[PARTIAL] += [R_iT]
    return { 'Consumption': C_pop, 'Welfare': W_pop, 'Rec': R_pop }



N = 1000
T = 50
nr_pop = 25
nr_ind = 25
sigma_ibar = 0.1
rho_ibar = 0

num_cores = 8
sim_results ={}
rho_vals = [0.1, 0.5, 0.9]
beta_vals = [0.1, 0.5, 0.9]
sigma_vals = [0.25, 1.0, 4.0]
for rho in rho_vals:
    for beta in beta_vals:
        for sigma in sigma_vals:
            print(rho, beta, sigma)
            sigma_i = sigma
            Sigma_V_i = cov_mat_fun(sigma_i,rho,N)
            Sigma_V = cov_mat_fun(sigma,rho,N)
            Sigma_V_ibar = cov_mat_fun(sigma_ibar,rho_ibar,N)
            sim_results[(N, T, rho, beta, sigma)] = Parallel(n_jobs=num_cores)(delayed(simulate)(N,T,sigma,sigma_i,sigma_ibar,beta,nr_ind,Sigma_V_i, Sigma_V, Sigma_V_ibar, seed=i+1) for i in range(nr_pop))
            print("finished a population run")
            with open(WRITE_FILE, 'wb') as fp:
                pickle.dump(sim_results, fp)

with open(WRITE_FILE, 'wb') as fp:
    pickle.dump(sim_results, fp)
