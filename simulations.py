import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pickle

from joblib import Parallel, delayed

import itertools
from scipy import stats
from scipy.stats import beta

mycolor1 = [0.368417, 0.506779, 0.709798]
mycolor2 = [0.880722, 0.611041, 0.142051]
mycolor3 = [0.560181, 0.691569, 0.194885]
mycolor4 = [0.922526, 0.385626, 0.209179]
mycolor5 = [0.528488, 0.470624, 0.701351]
mycolor6 = [0.772079, 0.431554, 0.102387]
mycolor7 = [0.363898, 0.618501, 0.782349]
mycolor8 = [1, 0.75, 0]
mycolor9 = [0.647624, 0.37816, 0.614037]
mycolor10 = [0.571589, 0.586483, 0.]


OMNI = 'rec'
PARTIAL = 'partial'
NO_REC = 'no_rec'

sigma = 1
sigma_i = 1
sigma_ibar = .1

rho_ibar = 0


N = 1000
T = 25
nr_pop = 1000
nr_ind = 1000

num_cores = 16



###SETUP FUNCTIONS
def d(i,j):
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))


def cov_mat_fun(sigm, rh, N):
    cov_mat = np.ones((N,N))
    for i in range(0,N):
        for j in range(0,N):
            cov_mat[i,j] = rh**d(i,j)
    cov_mat = sigm * cov_mat
    return cov_mat



### Welfare and Diversity Functions - Statistic Calculation Functions
def div_fun(CiT):
    div_score = 0
    for i in range(len(CiT)):
        for j in range(len(CiT)):
            div_score = div_score + d(i,j)
    return div_score*(T**(-2))


def w_fun(CiT,Ui):
    w_score = 0
    for i in range(len(CiT)):
        w_score = w_score + Ui[i]
    return w_score*(T**(-1))


## Bayesian updating

def update_Ui(Cit,Ui,mu_Ui, Sigma_Ui, Nset):
    x1 = Cit
    x2 = [n for n in Nset if n not in Cit]
    Nit = [n for n in Nset if n not in Cit]
    mu1 = np.array([mu_Ui[n] for n in x1]).reshape((1,len(x1)))
    mu2 = np.array([mu_Ui[n] for n in x2]).reshape((1,len(x2)))
    Sigma11 = np.ones((len(x1),len(x1)))
    Sigma22 = np.ones((len(x2),len(x2)))
    Sigma12 = np.ones((len(x1),len(x2)))
    for i in range(len(Cit)):
        for j in range(len(Cit)):
            Sigma11[i,j] = Sigma_Ui[Cit[i],Cit[j]]
    for i in range(len(Nit)):
        for j in range(len(Nit)):
            Sigma22[i,j] = Sigma_Ui[Nit[i],Nit[j]]
    for i in range(len(Cit)):
        for j in range(len(Nit)):
            Sigma12[i,j] = Sigma_Ui[Cit[i],Nit[j]]
    Sigma21 = np.transpose(Sigma12)
    a = np.array([Ui[n] for n in x1]).reshape((1,len(x1)))
    mubar = mu2 + (np.matmul(np.matmul(Sigma21, np.linalg.inv(Sigma11)),(a-mu1).T)).T
    mu_new = mu_Ui
    for i in range(len(x1)):
        mu_new[x1[i]] = mubar[0,i]
    return mu_new



## CHOICE FUNCTIONS 
def choice_helper(Cit,mu, Nset):
    x2 = [n for n in Nset if n not in Cit]
    cit = x2[np.argmax([mu[i] for i in x2])]
    return cit

def choice_ind(U_i,mu_U_i, Sigma_U_i,T,N, Nset):
    C_iT = []
    for t in range(T):
        mu_Uit = mu_U_i
        c_it = choice_helper(C_iT,mu_Uit, Nset)
        C_iT = C_iT + [c_it]
        mu_Uit = update_Ui(C_iT,U_i,mu_U_i, Sigma_U_i, Nset)
    return C_iT


def choice_omni(U_i,T,N, Nset):
    C_iT = []
    for t in range(T):
        c_it = choice_helper(C_iT,U_i, Nset)
        C_iT = C_iT + [c_it]
    return C_iT


def choice_part(U_i,mu_U_i, Sigma_U_i,V,T,N, Nset):
    C_iT = []
    R_iT = []
    for t in range(T):
        mu_Uit = mu_U_i
        c_it = choice_helper(C_iT,mu_Uit, Nset)
        Nit = [n for n in Nset if n not in C_iT]
        r_it = Nit[np.argmax([V[i] for i in Nit])]
        R_iT = R_iT + [r_it]
        C_iT = C_iT + [c_it]
        mu_Uit = update_Ui(C_iT,U_i,mu_U_i, Sigma_U_i, Nset)
        
    return C_iT, R_iT

def simulate(
	N,
	T,
	sigma,
	sigma_i,
	sigma_ibar,
	rho,
	rho_i,
	rho_ibar,
	beta,
	nr_ind,
	seed=1.0
):
	np.random.seed(int(seed))
	Nset = range(0,N)
	Sigma_V_i = cov_mat_fun(sigma_i,rho_i,N)
	Sigma_V = cov_mat_fun(sigma,rho,N)
	Sigma_V_ibar = cov_mat_fun(sigma_ibar,rho_ibar,N)
	mu_V = np.zeros(N)
	V = np.random.multivariate_normal(mu_V, Sigma_V)
	mu_V.reshape((1,N))
	C_pop = { NO_REC: [], OMNI: [], PARTIAL: []}
	W_pop = { NO_REC: [], OMNI: [], PARTIAL: []}
	R_pop = { NO_REC: [], OMNI: [], PARTIAL: []}

	for it_ind in range(nr_ind):
		mu_V_ibar = np.zeros(N)
		V_ibar = np.random.multivariate_normal(np.zeros(N), Sigma_V_ibar)
        
		mu_V_i = V_ibar
   
		V_i = np.random.multivariate_normal(mu_V_i, Sigma_V_i)
		mu_V_ibar.reshape((1,N))
		mu_V_i.reshape((1,N))

		U_i = beta * V_i + (1-beta) * V
		mu_U_i = beta * mu_V_i + (1-beta) * mu_V

		##No Rec Case
		Sigma_U_i = beta**2 * Sigma_V_ibar + (1-beta)**2 * Sigma_V
		C_iT = choice_ind(U_i,mu_U_i, Sigma_U_i,T,N, Nset)
		C_pop[NO_REC] += [C_iT]
		W_pop[NO_REC] += [w_fun(C_iT,U_i)]

		## OMNI CASE
		C_iT = choice_omni(U_i,T,N, Nset)
		C_pop[OMNI] += [C_iT]
		W_pop[OMNI] += [w_fun(C_iT,U_i)]

		## PARTIAL REC Case
		Sigma_U_i = beta**2 * Sigma_V_ibar
		C_iT, R_iT = choice_part(U_i,mu_U_i, Sigma_U_i,V,T,N, Nset)
		C_pop[PARTIAL] += [C_iT]
		W_pop[PARTIAL] += [w_fun(C_iT,U_i)]
		R_pop[PARTIAL] += [R_iT]
	return { 'Consumption': C_pop, 'Welfare': W_pop, 'Rec': R_pop }

sim_results ={}
for rho in [.1, .5, .9]:
	for beta in [.1, .5, .9]:
		sim_results[(rho, beta)] = Parallel(n_jobs=num_cores)(delayed(simulate)(N,T,sigma,sigma_i,sigma_ibar,rho,rho,rho_ibar,beta,nr_ind,seed=i+1) for i in range(nr_pop))


with open('sim_results.p', 'wb') as fp:
	pickle.dump(sim_results, fp)

