#!/usr/bin/python3
#coding=utf-8
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import itertools

import pandas as pd
import numpy as np
import pickle
import cProfile
from copy import copy

import numba

from joblib import Parallel, delayed

from scipy.stats import beta

mycolor1 = [0.368417, 0.506779, 0.709798]
mycolor2 = [0.880722, 0.611041, 0.142051]
mycolor3 = [0.560181, 0.691569, 0.194885]
mycolor5 = [0.528488, 0.470624, 0.701351]
mycolor6 = [0.772079, 0.431554, 0.102387]
mycolor7 = [0.363898, 0.618501, 0.782349]
mycolor8 = [1, 0.75, 0]
mycolor9 = [0.647624, 0.37816, 0.614037]
mycolor10 = [0.571589, 0.586483, 0.]


###SETUP FUNCTIONS
def d(i,j):
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))

def cov_mat_fun(sigm, rho, N):
    """ Create covariance matrix
    Returns:
        numpy.ndarray
    """
    cov_mat = np.ones((N,N))
    for i in range(0,N):
        for j in range(0,N):
            cov_mat[i,j] = rho**d(i,j)
    cov_mat = sigm * cov_mat
    return cov_mat

def certainty_equivalent(alpha, mu, sigma):
    """ CARA
    Args:
        alpha (int) :
        mu (numpy.ndarray): shape (N,)
        sigma (numpy.ndarray) : shape (N,N)
    Returns:
        (numpy.ndarray) of shape (N,)
    Notes:
        alpha: the coefficient of absolute risk aversion
        μ and σ2 are the mean and the variance of the distribution F
        https://ocw.mit.edu/courses/economics/14-123-microeconomic-theory-iii-spring-2015/lecture-notes-and-slides/MIT14_123S15_Chap3.pdf 
        pg 21
    """
    new_mu = mu - (.5 * alpha * sigma.diagonal()**2) 
    return new_mu

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

def update_Ui(Cit, Ui, ce_Ui, Sigma_Ui, Nset):
    """ Update beliefs using Baysian updating on Normal Distribution
    Args:
        Cit () :
        Ui () :
        ce_Ui (numpy.ndarray) : certainty equivelents. previously mu_Ui. Shape (N,)
        Sigma_Ui () :
        Nset (range) :
    """ 
    
    # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
    # μ_bar = μ_1 + Ε12 Ε22^-1 ( a - μ_2 )  

    x1 = Cit
    x2 = [n for n in Nset if n not in Cit]
    Nit = [n for n in Nset if n not in Cit]

    mu1 = np.array([ce_Ui[n] for n in x1]).reshape((1,len(x1)))
    mu2 = np.array([ce_Ui[n] for n in x2]).reshape((1,len(x2)))

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
    ce_new = Ui

    for i in range(len(x2)):
        ce_new[x2[i]] = mubar[0,i]
    return ce_new

## CHOICE FUNCTIONS
def choice_helper(Cit,mu, Nset):
    """
    Args:
        Cit (list):
        mu (numpy.ndarray):
        Nset (range) :
    """
    x2 = [n for n in Nset if n not in Cit]
    cit = x2[np.argmax([mu[i] for i in x2])]
    return cit

def choice_ind(U_i, mu_U_i, Sigma_U_i, T, N, Nset, alpha):
    C_iT = []
    for t in range(T):
        mu_Uit = mu_U_i
        if len(C_iT) > 0:
            # update beliefs
            mu_Uit = update_Ui(C_iT,U_i,mu_U_i, Sigma_U_i, Nset)
        # make choice
        gamma_mu_Uit = certainty_equivalent(alpha, mu_Uit, Sigma_U_i) # γ: uncertainty aversion
        c_it = choice_helper(C_iT,gamma_mu_Uit, Nset)
        C_iT = C_iT + [c_it]
    return C_iT


def choice_omni(U_i,T,N, Nset):
    C_iT = []
    for t in range(T):
        c_it = choice_helper(C_iT,U_i, Nset)
        C_iT = C_iT + [c_it]
    return C_iT


def choice_part(V_i, mu_V_i,Sigma_V_i,V,T,N, Nset, alpha):
    C_iT = []
    R_iT = []
    for t in range(T):
        mu_Vit = mu_V_i
        if len(C_iT) > 0:
            # update beliefs
            mu_Vit = update_Ui(C_iT,V_i,mu_V_i, Sigma_V_i, Nset)
        mu_Uit = beta*mu_Vit+(1-beta)*V
        # make choice
        gamma_mu_Uit = certainty_equivalent(alpha, mu_Uit, Sigma_V_i) # γ: uncertainty aversion
        c_it = choice_helper(C_iT,gamma_mu_Uit, Nset)
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
    alpha,
    epsilon,
    seed=1.0
):
    print("iteration: %s " % seed)

    np.random.seed(int(seed))

    Nset = range(0,N)   # set of N items I = {0, 1, ..., N − 1}


    # V = (v_n) n in I aka: common value component v_n in vector form
    mu_V = np.zeros(N)
    V = np.random.multivariate_normal(mu_V, Sigma_V)
    mu_V.reshape((1,N))

    C_pop = { NO_REC: [], OMNI: [], PARTIAL: []}
    W_pop = { NO_REC: [], OMNI: [], PARTIAL: []}
    R_pop = { NO_REC: [], OMNI: [], PARTIAL: []}

    for it_ind in range(nr_ind):

        # V_i = (v_in) n in I aka: consumer i’s idiosyncratic taste for good n in vector form
        mu_V_ibar = np.random.multivariate_normal(np.zeros(N), Sigma_V_ibar)
        mu_V_i = copy(mu_V_ibar)
        V_i = np.random.multivariate_normal(mu_V_i, Sigma_V_i)
        mu_V_i.reshape((1,N))

        # Utility in vector form
        U_i = V_i + (beta * V)
        mu_U_i = beta * mu_V_i + (1-beta) * mu_V
        Sigma_V = Sigma_V * (1/beta) # scale sigma V by 1/beta

        ## NO RECOMMENDATION CASE
        Sigma_U_i = beta**2 * Sigma_V_i + (1-beta)**2 * Sigma_V
        C_iT = choice_ind(U_i,copy(mu_U_i), Sigma_U_i,T,N, Nset, alpha)
        C_pop[NO_REC] += [C_iT]
        w_val = w_fun(C_iT,U_i)
        W_pop[NO_REC] += [w_val]
        print(C_iT)

        ## OMNISCIENT CASE
        C_iT = choice_omni(U_i,T,N, Nset)
        C_pop[OMNI] += [C_iT]
        w_val = w_fun(C_iT,U_i)
        W_pop[OMNI] += [w_val]
        print(C_iT)

        ## PARTIAL REC Case
        mu_V_i = copy(mu_V_ibar)
        mu_V_i.reshape((1,N))
        C_iT, R_iT = choice_part(V_i,mu_V_i, Sigma_V_i,V,T,N, Nset, alpha)
        C_pop[PARTIAL] += [C_iT]
        w_val = w_fun(C_iT,U_i)
        W_pop[PARTIAL] += [w_val]
        R_pop[PARTIAL] += [R_iT]

    return { 'Consumption': C_pop, 'Welfare': W_pop, 'Rec': R_pop }




OMNI = 'rec'
PARTIAL = 'partial'
NO_REC = 'no_rec'

N = 500      # Number of items to choose from
T = 10      # total items a consumer consumes / num time periods
nr_pop = 1
nr_ind = 5
sigma_ibar = .1
rho_ibar = 0

num_cores = 1
sim_results ={}

# Covariance structure
rho_vals = [0.1, 0.3, 0.5, 0.7, 0.9]

# utility idiosyncratic degree 
beta_vals = [0.1]

# absolute risk aversion
alpha_vals = [0, 1]

# action of the time for random exploration
epsilon_vals = [0, 1/10, 1/4]

sigma_vals = [0.25]

# itertools.product lists to get all permutations
vals = [rho_vals, beta_vals, sigma_vals, alpha_vals, epsilon_vals]
params = list(itertools.product(*vals))

for rho, beta, sigma, alpha, epsilon in params:
    print("STARTING")
    print("ρ: {} β: {} σ: {} α:{}  ε:{}".format(rho, beta, sigma, alpha, epsilon))
    sigma_i = sigma

    Sigma_V_i = cov_mat_fun(sigma_i,rho,N)
    Sigma_V = cov_mat_fun(sigma,rho,N)
    Sigma_V_ibar = cov_mat_fun(sigma_ibar,rho_ibar,N)

    sim_results[(N, T, rho, beta, sigma, alpha, epsilon)] = Parallel(n_jobs=num_cores)(delayed(simulate)(N,
                                                                T,
                                                                sigma,
                                                                sigma_i,
                                                                sigma_ibar,
                                                                beta,
                                                                nr_ind,
                                                                Sigma_V_i, 
                                                                Sigma_V, 
                                                                Sigma_V_ibar, 
                                                                alpha, 
                                                                epsilon, seed=i+1) for i in range(nr_pop))
    print("finished a population run")
    with open('sim_results_test.p', 'wb') as fp:
        pickle.dump(sim_results, fp)

with open('sim_results_test.p', 'wb') as fp:
    pickle.dump(sim_results, fp)
