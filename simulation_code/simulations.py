#!/usr/bin/python3
#coding=utf-8
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import itertools
import argparse

import pandas as pd
import numpy as np
import pickle
import cProfile

import numba
import multiprocessing

from joblib import Parallel, delayed

from scipy.stats import beta
from scipy import spatial 

import math
import datetime


###SETUP FUNCTIONS

#@numba.jit(nopython=True)
def iota(n,N):
    """ ι : N → R 2 associates with each index n a point in
            the unit circle, evenly spaced, with ι(n) = (cos(N/n) π, sin(N/n) * π)
    """
    pi = math.pi
    return ( math.cos(n/N) * pi, math.sin(n/N) * pi )

@numba.jit(nopython=True)
def hop_distance(i,j):
    """ the minimum # of hops from node i to node j
    """
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))

def cov_mat_fun(sigm, rho, N):
    """ Create covariance matrix
    Returns:
        numpy.ndarray
    """
    cov_mat = np.ones((N,N))
    for i,j in np.ndindex(cov_mat.shape):
        iota_i = iota(i,N)
        iota_j = iota(j,N)
        dist = spatial.distance.euclidean(iota_i, iota_j)
        cov_mat[i,j] = rho**dist
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

@numba.jit(nopython=True)
def init_sigma(x1,x2,Sigma_Ui, Cit, Nit):
    Sigma11 = np.ones((len(x1),len(x1)), dtype=np.float64)
    Sigma12 = np.ones((len(x1),len(x2)), dtype=np.float64)
    Sigma21 = np.ones((len(x2),len(x1)), dtype=np.float64)
    Sigma22 = np.ones((len(x2),len(x2)), dtype=np.float64)
    
    for i in range(len(Cit)):
        for j in range(len(Cit)):
            Sigma11[i,j] = Sigma_Ui[Cit[i],Cit[j]] 
        for j in range(len(Nit)):
            Sigma21[j,i] = Sigma_Ui[Cit[i], Nit[j]] 

    for i in range(len(Nit)):
        for j in range(len(Cit)):
            Sigma12[j,i] = Sigma_Ui[Nit[i],Cit[j]] 
        for j in range(len(Nit)):
            Sigma22[i,j] = Sigma_Ui[Nit[i], Nit[j]]

    return Sigma11, Sigma12, Sigma21, Sigma22

def get_mubar_sigmamu(Sigma_Ui, Ui, x1, Sigma11, Sigma12, Sigma21, Sigma22, mu1, mu2):
    a = np.array([Ui[n] for n in x1]).reshape((1,len(x1)))
    inv_mat = inv_nla_jit(Sigma11)
    inner = np.matmul(Sigma21, inv_mat)
    mubar = mu2 + (np.matmul(inner,(a-mu1).T)).T
    sigmabar = Sigma22 - np.matmul(inner, Sigma12)
    mu_new = Ui
    sigma_new = Sigma_Ui
    return mu_new, Sigma_Ui, sigmabar, mubar

    
@numba.jit(nopython=True)
def get_sigma_new_mu_new(x2, sigmabar, mu_new, sigma_new, mubar):
    len_x2 = len(x2)
    for i in range(len_x2):
        mu_new[x2[i]] = mubar[0,i]
        for j in range(len_x2):
            sigma_new[x2[i], x2[j]] = sigmabar[i, j]
    return mu_new, sigma_new

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

    mu1 = np.array([ce_Ui[n] for n in x1], dtype=np.float64).reshape((1,len(x1)))
    mu2 = np.array([ce_Ui[n] for n in x2], dtype=np.float64).reshape((1,len(x2)))
    
    Sigma11, Sigma12, Sigma21, Sigma22 = init_sigma(x1,x2, Sigma_Ui, Cit, Nit)
    
    #a = np.array([Ui[n] for n in x1]).reshape((1,len(x1)))
    #inv_mat = inv_nla_jit(Sigma11)
    #inner = np.matmul(Sigma21, inv_mat)
    #mubar = mu2 + (np.matmul(inner,(a-mu1).T)).T
    #sigmabar = Sigma22 - np.matmul(inner, Sigma12)
    #mu_new = Ui
    #sigma_new = Sigma_Ui
    mu_new, sigma_new, sigmabar, mubar = get_mubar_sigmamu(Sigma_Ui, Ui, x1, Sigma11, Sigma12, Sigma21, Sigma22, mu1, mu2)

    #for i in range(len(x2)):
    #    mu_new[x2[i]] = mubar[0,i]
    #    for j in range(len(x2)):
    #        sigma_new[x2[i], x2[j]] = sigmabar[i, j]

    mu_new, sigma_new =  get_sigma_new_mu_new(x2, sigmabar, mu_new, sigma_new, mubar)
    return mu_new, sigma_new

## CHOICE FUNCTIONS
def choice_helper(Cit,mu, choice_set):
    """
    Args:
        Cit (list):
        mu (numpy.ndarray):
        choice_set (list) :
    """
    cit = choice_set[np.argmax([mu[i] for i in choice_set])]
    return cit

def choice_ind(U_i, mu_U_i, Sigma_U_i, T, N, Nset, alpha, epsilon):
    C_iT = []
    for t in range(T):
        mu_Uit = mu_U_i
        Sigma_Uit = Sigma_U_i
        if len(C_iT) > 0:
            # update beliefs
            mu_Uit, Sigma_Uit = update_Ui(C_iT,U_i,mu_U_i, np.copy(Sigma_U_i), Nset)
        # make choice
        ce_Uit = certainty_equivalent(alpha, mu_Uit, Sigma_Uit)
        choice_set = [n for n in Nset if n not in C_iT]
        c_it = None
        if np.random.uniform() < epsilon:
            c_it = np.random.choice(choice_set)
        else:
            c_it = choice_helper(C_iT,ce_Uit, choice_set)
        C_iT = C_iT + [c_it]
    return C_iT


def choice_omni(U_i,T,N, Nset):
    C_iT = []
    for t in range(T):
        choice_set = [n for n in Nset if n not in C_iT]
        c_it = choice_helper(C_iT,U_i, choice_set)
        C_iT = C_iT + [c_it]
    return C_iT


def choice_part(V_i, mu_V_i,Sigma_V_i,V,T,N, Nset, alpha, epsilon):
    C_iT = []
    R_iT = []
    for t in range(T):
        mu_Vit = mu_V_i
        Sigma_Vit = Sigma_V_i
        if len(C_iT) > 0:
            # update beliefs
            mu_Vit, Sigma_Vit = update_Ui(C_iT,V_i,mu_V_i, np.copy(Sigma_V_i), Nset)
        mu_Uit = mu_Vit + beta * V
        # make choice
        ce_Uit = certainty_equivalent(alpha, mu_Uit, Sigma_Vit) # γ: uncertainty aversion
        choice_set = [n for n in Nset if n not in C_iT]
        c_it = None
        if np.random.uniform() < epsilon:
            c_it = np.random.choice(choice_set)
        else:
            c_it = choice_helper(C_iT,ce_Uit, choice_set)
        r_it = choice_set[np.argmax([V[i] for i in choice_set])]
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
        mu_V_i = np.copy(mu_V_ibar)
        V_i = np.random.multivariate_normal(mu_V_i, Sigma_V_i)
        mu_V_i.reshape((1,N))

        # Utility in vector form
        U_i = V_i + (beta * V)
        mu_U_i = mu_V_i + beta * mu_V
        
        ## NO RECOMMENDATION CASE
        if beta != 0:
            Sigma_U_i = Sigma_V_i + beta**2 * (Sigma_V)
        else:
            Sigma_U_i = Sigma_V_i

        C_iT = choice_ind(U_i,np.copy(mu_U_i), Sigma_U_i,T,N, Nset, alpha, epsilon)
        C_pop[NO_REC] += [C_iT]
        w_val = w_fun(C_iT,U_i)
        W_pop[NO_REC] += [w_val]

        ## OMNISCIENT CASE
        C_iT = choice_omni(U_i,T,N, Nset)
        C_pop[OMNI] += [C_iT]
        w_val = w_fun(C_iT,U_i)
        W_pop[OMNI] += [w_val]

        ## PARTIAL REC Case
        mu_V_i = np.copy(mu_V_ibar)
        mu_V_i.reshape((1,N))
        C_iT, R_iT = choice_part(V_i,mu_V_i, Sigma_V_i,V,T,N, Nset, alpha, epsilon)
        C_pop[PARTIAL] += [C_iT]
        w_val = w_fun(C_iT,U_i)
        W_pop[PARTIAL] += [w_val]
        R_pop[PARTIAL] += [R_iT]

    return { 'Consumption': C_pop, 'Welfare': W_pop, 'Rec': R_pop }




OMNI = 'rec'
PARTIAL = 'partial'
NO_REC = 'no_rec'

nr_pop = 50
nr_ind = 100
sigma_ibar = .1
rho_ibar = 0

num_cores = multiprocessing.cpu_count() - 1
sim_results ={}

N_vals = [200]

T_vals = [20]

# Covariance structure
rho_vals = [0.1, 0.3, 0.5, 0.7, 0.9]

# utility idiosyncratic degree 
beta_vals = [0, 1, 2, 10]

# absolute risk aversion
alpha_vals = [0, 1]

# action of the time for random exploration
epsilon_vals = [0, 0.1]

sigma_vals = [0.25, 1]

if __name__ == '__main__':
    # itertools.product lists to get all permutations
    vals = [N_vals, T_vals, rho_vals, beta_vals, sigma_vals, alpha_vals, epsilon_vals]
    params = list(itertools.product(*vals))

    print(len(params))

    parser = argparse.ArgumentParser()
    parser.add_argument('--no-joblib', dest='no_joblib', action='store_true',  help=' turn off joblib ')
    args = parser.parse_args()
    print(args)

    for N, T, rho, beta, sigma, alpha, epsilon in params:
        print("STARTING")
        print("N: {}, T: {}, ρ: {} β: {} σ: {} α:{}  ε:{}".format(N, T, rho, beta, sigma, alpha, epsilon))
        print("@ {}".format(datetime.datetime.now()))
        sigma_i = sigma


        Sigma_V_i = cov_mat_fun(sigma_i,rho,N)
        Sigma_V = cov_mat_fun(sigma,rho,N)
        if beta != 0:
            Sigma_V = Sigma_V / beta**2
        Sigma_V_ibar = cov_mat_fun(sigma_ibar,rho_ibar,N)


        if args.no_joblib:
            sim_results[(N, T, rho, beta, sigma, alpha, epsilon)] = [ simulate( N,
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
                                                                        epsilon, seed=i+1) for i in range(1) ]
        else:
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
        with open('data/sim_results.p', 'wb') as fp:
            pickle.dump(sim_results, fp)
    with open('data/sim_results.p', 'wb') as fp:
        pickle.dump(sim_results, fp)
