
using JSON
using IterTools;
using Dates;
using Distances;
using Random;
using Distributions;
using LinearAlgebra;

""" the minimum # of hops from node i to node j
"""

function hop_distance(
        i::Int64,
        j::Int64,
        N::Int64
    )::Int64
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))
end

"""
    cov_mat_fun(sigma::Float64, rho::Float64, N::Int64)
Create covariance matrix
# Arguments
- `sigma::Float64`: 
- `rho::Float64`: Covariance coefficient
- `N::Int`: Number of goods
# Returns
- Covariance matrix 
"""

function cov_mat_fun(sigma::Float64, rho::Float64, N::Int64)::Array{Float64,2}
    cov_mat = zeros(Float64, N, N)
    for i in 1:N 
        for j in 1:N
            dist = hop_distance( i, j , N)
            #dist = Distances.euclidean( iota(i,N), iota(j,N) )
            cov_mat[i,j] = sigma * rho^dist
        end
    end
    return cov_mat
end

""" 
CARA
# Arguments:
- alpha (int) :
- mu (numpy.ndarray): shape (N,)
- sigma (numpy.ndarray) : shape (N,N)
# Returns:
    (numpy.ndarray) of shape (N,)
# Notes:
    alpha: the coefficient of absolute risk aversion
    μ and σ2 are the mean and the variance of the distribution F
    https://ocw.mit.edu/courses/economics/14-123-microeconomic-theory-iii-spring-2015/lecture-notes-and-slides/MIT14_123S15_Chap3.pdf 
    pg 21
"""

function certainty_equivalent(
        alpha::Float64, 
        mu::Array{Float64,1},
        sig::Array{Float64,2}
    )
    new_mu = mu - (.5 * alpha * diag(sig)) 
    return new_mu
end

### Welfare Functions - Statistic Calculation Functions


function init_sigma(Cit::Array{Int64,1},
            Nit::Array{Int64,1},
            Sigma_Ui::Array{Float64,2})

    Sigma11::Array{Float64,2} = ones(Float64, length(Cit), length(Cit))
    Sigma12::Array{Float64,2} = ones(Float64, length(Cit), length(Nit))
    Sigma21::Array{Float64,2} = ones(Float64, length(Nit), length(Cit))
    Sigma22::Array{Float64,2} = ones(Float64, length(Nit), length(Nit))

    for i in 1:length(Cit)
        for j in 1:length(Cit)
            Sigma11[i,j] = Sigma_Ui[Cit[i],Cit[j]] 
        end
        for j in 1:length(Nit)
            Sigma21[j,i] = Sigma_Ui[Cit[i], Nit[j]] 
        end
    end

    for i in 1:length(Nit)
        for j in 1:length(Cit)
            Sigma12[j,i] = Sigma_Ui[Nit[i],Cit[j]] 
        end
        for j in 1:length(Nit)
            Sigma22[i,j] = Sigma_Ui[Nit[i], Nit[j]]
        end
    end
    return Sigma11, Sigma12, Sigma21, Sigma22
end

function get_mubar_sigmamu(
        Sigma_Ui::Array{Float64,2}, 
        Ui::Array{Float64,1}, 
        x1::Array{Int64,1}, 
        Sigma11::Array{Float64,2}, 
        Sigma12::Array{Float64,2}, 
        Sigma21::Array{Float64,2}, 
        Sigma22::Array{Float64,2}, 
        mu1::Array{Float64,1}, 
        mu2::Array{Float64,1}
    )
   
    a = [Ui[n] for n in x1]

    inv_mat = inv(Sigma11)

    inner = Sigma21 * inv_mat

    mubar = mu2 + inner * (a-mu1)
    
    sigmabar = Sigma22 - (inner * Sigma12)

    return mubar, sigmabar
end

"""
    update_Ui()

Bayesian Update

# Arguments

"""

function update_Ui(
            Cit::Array{Int64,1}, 
            Ui::Array{Float64,1}, 
            mu_Ui::Array{Float64,1}, 
            Sigma_Ui::Array{Float64,2}, 
            Nset::Array{Int64,1}
    )
    
    # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
    # μ_bar = μ_1 + Ε12 Ε22^-1 ( a - μ_2 )  

    x1 = Cit
    x2 = [n for n in Nset if n ∉ Cit]

    mu1 = [mu_Ui[n] for n in x1]
    mu2 = [mu_Ui[n] for n in x2]
   
    Sigma11, Sigma12, Sigma21, Sigma22 = init_sigma(x1,x2, Sigma_Ui)

    mubar, sigmabar = get_mubar_sigmamu(Sigma_Ui, Ui, x1, Sigma11, Sigma12, Sigma21, Sigma22, mu1, mu2)
    return mubar, sigmabar
end

sigma = 1.0
rho = 0.5
N = 4
p = -0.5
init_item = 1
cov_mat = cov_mat_fun(sigma, rho, N)
Nset = collect(1:N)
mus = zeros(N)
mubar, sigmabar = update_Ui([init_item], [p, p, p, p], mus, cov_mat, Nset)
print(mubar, sigmabar)


