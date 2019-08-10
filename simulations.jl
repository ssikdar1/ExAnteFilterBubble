using IterTools;
using Dates;
using Distances;
using Random;
using Distributions;

"""
    iota(n::Int64, N::Int64)
ι : N → R 2 associates with each index n a point in
        the unit circle, evenly spaced, with ι(n) = (cos(n/N) π, sin(n/N) * π)
# Arguments
- `n::Int64`: point on circle  
- `N::Int64`: # of goods so length of circle 
# Returns
- the cos and sin components as Array
"""
function iota(n,N)
    return [cos(n/N) * pi, sin(n/N) * pi]
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
function cov_mat_fun(sigma::Float64, rho::Float64, N::Int64)
    cov_mat = zeros(Float64, N, N)
    for i in 1:N 
        for j in 1:N
            dist = Distances.euclidean( iota(i,N), iota(j,N) )
            cov_mat[i,j] = rho^dist
        end
    end
    return cov_mat
end

"""
    simulate(
        N::Int64
        T::Int64, 
        sigma::Float64,
        sigma_i::Float64, 
        sigma_ibar::Float64,
        beta::Int64, 
        nr_ind::Int64,
        Sigma_V_i::Array{Float64,2},
        Sigma_V::Array{Float64,2},
        Sigma_V_ibar::Array{Float64,2},
        alpha::Int64,
        epsilon::Float64,
        seed::Float64
    )

# Arguments
# Returns
"""
function simulate(N::Int64,
    T::Int64, 
    sigma::Float64,
    sigma_i::Float64, 
    sigma_ibar::Float64,
    beta::Int64, 
    nr_ind::Int64,
    Sigma_V_i::Array{Float64,2},
    Sigma_V::Array{Float64,2},
    Sigma_V_ibar::Array{Float64,2},
    alpha::Int64,
    epsilon::Float64,
    seed::Int64
    )

    
    println("iteration: $seed")
    Random.seed!(seed);

    Nset = [ n for n=0:N]   # set of N items I = {1, ..., N}

    # V = (v_n) n in I aka: common value component v_n in vector form

    # MvNormal(mu, sig) 
    #Construct a multivariate normal distribution with mean mu and covariance represented by sig.
    # https://juliastats.github.io/Distributions.jl/stable/multivariate/#Distributions.MvNormal
    mu_V = zeros(Float64, N)
    V = MvNormal(mu_V, Sigma_V).Σ.mat
    mu_V = mu_V'

    C_pop = Dict( "rec"  => [], "omni"  => [], "partial" => [])
    W_pop = Dict( "rec"  => [], "omni"  => [], "partial" => [])
    R_pop = Dict( "rec"  => [], "omni"  => [], "partial" => [])

    for it_ind=1:nr_ind
        # V_i = (v_in) n in I aka: consumer i’s idiosyncratic taste for good n in vector form
        mu_V_ibar = MvNormal(zeros(Float64, N), Sigma_V_ibar).μ # TODO why is this always 0?
        mu_V_i = mu_V_ibar
        V_i = MvNormal(mu_V_i, Sigma_V_i).Σ.mat
        mu_V_i = mu_V_i'

        # Utility in vector form
        U_i = V_i + (beta * V)
        mu_U_i = mu_V_i + beta * mu_V

        ## NO RECOMMENDATION CASE
        if beta != 0
            Sigma_U_i = Sigma_V_i + beta^2 * (Sigma_V)
        else
            Sigma_U_i = Sigma_V_i
        end

        # TODO 
        #C_iT = choice_ind(U_i,np.copy(mu_U_i), Sigma_U_i,T,N, Nset, alpha, epsilon)
        #C_pop[NO_REC] += [C_iT]
        #w_val = w_fun(C_iT,U_i)
        #W_pop[NO_REC] += [w_val]


        
    
    end

    return Dict( "Consumption" => C_pop, "Welfare" => W_pop, "Rec" => R_pop )
end

#
nr_pop = 50
#
nr_ind = 100
#
sigma_ibar = .1
#
rho_ibar = 0.0

N_vals = [20]

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

params = Iterators.product(N_vals, T_vals, rho_vals, beta_vals, sigma_vals, alpha_vals, epsilon_vals)

println(length(collect(params)))



for (N, T, rho, beta, sigma, alpha, epsilon) in params
    println("STARTING")
    println("N: $N, T: $T, ρ: $rho β: $beta σ: $sigma α: $alpha  ε: $epsilon")
    println(Dates.now())

    sigma_i = sigma

    Sigma_V_i = cov_mat_fun(sigma, rho, N)
    Sigma_V = cov_mat_fun(sigma,rho,N)

    if beta != 0
        Sigma_V = Sigma_V / beta^2
    end

    Sigma_V_ibar = cov_mat_fun(sigma_ibar,rho_ibar,N)

    # TODO parallelism?
    # TODO list comprehension
    for i= 1:nr_pop
        simulate(N, T,sigma, sigma_i, sigma_ibar, beta, nr_ind, Sigma_V_i,  Sigma_V,  Sigma_V_ibar,  alpha, epsilon, i)
    end

    break
end
