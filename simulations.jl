using Distributed
using JSON
addprocs(7)
@everywhere using IterTools;
@everywhere using Dates;
@everywhere using Distances;
@everywhere using Random;
@everywhere using Distributions;
@everywhere using LinearAlgebra;

""" the minimum # of hops from node i to node j
"""

@everywhere function hop_distance(
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

@everywhere function cov_mat_fun(sigma::Float64, rho::Float64, N::Int64)::Array{Float64,2}
    cov_mat = zeros(Float64, N, N)
    for i in 1:N 
        for j in 1:N
            dist = hop_distance( i, j , N)
            cov_mat[i,j] = rho^dist
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

@everywhere function certainty_equivalent(
        alpha::Int64, 
        mu, 
        sigma::Array{Float64,2}
    )
    new_mu = mu - (.5 * alpha * diag(sigma).^2) 
    return new_mu
end

### Welfare Functions - Statistic Calculation Functions

@everywhere function w_fun(
        CiT::Array{Int64, 1},
        Ui::Array{Float64,1},
        T::Int64
    )
    w_score = 0.0
    for i in 1:length(CiT)
        w_score = w_score + Ui[CiT[i]]
    end
    return w_score*(T^(-1))
end



"""
    init_sigma
init for bayseian update

"""

@everywhere function init_sigma(x1::Array{Int64,1},
            x2::Array{Int64,1},
            Sigma_Ui::Array{Float64,2}, 
            Cit::Array{Int64,1}, 
            Nit::Array{Int64,1})

    Sigma11::Array{Float64,2} = ones(Float64, length(x1), length(x1))
    Sigma12::Array{Float64,2} = ones(Float64, length(x1), length(x2))
    Sigma21::Array{Float64,2} = ones(Float64, length(x2), length(x1))
    Sigma22::Array{Float64,2} = ones(Float64, length(x2), length(x2))

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

@everywhere function get_mubar_sigmamu(
        Sigma_Ui, 
        Ui, 
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

    #mubar = mu2 + (np.matmul(inner,(a-mu1).T)).T
    mubar = mu2 + inner * (a-mu1)
    
    sigmabar = Sigma22 - (inner * Sigma12)

    mu_new = Ui
    sigma_new = Sigma_Ui

    return mu_new, sigma_new, sigmabar, mubar
end

@everywhere function get_sigma_new_mu_new(
        x2, 
        sigmabar, 
        mu_new, 
        sigma_new, 
        mubar
    )
    for i in 1:length(x2)
        mu_new[x2[i]] = mubar[i,1]
        for j in 1:length(x2)
            sigma_new[x2[i], x2[j]] = sigmabar[i, j]
        end
    end
    return mu_new, sigma_new
end

"""
    update_Ui()

Bayesian Update

# Arguments

"""

@everywhere function update_Ui(
            Cit, 
            Ui, 
            ce_Ui, 
            Sigma_Ui::Array{Float64,2}, 
            Nset::Array{Int64,1}
    )
    
    # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
    # μ_bar = μ_1 + Ε12 Ε22^-1 ( a - μ_2 )  

    #println("update Ui")

    x1 = Cit
    x2 = [n for n in Nset if n ∉ Cit]
    Nit = [n for n in Nset if n ∉ Cit]

    mu1 = convert(Array{Float64,1},[ce_Ui[n] for n in x1])
    mu2 = convert(Array{Float64,1},[ce_Ui[n] for n in x2])
   
    Sigma11, Sigma12, Sigma21, Sigma22 = init_sigma(x1,x2, Sigma_Ui, Cit, Nit)

    mu_new, sigma_new, sigmabar, mubar = get_mubar_sigmamu(Sigma_Ui, Ui, x1, Sigma11, Sigma12, Sigma21, Sigma22, mu1, mu2)
    mu_new, sigma_new =  get_sigma_new_mu_new(x2, sigmabar, mu_new, sigma_new, mubar)
    return mu_new, sigma_new
end

@everywhere function choice_helper(
        Cit::Array{Int64, 1},
        mu, 
        choice_set::Array{Int64, 1}
    )
    cit = choice_set[argmax([mu[i] for i in choice_set])]
    return cit
end 


@everywhere function choice_part(V_i, mu_V_i, Sigma_V_i, V, T, N, Nset, alpha, epsilon, beta)
    C_iT::Array{Int64,1} = []
    R_iT::Array{Int64,1} = []

    for t=1:T
        mu_Vit = mu_V_i
        Sigma_Vit = Sigma_V_i
        if length(C_iT) > 0
            # update beliefs
            mu_Vit, Sigma_Vit = update_Ui(C_iT, V_i, mu_V_i, copy(Sigma_V_i), Nset)
        end
        mu_Uit = mu_Vit + beta * V
        # make choice
        ce_Uit = certainty_equivalent(alpha, mu_Uit, Sigma_Vit) # γ: uncertainty aversion
        choice_set = [n for n in Nset if n  ∉ C_iT]
        c_it = nothing
        if rand() < epsilon
            c_it = rand(choice_set)
        else
            c_it = choice_helper(C_iT,ce_Uit, choice_set)
        end
        r_it = choice_set[argmax([V[i] for i in choice_set])]
        append!(R_iT,r_it)
        append!(C_iT,c_it)
    end

    return C_iT, R_iT
end

@everywhere function choice_omni(U_i,T,N, Nset)
    C_iT::Array{Int64,1} = []
    for t=1:T
        choice_set = [n for n in Nset if n ∉ C_iT]
        c_it = choice_helper(C_iT,U_i, choice_set)
        append!(C_iT, c_it)
    end
    return C_iT
end

@everywhere function choice_ind(U_i, 
			mu_U_i, 
			Sigma_U_i::Array{Float64,2}, 
			T::Int64, 
			N::Int64, 
			Nset, 
			alpha::Int64, 
			epsilon::Float64)

    #println("choice_ind")
    C_iT::Array{Int64,1} = []
    for t=1:T
        #println(t)
        if length(C_iT) > 0
            mu_Uit, Sigma_Uit = update_Ui(C_iT, U_i, mu_U_i, Sigma_U_i, Nset)
        else
            mu_Uit = copy(mu_U_i)'
            Sigma_Uit = copy(Sigma_U_i)
        end
        
        mu_Uit = Array(mu_Uit)
        # make choice
        ce_Uit = certainty_equivalent(alpha, mu_Uit, Sigma_Uit)
        choice_set = [n for n in Nset if n ∉ C_iT]
        c_it = nothing
        if rand() < epsilon
            c_it = rand(choice_set)
        else
            c_it = choice_helper(C_iT,ce_Uit, choice_set)
        end
        append!(C_iT, c_it)
    end
		
    return C_iT
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

@everywhere function simulate(N::Int64,
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

    
    print("iteration: $seed ")
    Random.seed!(seed);

    Nset = [ n for n=1:N]   # set of N items I = {1, ..., N}

    # V = (v_n) n in I aka: common value component v_n in vector form

    # MvNormal(mu, sig) 
    #Construct a multivariate normal distribution with mean mu and covariance represented by sig.
    # https://juliastats.github.io/Distributions.jl/stable/multivariate/#Distributions.MvNormal
    mu_V = zeros(Float64, N)
    V = rand(MvNormal(mu_V, Sigma_V))
    mu_V = mu_V'

    C_pop = Dict( "no_rec"  => [], "omni"  => [], "partial" => [])
    W_pop = Dict( "no_rec"  => [], "omni"  => [], "partial" => [])
    R_pop = Dict( "no_rec"  => [], "omni"  => [], "partial" => [])

    for it_ind=1:nr_ind
        #println("nr_ind: $it_ind")
        # V_i = (v_in) n in I aka: consumer i’s idiosyncratic taste for good n in vector form

        mu_V_i = mu_V_ibar = rand(MvNormal(zeros(Float64, N), Sigma_V_ibar))
        #@show mu_V_i
        #@show Sigma_V_i
        V_i = rand(MvNormal(mu_V_i, Sigma_V_i))

        # Utility in vector form
        U_i = V_i + (beta * V)
        mu_U_i = mu_V_i' + beta * mu_V

        ## NO RECOMMENDATION CASE
        Sigma_U_i = Sigma_V_i + beta^2 * (Sigma_V)

        C_iT = choice_ind(U_i, mu_U_i, Sigma_U_i,T,N, Nset, alpha, epsilon)
        append!(C_pop["no_rec"], C_iT)
        w_val = w_fun(C_iT, U_i, T)
        append!(W_pop["no_rec"], w_val)

        
        ## OMNISCIENT CASE
        C_iT = choice_omni(U_i,T,N, Nset)
        append!(C_pop["omni"], C_iT)
        w_val = w_fun(C_iT, U_i, T)
        append!(W_pop["omni"],  w_val)

 
        ## PARTIAL REC Case
        C_iT, R_iT = choice_part(V_i, copy(mu_V_i), copy(Sigma_V_i), V, T, N, Nset, alpha, epsilon, beta)
        append!(C_pop["partial"], C_iT)
        w_val = w_fun(C_iT, U_i, T)
        append!(W_pop["partial"], w_val)
        append!(R_pop["partial"], R_iT)
  
 
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

N_vals = [200]

T_vals = [20]

# Covariance structure
rho_vals = [0.1, 0.3, 0.5, 0.7, 0.9]

# utility idiosyncratic degree 
beta_vals = [ 0, 2, 10, 0]

# absolute risk aversion
alpha_vals = [0, 1]

# action of the time for random exploration
epsilon_vals = [0, 0.1]

sigma_vals = [0.25, 1]

params = Iterators.product(N_vals, T_vals, rho_vals, beta_vals, sigma_vals, alpha_vals, epsilon_vals)

println(length(collect(params)))


sim_results = Dict()
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

    sim_results[(N, T, rho, beta, sigma, alpha, epsilon, nr_pop, nr_ind)] = @sync @distributed vcat for i= 1:nr_pop
        simulate(N, T,sigma, sigma_i, sigma_ibar, beta, nr_ind, Sigma_V_i,  Sigma_V,  Sigma_V_ibar,  alpha, epsilon, i)
    end
    break
end

#WORKING_DIR = "/Users/guyaridor/Desktop/"
WORKING_DIR = "/home/guyaridor/ExAnteFilterBubble/"
open(string(WORKING_DIR, "new_sim.json"),"w") do f
    JSON.print(f, sim_results)
end
