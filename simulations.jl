using IterTools;
using Dates;
using Distances;
using Random;
using Distributions;
using LinearAlgebra;

"""
    iota(n::Int64, N::Int64)
ι : N → R 2 associates with each index n a point in
        the unit circle, evenly spaced, with ι(n) = (cos(n/N) π, sin(n/N) * π)
# Arguments
- `n::Int64`: point on circle  
- `N::Int64`: # of goods so length of circle 
# Returns
- the cos and sin components as Array{Float64,1}
"""
function iota(n,N)::Array{Float64,1}
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
function cov_mat_fun(sigma::Float64, rho::Float64, N::Int64)::Array{Float64,2}
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
        alpha::Int64, 
        mu::Array{Float64,2}, 
        sigma::Array{Float64,2}
    )
    if size(mu) == (200, 200)
        mu = diag(mu)
    end 

    if size(mu) == (1, 200)
        mu = mu'
    end
    @show size(mu)
    @show size (.5 * alpha * diag(sigma).^2)
    new_mu = mu - (.5 * alpha * diag(sigma).^2) 
    return new_mu
end

### Welfare Functions - Statistic Calculation Functions

function w_fun(
        CiT::Array{Int64, 1},
        Ui::Array{Float64,2}
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
function init_sigma(x1::Array{Int64,1},
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

function get_mubar_sigmamu(
        Sigma_Ui::Array{Float64,2}, 
        Ui::Array{Float64,2}, 
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
    mubar = mu2 + inner * transpose(a-mu1)
    
    sigmabar = Sigma22 - (inner * Sigma12)

    mu_new = Ui
    sigma_new = Sigma_Ui

    return mu_new, sigma_new, sigmabar, mubar
end

function get_sigma_new_mu_new(
        x2::Array{Int64,1}, 
        sigmabar::Array{Float64,2}, 
        mu_new::Array{Float64,2}, 
        sigma_new::Array{Float64,2}, 
        mubar::Array{Float64,2}
    )
    @show size(mu_new)
    for i in 1:length(x2)
        mu_new[x2[i]] = mubar[i,1]
        for j in 1:length(x2)
            sigma_new[x2[i], x2[j]] = sigmabar[i, j]
        end
    end
    @show mubar
    return mu_new, sigma_new
end

"""
    update_Ui()

Bayesian Update

# Arguments

"""
function update_Ui(
            Cit::Array{Int64,1}, 
            Ui::Array{Float64,2}, 
            ce_Ui::LinearAlgebra.Adjoint{Float64,Array{Float64,1}}, 
            Sigma_Ui::Array{Float64,2}, 
            Nset::Array{Int64,1}
    )
    
    # https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions
    # μ_bar = μ_1 + Ε12 Ε22^-1 ( a - μ_2 )  

    println("update Ui")

    x1 = Cit
    x2 = [n for n in Nset if n ∉ Cit]
    Nit = [n for n in Nset if n ∉ Cit]

    mu1 = convert(Array{Float64,1},[ce_Ui[n] for n in x1])
    mu2 = convert(Array{Float64,1},[ce_Ui[n] for n in x2])
   
    Sigma11, Sigma12, Sigma21, Sigma22 = init_sigma(x1,x2, Sigma_Ui, Cit, Nit)

    mu_new, sigma_new, sigmabar, mubar = get_mubar_sigmamu(Sigma_Ui, Ui, x1, Sigma11, Sigma12, Sigma21, Sigma22, mu1, mu2)
    mu_new, sigma_new =  get_sigma_new_mu_new(x2, sigmabar, mu_new, sigma_new, mubar)
    println("sjkdakdjhaskjhk")
    @show size(mu_new)
    return mu_new, sigma_new
end

function choice_helper(
        Cit::Array{Int64, 1},
        mu::Array{Float64, 2}, 
        choice_set::Array{Int64, 1}
    )
    cit = choice_set[argmax([mu[i] for i in choice_set])]
    return cit
end 

"""
    choice_ind()

the no recommendation case

# Arguments

"""
function choice_ind(U_i::Array{Float64,2}, 
			mu_U_i::LinearAlgebra.Adjoint{Float64,Array{Float64,1}}, 
			Sigma_U_i::Array{Float64,2}, 
			T::Int64, 
			N::Int64, 
			Nset, 
			alpha::Int64, 
			epsilon::Float64)

    println("choice_ind")
    C_iT::Array{Int64,1} = []
    for t=1:T
        mu_Uit = mu_U_i
        Sigma_Uit = Sigma_U_i
        if length(C_iT) > 0
            # update beliefs
            mu_Uit, Sigma_Uit = update_Ui(C_iT, U_i, mu_U_i, Sigma_U_i, Nset)
            println("ererere")
            @show size(mu_Uit)
        end
        
        mu_Uit = Array(mu_Uit)
        # make choice
        @show typeof(Array(mu_Uit))
        @show size(mu_Uit)
        ce_Uit = certainty_equivalent(alpha, mu_Uit, Sigma_Uit)
        choice_set = [n for n in Nset if n ∉ C_iT]
        c_it = nothing
        if rand() < epsilon
            c_it = rand(choice_set)
        else
            c_it = choice_helper(C_iT,ce_Uit, choice_set)
        end
        C_iT = append!(C_iT, c_it)
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

    
    print("iteration: $seed ")
    Random.seed!(seed);

    Nset = [ n for n=1:N]   # set of N items I = {1, ..., N}

    # V = (v_n) n in I aka: common value component v_n in vector form

    # MvNormal(mu, sig) 
    #Construct a multivariate normal distribution with mean mu and covariance represented by sig.
    # https://juliastats.github.io/Distributions.jl/stable/multivariate/#Distributions.MvNormal
    mu_V = zeros(Float64, N)
    V = MvNormal(mu_V, Sigma_V).Σ.mat
    mu_V = mu_V'

    C_pop = Dict( "no_rec"  => [], "omni"  => [], "partial" => [])
    W_pop = Dict( "no_rec"  => [], "omni"  => [], "partial" => [])
    R_pop = Dict( "no_rec"  => [], "omni"  => [], "partial" => [])

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
        C_iT = choice_ind(U_i, mu_U_i, Sigma_U_i,T,N, Nset, alpha, epsilon)
        append!(C_pop["no_rec"], C_iT)
        w_val = w_fun(C_iT,U_i)
        append!(W_pop[NO_REC], w_val)

        
    
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
