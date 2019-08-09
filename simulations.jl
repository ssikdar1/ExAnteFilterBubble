using IterTools;
using Dates;

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

    

end
