
using Distributions;
using Distances;
using Random;

function hop_distance(
        i::Int64,
        j::Int64,
        N::Int64
    )::Int64
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))
end

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

N = 20
sigma = 0.25
rho = 0.1

Sigma_V = cov_mat_fun(sigma,rho,N)
Nset = [ n for n=1:N]   # set of N items I = {1, ..., N}
mu_V = zeros(Float64, N)
V = rand(MvNormal(mu_V, Sigma_V))

@show V
x = rand(V)
@show x
@show typeof(x)
y = argmax(x)
@show y

x = zeros(Float64, N)
@show Random.rand!(V, x)
@show x
