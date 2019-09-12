import JSON;


OMNI = "omni"
PARTIAL = "partial"
NO_REC = "no_rec"
rec_policy_keys = [OMNI, PARTIAL, NO_REC]

function load()
    s::String = ""
    open("data/new_sim.json") do f
        s = read(f, String)
    end

    j = JSON.parse(s)

    return j
end


function iota(n::Int64,N::Float64)::Array{Float64,1}
    return [cos(n/N) * pi, sin(n/N) * pi]
end

function hop_distance(
        i::Int64,
        j::Int64,
        N::Int64
    )::Int64
    return min(abs(i-j), abs(j-i), abs(j-i-N), abs(i-j-N))
end


function div_fun(CiT, T, N)
    div_score = 0.0
    for i in 1:length(CiT)
        for j in 1:length(CiT)
            if i == j
                continue
            end
            div_score = div_score + euclidean( iota(CiT[i], N), iota(CiT[j], N) )
        end
    end
    return div_score*(1.0/(T*(T-1)))
end


data = load()

stripChar = (s, r) -> replace(s, Regex("[$r]") => "")

INDIVIDUAL_FIELD_NAMES =["pop_idx", "indiv_idx", "regime", "welfare", "diversity_score", "rho", "beta", "epsilon", "follow_recommendation", "N", "T", "sigma", "beta", "alpha", "epsilon", "nr_pop", "nr_ind"]

@show typeof(data)
for (k, value) in data
    for pop_idx in 1:length(value)
        key = split(stripChar(k, "(),"))
        dat = Dict(
            "N" =>  parse(Float64,key[1]),
            "T" =>  parse(Float64, key[2]),
            "rho" =>  key[3],
            "beta" =>  key[4],
            "sigma" =>  key[5],
            "alpha" =>  key[6],
            "epsilon" =>  key[7]
            #"pop_idx" =>  pop_idx 
        )
       
        T = dat["T"]
        N = dat["N"]

        cur = value[pop_idx]
        welfare = cur["Welfare"]
        consumption = cur["Consumption"]
        
        for policy in rec_policy_keys
            dat["regime"] = policy
            for indiv_idx in 1:length(welfare[policy])
                    dat["indiv_idx"] = indiv_idx
                    dat["welfare"] = welfare[policy][indiv_idx]
                    consumption_arr = Array(consumption[policy])
                    dat["diversity_score"] = div_fun(consumption_arr[:,indiv_idx], T, N)
                    dat["follow_recommendation"] = False
                    if policy == PARTIAL
                        follow_rec_arr = Array(cur["Rec"][policy])
                        dat["follow_recommendation"] = follow_rec(consumption_arr[:,indiv_idx], follow_rec_arr[:,indiv_idx], T, N)
                    end
                    cur_dat = copy(dat)
                    data_writer.writerow(cur_dat)
            end
        end

        break
    end
end 
