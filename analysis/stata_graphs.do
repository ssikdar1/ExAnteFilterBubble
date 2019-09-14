import delimited "/Users/guyaridor/Desktop/ExAnteFilterBubble/data/data_for_stata.csv", clear


keep if formatted_regime == "Omniscient" & alpha == 5.0
gen welfare = pop_welfare_avg
gen diversity = pop_diversity_avg

binscatter diversity welfare

import delimited "/Users/guyaridor/Desktop/ExAnteFilterBubble/data/time_path_n_200_t_20.csv", clear

keep if rho == 0.0
