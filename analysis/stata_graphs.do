import delimited "/Users/guyaridor/Desktop/ExAnteFilterBubble/data/data_for_stata.csv", clear


keep if formatted_regime == "Partial"
gen welfare = pop_welfare_avg
gen diversity = pop_diversity_avg

binscatter diversity welfare
