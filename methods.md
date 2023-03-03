# Methods

Notes in methods.

- note Elith Kearney etc on range expansion MEE
-


## Data

## WHO Malaria Threats Map

Invasive vector species data were sourced from the [WHO Malaria Threats Map](https://apps.who.int/malaria/maps/threats/)

[data/MTM_INVASIVE_VECTOR_SPECIES_20230303.xlsx]([data/MTM_INVASIVE_VECTOR_SPECIES_20230303.xlsx])

## Spread model 

[R/toy_spread_model.R](R/toy_spread_model.R)

essentially a dynamic range model a la Pagel: https://onlinelibrary.wiley.com/doi/10.1111/j.1466-8238.2011.00663.x but on a network), so I coded up this toy example. It doesn't run, since there are some bits of the toolkit for these models that haven't yet been ported to TF2
