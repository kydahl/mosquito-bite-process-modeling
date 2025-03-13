# mosquito-bite-process-modeling
Code for "Once bitten, twice shy: A modeling framework for incorporating heterogeneous mosquito biting into transmission models"

The scripts below are listed in order of dependency: run the first scripts before the latter ones.

To obtain the data used in the manuscript:

Run `Compare_GCD_R0.r` to obtain values of the basic reproduction number (R0) as a function of gonotrophic cycle duration (GCD) or biting state duration (theta) for each of the five model types: standard, exponential, empirical, phenomenological, mechanistic.

Run `Run_LHS.jl` to obtain the Latin Hypercube samples needed for the PRCC analysis of the mechanistic parameters against the three output variables: gonotrophic cycle duration (GCD), basic offspring number (N0), and basic reproduction number (R0).

To generate the figures in the manuscript:

Run `Generate_Figures.R`. Be sure to set `calculate_PRCCs_bool` to `TRUE` for the initial calculation of the PRCC values.