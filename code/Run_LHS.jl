# Set up PRCC analysis
include("Calculate_GCD_R0_N0.jl")

# Parameters for the `flighty' and `persistent' mosquito data sets
# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG)
const base_params_flighty = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.9f0, 1.0f0, 0.5f0, 0.5f0, 0.5f0, 0.5f0]
const base_params_persistent = [1/480f0, 1/10f0, 1/5f0, 1/1f0, 1f0 - 0.66f0, 1.0f0, 0.7f0, 0.8f0, 0.9f0, 0.5f0]
    
# Function to set up uniform parameter ranges
function parameter_setup(baseline_vals, stretch_val) 
    ubs = (1f0 + stretch_val) .* baseline_vals
    lbs = max(0,(1f0 - stretch_val)) .* baseline_vals
    ubs[5:end] = [min(v,1) for v in ubs[5:end]]

    return lbs, ubs
end

persistent_lbs, persistent_ubs = parameter_setup(base_params_persistent, 0.1)
flighty_lbs, flighty_ubs = parameter_setup(base_params_flighty, 0.1)

# Set up LHC sampling
using LatinHypercubeSampling

# (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in
min_lbs = [1/(3*1440.0f0), 1/(3*1440.0f0), 1/(3*1440.0f0), 1/(3*1440.0f0), 0.2f0, 0.2f0, 0.2f0, 0.2f0, 0.2f0, 0.0f0]
max_ubs = [60.0f0, 60.0f0, 60.0f0, 60.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0, 100.0f0]

# Set number of LHC samples
n_samples = 10_000::Int

# Set up initial grid
using QuasiMonteCarlo
max_scaled_plan = QuasiMonteCarlo.sample(n_samples, min_lbs, max_ubs, LatinHypercubeSample())
flighty_scaled_plan = QuasiMonteCarlo.sample(n_samples, flighty_lbs, flighty_ubs, LatinHypercubeSample())
persistent_scaled_plan = QuasiMonteCarlo.sample(n_samples, persistent_lbs, persistent_ubs, LatinHypercubeSample())

# Set up grid of parameter combinations
using Base.Threads
using IterTools

function output_calc(LHS_samples)
    # Prepare to store results
    n_samples = size(LHS_samples)[2]
    GCD_results = Vector{Float64}(undef, n_samples::Int)
    N_offspring_results = Vector{Float64}(undef, n_samples::Int)
    R0_results = Vector{Float64}(undef, n_samples::Int)

    # Evaluate the function across all input combinations in parallel

    @threads for idx in ProgressBar(1:n_samples)#(i, (sigma, lQ, lL, lP, lG, pQ, pL, pP, pG)) in ProgressBar(enumerate(parameter_grid))
        B_vals = LHS_samples[1:9,idx]
        # GCD values
        GCD_results[idx] = GCD_func(B_vals)
        # Basic offspring number values
        curr_N_offspring = N_offspring_func(B_vals)
        N_offspring_results[idx] = curr_N_offspring

        # R0 values
        R0_results[idx] = R0_func(B_vals)

    end
    scaled_plan_df = DataFrame(transpose(LHS_samples), [:lQ, :lL, :lP, :lG, :sigma, :pQ, :pL, :pP, :pG, :dummy])
    output_df = scaled_plan_df
    output_df[!,:GCD] = GCD_results
    output_df[!,:N_offspring] = N_offspring_results
    output_df[!,:R0] = R0_results

    return output_df
end

# Load output function
output_calc(QuasiMonteCarlo.sample(100, min_lbs, max_ubs, LatinHypercubeSample()))

# Get outputs for each type then join them
max_results = output_calc(max_scaled_plan)
# flighty_results = output_calc(flighty_scaled_plan)
# persistent_results = output_calc(persistent_scaled_plan)

max_results[!, :type] .= "max"
# flighty_results[!, :type] .= "flighty"
# persistent_results[!, :type] .= "persistent"

all_results = vcat(max_results) #, flighty_results, persistent_results)

# Save outputs ----
# Parameter grid and outputs
using CodecZlib
open(joinpath(dirname(dirname(pwd())), "data", "LHS_samples.csv.gz"), "w") do io
    gzip_io = GzipCompressorStream(io)
    CSV.write(gzip_io, all_results)
    close(gzip_io)
end