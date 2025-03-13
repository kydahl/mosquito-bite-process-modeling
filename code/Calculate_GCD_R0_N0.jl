# Numerical calculators for key outputs
# Load packages
using IntervalSets
using LinearAlgebra
using Symbolics
using Serialization
using DataFrames
using CSV
using ProgressBars

# Load parameters
include("Parameters.jl")

# Define output functions

# Gonotrophic cycle duration (GCD) calculator
function GCD_func(B_vals_in)
    (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in

    # one_vec_five = [1.0f0 1.0f0 1.0f0 1.0f0 1.0f0]
    # alpha_vec_five = [1.0f0; 0.0f0; 0.0f0; 0.0f0; 0.0f0]
    one_vec_four = [1.0f0; 1.0f0; 1.0f0; 1.0f0]
    alpha_vec_four = [1.0f0; 0.0f0; 0.0f0; 0.0f0]

    # Define subintensity matrices
    A_mat = [
        -lQ+(1-pQ)*lQ       pQ*lQ               0f0   0f0   
        (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0f0 
        (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
        (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0f0   -lG
    ]

    GCD = transpose(alpha_vec_four) * transpose(A_mat) * one_vec_four

    return(GCD)
end

# Basic offspring number (N0) calculator
function N_offspring_func(B_vals_in)
    (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in

    # one_vec_five = [1.0f0 1.0f0 1.0f0 1.0f0 1.0f0]
    # alpha_vec_five = [1.0f0; 0.0f0; 0.0f0; 0.0f0; 0.0f0]
    one_vec_four = [1.0f0; 1.0f0; 1.0f0; 1.0f0]
    alpha_vec_four = [1.0f0; 0.0f0; 0.0f0; 0.0f0]

    # Define subintensity matrices
    A_mat = [
        -lQ+(1-pQ)*lQ       pQ*lQ               0f0   0f0   
        (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0f0 
        (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
        (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0f0   -lG
    ]

    tau = transpose(-A_mat * one_vec_four) * inv(mu * I - transpose(A_mat)) * alpha_vec_four
    tau = tau[1]
    rho = (gV / (mu + gV)) * (gR / (mu + gR)) * tau
    nG = 1.0f0 / (1.0f0 - rho)
    # Basic offspring number
    N_offspring = tau * (varPhi / (mu + gV)) * (rhoJ / (rhoJ + muJ)) * nG
    return(N_offspring)
end

# Basic reproduction number (R0) calculator
function R0_func(B_vals_in)
    (lQ, lL, lP, lG, sigma, pQ, pL, pP, pG) = B_vals_in

    one_vec_four = [1.0f0; 1.0f0; 1.0f0; 1.0f0]
    alpha_vec_four = [1.0f0; 0.0f0; 0.0f0; 0.0f0]

    # Define subintensity matrices
    A_mat = [
        -lQ+(1-pQ)*lQ       pQ*lQ               0f0   0f0   
        (1-sigma)*(1-pL)*lL -lL+sigma*(1-pL)*lL pL*lL 0f0 
        (1-sigma)*(1-pP)*lP sigma*(1-pP)*lP     -lP   pP*lP
        (1-sigma)*(1-pG)*lG sigma*(1-pG)*lG     0f0   -lG
    ]

    tau = transpose(-A_mat * one_vec_four) * inv(mu * I - transpose(A_mat)) * alpha_vec_four
    tau = tau[1]
    rho = (gV / (mu + gV)) * (gR / (mu + gR)) * tau
    nG = 1.0f0 / (1.0f0 - rho)
    # Basic offspring number
    N_offspring = tau * (varPhi / (mu + gV)) * (rhoJ / (rhoJ + muJ)) * nG
    
    if N_offspring < 1
        return(0)
    else 
        # Distribution of mosquitoes across states at equilibrium
        B_prefactor = rhoJ * KJ *(N_offspring - 1) * nG / N_offspring
        temp_inv = inv(mu * I - transpose(A_mat))
        B_postfactor = temp_inv * alpha_vec_four

        B_star = B_prefactor * B_postfactor
        KB = sum(B_star)

        betaH_mat = zeros(Float64, 4,4); betaV_mat = zeros(Float64, 4,4); LambdaH = zeros(Float64, 4,4); LambdaV = zeros(Float64, 4,4)
        
        LambdaH[3,3] = lP * KB / KH # Ross-Macdonald contact rate assumption
        LambdaV[4,4] = lG
        betaH_mat[3,3] = bH
        betaV_mat[4,4] = bB

        spec_mat = alpha_vec_four * transpose(-A_mat * one_vec_four)

        GammaI = inv(mu * I - transpose(A_mat) + (gV / (mu + gV)) * (gR / (mu + gR)) * spec_mat)
        GammaE = inv((eta + mu) * I - transpose(A_mat) + (gR / (mu + gR + eta)) * (gV / (mu + gV + eta)) * spec_mat)

        complicated_probability = (gV/(mu+gV)) * ((eta / (mu+gV+eta)) * (gR/(mu+gR+eta)) + (eta/(mu+gR+eta) * (gR/(mu+gR))))
        tauE = (eta * I + complicated_probability * spec_mat) * GammaE

        sum_B_star = sum(B_star)

        R02 = (1 / (gH + muH)) * transpose(one_vec_four) * betaH_mat * LambdaH * GammaI * tauE * betaV_mat * LambdaV * (B_star / sum_B_star)
        R0 = sqrt(R02[1])
        return(R0)
    end
end

# Run each function once to pre-load it into Julia

GCD_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
N_offspring_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
temp = R0_func([1.0f0, 1.0f0, 1.0f0, 1.0f0, 0.5f0, 1.0f0, 1.0f0, 1.0f0, 1.0f0])
