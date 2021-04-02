# sample-greedy.jl
# Chris Harshaw, October 2020
#
# This file contains functions for greedy and sample greedy algorithms.
#

using StatsBase
# include("simultaneous-greedys.jl")

"""
    greedy([gnd, pq], f_diff, ind_add_oracle; knapsack_constraints=nothing, density_ratio=0.0)

The greedy algorithm, possibly with density ratio thresholding for knapsack constraints.

# Arguments
- `gnd`: a integer array of elements, or a priority queue of elements
- `f_diff`: a marginal difference oracle
- `ind_add_oracle`: an independence oracle

# Optional Arguments
- `knapsack_constraints`: a 2D array of knapsack constraints (default: nothing)
- `density_ratio`: the maximum density ratio threshold (default: 0.0)
- `epsilon`: a number in [0,1] which controls the approximation / speed trade off. Set to `0.0` for exact algorithm.
- `opt_size_ub`: an upper bound on the size of the optimal set

# Output 
- `best_sol`: the best solution 
- `best_f_val`: value of the objective at the solution
- `num_fun`: the number of function evaluations
- `num_oracle`: the number of independence oracle evaluations
- `knap_reject`: an indicator whether the knapsack feasibility was rejected or not
"""
function greedy(pq::PriorityQueue, f_diff, ind_add_oracle; knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, density_ratio::AbstractFloat=0.0, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(pq), verbose::Bool=false)
    best_sol, best_f_val, num_fun, num_oracle, knap_reject = simultaneous_greedy_alg(pq, 1, epsilon, f_diff, ind_add_oracle, knapsack_constraints, density_ratio, opt_size_ub, verbose=verbose)
    return best_sol, best_f_val, num_fun, num_oracle, knap_reject
end

function greedy(gnd::Array{<:Integer}, f_diff, ind_add_oracle; knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, density_ratio::AbstractFloat=0.0, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(gnd), verbose::Bool=false)
    
    # initialize a priority queue, initialize algorithm parameters
    pq, num_fun, num_oracle = initialize_pq(gnd, f_diff, ind_add_oracle, 1, knapsack_constraints)

    # run greedy algorithm
    best_sol, best_f_val, num_f, num_or, knap_reject = greedy(pq, f_diff, ind_add_oracle; knapsack_constraints=knapsack_constraints, density_ratio=density_ratio, epsilon=epsilon, opt_size_ub=opt_size_ub, verbose=verbose)

    # update the number of oracle queries
    num_fun += num_f
    num_oracle += num_or

    return best_sol, best_f_val, num_fun, num_oracle, knap_reject
end

"""
    sample_greedy(gnd, sample_prob, f_diff, ind_add_oracle)

Independently subsample elements from the ground set, then run the greedy algorithm.

# Arguments
- `gnd`: a integer array of elements
- `sample_prob`: the probability to sample an element of the ground set
- `f_diff`: a marginal difference oracle
- `ind_add_oracle`: an independence oracle

# Output 
- `best_sol`: the best solution 
- `best_f_val`: value of the objective at the solution
- `num_fun`: the number of function evaluations
- `num_oracle`: the number of independence oracle evaluations
"""
function sample_greedy(gnd::Array{<:Integer}, sample_prob::AbstractFloat, f_diff, ind_add_oracle; verbose::Bool=false)

    @assert 0 <= sample_prob <= 1

    # subsample the elements
    sample_gnd = Int64[]
    for elm in gnd 
        if rand() <= sample_prob
            push!(sample_gnd, elm)
        end
    end

    if verbose
        print("The subsampled elements are ")
        print(sample_gnd)
    end

    # run greedy with density ratio of 0.0
    best_sol, best_f_val, num_fun, num_oracle, _ = greedy(sample_gnd, f_diff, ind_add_oracle; knapsack_constraints=nothing, density_ratio=0.0, verbose=verbose)
    return best_sol, best_f_val, num_fun, num_oracle
end