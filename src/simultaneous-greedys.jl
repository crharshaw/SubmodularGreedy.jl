# simultaneous-greedys.jl
# Chris Harshaw, October 2020
#
# This file contains functions for simultaneous greedy algorithm and its variants.
#

# include("helper-funs.jl")

"""
    simultaneous_greedy_alg()

Construct `num_sol` disjoint feasible sets in simultaneous greedy fashion.

This implementation offers approximate greedy search and lazy evaluations.

# Arguments
- `pq`: a priority queue of element / solution pairs and marginal gains
- `num_sol`: the number of solutions to sequentially construct 
- `epsilon`: a number in [0,1] which controls the approximation / speed trade off. Set to `0.0` for exact algorithm.
- `f_diff`: a marginal difference oracle
- `ind_add_oracle`: an independence oracle
- `knapsack_constraints`: a 2D array of knapsack constraints or nothing
- `density_ratio`: the maximum density ratio threshold
- `opt_ub_size`: an upper bound on the size of the optimal set

# Output 
- `best_sol`: the best solution 
- `best_f_val`: value of the objective at the solution
- `num_fun`: the number of function evaluations
- `num_oracle`: the number of independence oracle evaluations
- `knap_reject`: an indicator whether the knapsack feasibility was rejected or not
"""
function simultaneous_greedy_alg(pq::PriorityQueue, num_sol::Integer, epsilon::AbstractFloat, f_diff, ind_add_oracle, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, density_ratio::AbstractFloat, opt_size_ub::Integer; verbose=false)

    @assert num_sol > 0
    @assert 0 <= epsilon <= 1

    # get dimensions
    n = length(pq)

    # initialize the list of k solutions, their function values
    sol_list = [Set{Int64}() for i=1:num_sol]
    sol_vals_list = [0.0 for i=1:num_sol]

    # create solution costs (knapsack) and create densities of each element
    sol_costs = init_knapsack_costs(num_sol, knapsack_constraints)

    # initialize counts of function / oracle queries / knapsack reject status
    num_fun = 0
    num_oracle = 0
    knap_reject = false

    # get the maximum gain element
    best_elm_info, max_gain = peek(pq)

    # set the threshold
    threshold = (1 - epsilon) * max_gain
    min_threshold = epsilon * max_gain / opt_size_ub
    prev_threshold = threshold

    iter = 1
    print_iter = true

    # repeat until the threshold becomes too small, or we run out of elements
    while (threshold > min_threshold) && (length(pq) > 0)

        if verbose && (epsilon > 0.0) && (prev_threshold > threshold)
            println("The updated threshold is $threshold \n")
        end

        # print iteration info
        if verbose && print_iter
            println("\nIteration $iter")
            for i=1:num_sol
                print("\tSet $i is ")
                printset(sol_list[i])
                val = sol_vals_list[i]
                println(" with value $val")
            end
            println("")
            print_iter = false
        end

        # take out the best element / solution pair 
        (elm, sol_ind, prev_size, density), f_gain = dequeue_pair!(pq)
        add_this_elm = false

        if verbose
            println("\tLooking at element $elm with solution $sol_ind, previously queried at $prev_size with gain $f_gain and density $density")
        end

        # if prev_size = |sol|, then this element / solution pair has maximal gain -- so, add it!
        if prev_size == length(sol_list[sol_ind])

            # check that knapsack constraint is satisfied, then update (if knapsack isn't satisfied, we just go back to the pq)
            if knapsack_feasible_to_add(elm, sol_ind, sol_costs, knapsack_constraints)

                if verbose
                    println("\t\tWe will be adding this element $elm to solution $sol_ind")
                end

                # update solution & solution info
                union!(sol_list[sol_ind], elm)
                sol_vals_list[sol_ind] += f_gain
                update_sol_costs!(elm, sol_ind, sol_costs, knapsack_constraints)

                # remove added element from all element / solution pairs
                keys_to_remove = [k for k in keys(pq) if k[1] == elm]
                for k in keys_to_remove
                    delete!(pq, k)
                end

                iter += 1
                print_iter = true
            end

            # if the gain is less than the threshold, update the threshold appropriately (if no  more elements, then continue)
            if (f_gain < threshold) && (length(pq) > 0)
                prev_threshold = threshold
                threshold = min( (1 - epsilon) * threshold, peek(pq)[2])
            end

        # otherwise, we need to re-evaluate f(elm | sol) and decide whether to re-enqueue
        else

            # first test if S + e remains independent (otherwise throw away element)
            num_oracle += 1
            if ind_add_oracle(elm, sol_list[sol_ind])
                
                # evaluate the new marginal gain
                num_fun += 1
                f_gain = f_diff(elm, sol_list[sol_ind])

                # next, only keep considering the element if marginal gain is above the density ratio (otherwise throw away element)
                if (f_gain > density_ratio * density)

                    # re-enqueue this element
                    prev_size = length(sol_list[sol_ind]) # this is the current set size
                    enqueue!(pq, (elm, sol_ind, prev_size, density), f_gain)

                else
                    # record if the density condition is not satisfied
                    knap_reject = true
                end
            end
        end # end decision of whether to add or re-enqueue
    end # end while -- done creating solutions

    # add the single best element to the solution list
    push!(sol_list, Set([best_elm_info[1]]))
    push!(sol_vals_list, max_gain)
    
    # return the best solution
    best_sol_ind = argmax([sol_vals_list])
    best_sol = sol_list[best_sol_ind]
    best_f_val = sol_vals_list[best_sol_ind]

    return best_sol, best_f_val, num_fun, num_oracle, knap_reject
end 

"""
    density_search()

A binary search routine for finding a good density ratio.

# Arguments
- `pq`: a priority queue of elements and marginal gains
- `num_sol`: the number of solutions to sequentially construct 
- `beta_scaling`: scaling term for the density ratio
- `delta`: determines the granularity of the density grid search 
- `f_diff`: a marginal difference oracle
- `ind_add_oracle`: an independence oracle
- `knapsack_constraints`: a 2D array of knapsack constraints or nothing
- `epsilon`: a number in [0,1] which controls the approximation / speed trade off. Set to `0.0` for exact algorithm.
- `opt_ub_size`: an upper bound on the size of the optimal set

# Output 
- `best_sol`: the best solution 
- `best_f_val`: value of the objective at the solution
- `num_fun`: the number of function evaluations
- `num_oracle`: the number of independence oracle evaluations
"""
function density_search(pq::PriorityQueue, num_sol::Integer, beta_scaling::AbstractFloat, delta::AbstractFloat, f_diff, ind_add_oracle, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, epsilon::AbstractFloat, opt_size_ub::Integer; verbose=true)

    # get the max gain 
    _, max_gain = peek(pq)

    # initialize counts of function / oracle queries
    num_fun = 0
    num_oracle = 0

    # begin by specifying upper and lower bounds for the density ratio
    lower_density_ratio = beta_scaling * max_gain * 1
    upper_density_ratio = beta_scaling * max_gain * opt_size_ub

    # initialize best solutions
    best_sol = nothing
    best_val = -Inf

    iter = 1
    while upper_density_ratio > (1+delta) * lower_density_ratio

        if verbose
            println("\nIteration $iter ")
            println("\tUpper density ratio is $upper_density_ratio and lower density ratio is $lower_density_ratio")
            println("\tBest value seen is $best_val")
            println("\tSo far, we have $num_fun function evaluations and $num_oracle independence queries")
        end

        # compute the new density ratio (geometric average)
        density_ratio = sqrt( lower_density_ratio * upper_density_ratio )

        # run the algorithm 
        sol, fval, num_f, num_or, knap_reject = simultaneous_greedy_alg(deepcopy(pq), num_sol, epsilon, f_diff, ind_add_oracle, knapsack_constraints, density_ratio, opt_size_ub)

        if verbose
            println("\tTesting the middle density threshold $density_ratio,")
            print("\tSolution is ")
            printlnset(sol)
            println("\tWith value $fval")
            println("\tand the knap_reject flag was set to $knap_reject")
        end

        # update the best set 
        if fval > best_val 
            best_val = fval 
            best_sol = sol
        end

        # update the number of function and oracle queries 
        num_fun += num_f 
        num_oracle += num_or

        # update the binary search based on knap reject status
        if knap_reject
            upper_density_ratio = density_ratio
        else
            lower_density_ratio = density_ratio
        end

        iter += 1
    end # end binary search
    
    return best_sol, best_val, num_fun, num_oracle
end

"""
    init_sgs_params()

Initializes the best parameters for simultaneous greedys based on our analysis.

# Arguments
- `num_sol`: integer for number of solutions (if `0`, then this is set automatically based on worst-case analysis)
- `k`: parameter of the indepence set (`k`-extendible or `k`-system)
- `extendible`: set true if the indepence system is `k`-extendible
- `monotone`: set true if the objective is monotone 
- `knapsack_constraints`: a 2D array of knapsack constraints or nothing

# Output 
- `num_sol`: the number of solutions to sequentially construct 
- `run_density_search`: bool, `true` if a density search should be run, `false` otherwise
- `beta_scaling`: a scaling parameter for the density ratio when using density search
"""
function init_sgs_params(num_sol::Integer, k::Integer, extendible::Bool, monotone::Bool, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, epsilon::AbstractFloat)

    # test whether all sets are feasible wrt knapsack constraints
    if ignorable_knapsack(knapsack_constraints)
        run_density_search = false
        m = 0
    else
        run_density_search = true
        m = num_knapsack(knapsack_constraints)
        @assert k > 0
    end

    # compute number of solutions and beta_scaling 
    if monotone

        if num_sol == 0
            num_sol = 1
        end
        
        if extendible
            p = max(num_sol-1, k)
        else
            p = k + num_sol - 1
        end
        
        beta_scaling = 2 * (1 - epsilon)^2 / (p + 1 + 2*m)

    else

        if extendible
            M = max( Integer(ceil(sqrt(1 + 2*m))), k )

            # only set num_sol if it was given as 0
            if num_sol == 0
                num_sol = M + 1
            end
            p = M 
        else

            # only set num_sol if it was given as 0
            if num_sol == 0
                num_sol = Integer(floor(2 + sqrt(k + 2*m + 2)))
            end
            p = k + num_sol - 1
        end

        beta_scaling = 2 * (1 - epsilon) * (1 - 1/num_sol - epsilon) / (p + 1 + 2*m)
    end

    return num_sol, run_density_search, beta_scaling
end

"""
    simultaneous_greedys(gnd, f_diff, ind_add_oracle, k, ...)

A fast implementation of simultaneous greedys algorithm using approximate greedy search and lazy evaluations.

This implementation allows the user to either set some parameters or the code automatically sets them based on worst-case analysis.
For instance, explicitly set the number of solutions by setting the the keyword argument `num_sol` directly; otherwise, `num_sol`
will be set automatically based on the parameter `k` which describes the independence system. At least `num_sol` or `k` must be 
supplied as a keyword argument. 

If knapsack constraints are supplied, then one needs to specify `k` for hyper-parameters used in the density ratio search.
In both of these cases, setting `extendible` and `monotone` further improves the automatic setting of parameters.

# Arguments
- `gnd`: a integer array of elements
- `f_diff`: a marginal difference oracle
- `ind_add_oracle`: an independence oracle

# Optional Keyword Arguments
- `num_sol`: the number of solutions to use (if `0`, then this is set automatically based on worst-case analysis)
- `k`: parameter of the indepence set, e.g. `k`-extendible or `k`-system. If set to `0` then this is not used.
- `knapsack_constraints`: a 2D array of knapsack constraints (default: nothing)
- `extendible`: set true if the indepence system is `k`-extendible
- `monotone`: set true if the objective is monotone
- `epsilon`: a number in [0,1] which controls the approximation / speed trade off. Set to `0.0` for exact algorithm.
- `opt_ub_size`: an upper bound on the size of the optimal set
- `verbose_lvl`: set `0` to silence output, `1` to have mild output of parameters and `2` for full algorithm output.

# Output 
- `best_sol`: the best solution 
- `best_f_val`: value of the objective at the solution
- `num_fun`: the number of function evaluations
- `num_oracle`: the number of independence oracle evaluations
"""
function simultaneous_greedys(gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, extendible::Bool=false, monotone::Bool=false, epsilon::AbstractFloat=0.0, opt_size_ub::Integer=length(gnd), verbose_lvl::Integer=1)

    # check that inputs make sense
    @assert (num_sol > 0) || (k > 0) "At least num_sol or k need to be specified."
    @assert dimension_check(gnd, knapsack_constraints) "There are more elemnents in knapsack constraints than in the ground set."
    @assert ignorable_knapsack(knapsack_constraints) || (epsilon > 0) "Because knapsack constraints are provided, please specify a non-zero value of epsilon for running the density search."
    @assert ignorable_knapsack(knapsack_constraints) || (k > 0) "Because knapsack constraints are provided, please specify the independence system parameter k for running the density search."
    @assert 0.0 <= epsilon <= 1.0 "Epsilon needs to be set in the range [0, 1]"

    # initialize a priority queue, initialize algorithm parameters
    num_sol, run_density_search, beta_scaling  = init_sgs_params(num_sol, k, extendible, monotone, knapsack_constraints, epsilon)
    pq, num_fun, num_oracle = initialize_pq(gnd, f_diff, ind_add_oracle, num_sol, knapsack_constraints)

    # verbosity levels
    info_verbose = verbose_lvl >= 1
    alg_verbose = verbose_lvl >= 2

    if info_verbose

        # print problem info
        n = length(gnd)
        println("Running simultaneous greedys\n============================")
        println("Ground set has $n elements")
        non_str = !monotone ? "non-" : ""
        println("Objective function is " * non_str * "monotone")
        if k > 0
            ext_str = extendible ? "extendible " : ""
            println("Independence system is $k-" * ext_str * "system")
        else
            println("The independence system parameter k is not specified")
        end

        if (knapsack_constraints !== nothing) && ignorable_knapsack(knapsack_constraints)
            println("Knapsack constraints are always satisfied and thus will be ignored")
        elseif !ignorable_knapsack(knapsack_constraints)
            m, n_k = size(knapsack_constraints)
            println("There are $m knapsack constraints")
        end

        # print algorithm info
        println("\nConstructing $num_sol solutions")
        
        if epsilon > 0
            println("An faster (approximate) thresholding greedy search will be run with parameters:")
            println("\terror term is $epsilon")
            println("\tbound on largest set is $opt_size_ub")
        end

        if run_density_search
            println("A binary search for the best density ratio will be run with the parameters:")
            println("\tbeta scaling term is $beta_scaling")
            println("\terror term is $epsilon")
            println("\tbound on largest set is $opt_size_ub")
        end
    end

    # run the algorithms
    if run_density_search
        delta = epsilon
        best_sol, best_val, num_f, num_or = density_search(pq, num_sol, beta_scaling, delta, f_diff, ind_add_oracle, knapsack_constraints, epsilon, opt_size_ub; verbose=alg_verbose)

    else
        # run the plain simultaneous greedy with density ratio = 0
        density_ratio = 0.0
        best_sol, best_val, num_f, num_or, _ = simultaneous_greedy_alg(pq, num_sol, epsilon, f_diff, ind_add_oracle, knapsack_constraints, density_ratio, opt_size_ub; verbose=alg_verbose)
    end

    # update the number of oracle queries
    num_fun += num_f 
    num_oracle += num_or

    # print final info
    if info_verbose
        print("The best solution is ")
        println(best_sol)
        if length(best_sol) <= 10
            print("\n\nObtained solution S = ")
            printlnset(best_sol)
        end
        println("Obtained solution has value $best_val")
        println("Required $num_fun function evaluations and $num_oracle independence queries")
    end

    return best_sol, best_val, num_fun, num_oracle
end