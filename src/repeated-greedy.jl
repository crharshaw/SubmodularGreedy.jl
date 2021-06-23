# repeated-greedy.jl
# Chris Harshaw, October 2020
#
# This file contains functions for repeated greedy algorithm and its variants.
#

# include("sample-greedy.jl")

"""
    deterministic_usm(gnd, f_diff)

An deterministic algorithm for unconstrained submodular maximization.

# Arguments
- `gnd`: a integer array of elements
- `f_diff`: a marginal difference oracle

# Output 
- `X`: the solution
- `f_val`: objective value of solution
- `num_fun`: the number of function evaluations
"""
function deterministic_usm(gnd::Array{<:Integer}, f_diff; verbose::Bool=false)

    # initialize sets - one empty and one is ground
    Y = Set(gnd)
    X = Set{Int64}()

    # initialize function value and oracle count
    f_val = 0
    num_fun = 0

    # go through each element
    iter = 1
    for elm in gnd 

        # determine marginal gain of adding to X or removing from Y
        a = f_diff(elm, X)
        b = -f_diff(elm, setdiff(Y,elm))
        num_fun += 2

        if verbose
            println("\nIteration $iter looking at element $elm")
            print("Set X is ")
            printlnset(X)
            print("Set Y is ")
            printlnset(Y)
            println("The value a is $a \nThe value b is $b")
        end

        # greedily choose the better option
        if a >= b 
            union!(X,elm)
            f_val += a
        else
            setdiff!(Y, elm)
        end

        iter += 1
    end
    return X, f_val, num_fun
end

"""
    repeated_greedy_alg()

Construct `num_sol` disjoint feasible sets in repeated greedy fashion.

This implementation offers approximate greedy search and lazy evaluations.

# Arguments
- `pq`: a priority queue of elements and marginal gains
- `num_sol`: the number of solutions to sequentially construct 
- `f_diff`: a marginal difference oracle
- `ind_add_oracle`: an independence oracle
- `knapsack_constraints`: a 2D array of knapsack constraints or nothing
- `density_ratio`: the maximum density ratio threshold
- `epsilon`: a number in [0,1] which controls the approximation / speed trade off. Set to `0.0` for exact algorithm.
- `opt_size_ub`: an upper bound on the size of the optimal set

# Output 
- `best_sol`: the best solution 
- `best_f_val`: value of the objective at the solution
- `num_fun`: the number of function evaluations
- `num_oracle`: the number of independence oracle evaluations
- `knap_reject`: an indicator whether the knapsack feasibility was rejected or not
"""
function repeated_greedy_alg(pq::PriorityQueue, num_sol::Integer, f_diff, ind_add_oracle, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, density_ratio::AbstractFloat, epsilon::AbstractFloat, opt_size_ub::Integer; verbose::Bool=false)

    # initialize solution + oracle info
    best_sol = nothing
    best_f_val = -Inf 
    num_fun = 0
    num_oracle = 0
    knap_reject = false

    # repeat the greedy procedure num_sol times
    iter = 1
    for i=1:num_sol

        # run ratio greedy procedure to obtain solution
        sol, f_val, num_f, num_or, greedy_kr = greedy(deepcopy(pq), f_diff, ind_add_oracle; knapsack_constraints=knapsack_constraints, density_ratio=density_ratio, epsilon=epsilon, opt_size_ub=opt_size_ub)

        # update the knapsack reject flag
        knap_reject = knap_reject || greedy_kr

        if verbose
            println("\nIteration $iter")
            print("\tGreedy returned set ")
            printset(sol)
            println(" with value $f_val")
        end

        # update solution + oracle info
        if f_val > best_f_val 
            best_sol = sol 
            best_f_val = f_val
        end
        num_fun += num_f
        num_oracle += num_or

        # remove elements in the greedy solution from the ground set
        keys_to_remove = [k for k in keys(pq) if k[1] in sol]
        for k in keys_to_remove
            delete!(pq, k)
        end

        # run an unconstrained maximization on the given solution
        sol, f_val, num_f = deterministic_usm(collect(sol), f_diff)

        # update solution + oracle info
        if f_val > best_f_val 
            best_sol = sol 
            best_f_val = f_val
        end
        num_fun += num_f

        if verbose
            print("\tUnconstrained returned set ")
            printset(sol)
            println(" with value $f_val")
        end

        # update iteration counter
        iter += 1

        # break if no more elements
        if length(pq) == 0
            break
        end
    end

    return best_sol, best_f_val, num_fun, num_oracle, knap_reject
end

"""
    repeated_greedy_density_search()

A grid search routine for finding a good density ratio.

# Arguments
- `pq`: a priority queue of elements and marginal gains
- `num_sol`: the number of solutions to sequentially construct
- `f_diff`: a marginal difference oracle
- `ind_add_oracle`: an independence oracle
- `knapsack_constraints`: a 2D array of knapsack constraints or nothing 
- `beta_scaling`: scaling term for the density ratio
- `delta`: a number in [0,1] which controls the approximation / speed trade off in the multiplicative search for best density ratio
- `epsilon`: a number in [0,1] which controls the approximation / speed trade off in greedy search
- `opt_size_ub`: an upper bound on the size of the optimal set

# Output 
- `best_sol`: the best solution 
- `best_f_val`: value of the objective at the solution
- `num_fun`: the number of function evaluations
- `num_oracle`: the number of independence oracle evaluations
"""
function repeated_greedy_density_search(pq::PriorityQueue, num_sol::Integer, f_diff, ind_add_oracle, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, beta_scaling::AbstractFloat, delta::AbstractFloat, epsilon::AbstractFloat, opt_size_ub::Integer; verbose=false)

    # get the max gain 
    _, max_gain = peek(pq)

    # initialize solution + oracle info 
    best_sol = nothing
    best_val = -Inf
    num_fun = 0
    num_oracle = 0

    # begin by specifying upper and lower bounds for the density ratio
    lower_density_ratio = beta_scaling * max_gain * 1
    upper_density_ratio = beta_scaling * max_gain * opt_size_ub

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

        # compute the new density ratio, run repeated greedy
        sol, f_val, num_f, num_or, knap_reject = repeated_greedy_alg(deepcopy(pq), num_sol, f_diff, ind_add_oracle, knapsack_constraints, density_ratio, epsilon, opt_size_ub)

        # update solution + oracle info
        if f_val > best_f_val 
            best_sol = sol 
            best_f_val = f_val
        end
        num_fun += num_f
        num_oracle += num_or

       # update the binary search based on knap reject status
       if knap_reject
            upper_density_ratio = density_ratio
        else
            lower_density_ratio = density_ratio
        end

        iter += 1

    end # end grid search

    return best_sol, best_val, num_fun, num_oracle
end

"""
    init_rg_params()

Initializes the best parameters for repeated greedys based on analyses in the literature.

# Arguments
- `num_sol`: number of solutions to construct (might get overridden)
- `k`: parameter of the indepence set, i.e. `k`-system
- `knapsack_constraints`: a 2D array of knapsack constraints or nothing

# Output 
- `num_sol`: the number of solutions to sequentially construct 
- `run_density_search`: bool, `true` if a density search should be run, `false` otherwise
- `beta_scaling`: a scaling parameter for the density ratio when using density search
"""
function init_rg_params(num_sol::Integer, k::Integer, monotone::Bool, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}, epsilon::AbstractFloat)

    # this is the approximation ratio of usm routine
    usm_a = 3.0

    # TODO: this needs to be updated to include usm_a when we do that

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
        beta_scaling = 2 * (1 - epsilon)^2 / (k + 2*m + 1 + usm_a*(num_sol - 1)/2)
        
    else
        if num_sol == 0
            num_sol = Integer( floor(1 + sqrt(2*(k + 2*m + 1) / usm_a) )) 
        end
        beta_scaling = 2 * (1 - epsilon) * (1 - 1/num_sol - epsilon) / (k + 2*m + 1 + usm_a*(num_sol - 1)/2)
        
    end

    return num_sol, run_density_search, beta_scaling
end

"""
    repeated_greedy(gnd, f_diff, ind_add_oracle, k, ...)

A fast implementation of repeated greedy algorithm.

This implementation allows the user to either set some parameters or the code automatically sets them based on worst-case analysis.
For instance, explicitly set the number of solutions by setting the the keyword argument `num_sol` directly; otherwise, `num_sol`
will be set automatically based on the parameter `k` which describes the independence system. At least `num_sol` or `k` must be 
supplied as a keyword argument. 

If knapsack constraints are supplied, then one needs to specify `k` for hyper-parameters used in the density ratio search.
In both of these cases, setting `monotone` further improves the automatic setting of parameters.

# Arguments
- `gnd`: a integer array of elements
- `f_diff`: a marginal difference oracle
- `ind_add_oracle`: an independence oracle

# Optional Keyword Arguments
- `num_sol`: the number of solutions to use (if `0`, then this is set automatically based on worst-case analysis)
- `k`: parameter of the indepence set, i.e. `k`-system
- `knapsack_constraints`: a 2D array of knapsack constraints (default: nothing)
- `epsilon`: a number in [0,1] which controls the approximation / speed trade off
- `opt_size_ub`: an upper bound on the size of the optimal set
- `verbose_lvl`: set `0` to silence output, `1` to have mild output of parameters amnd `2` for full algorithm output.

# Output 
- `best_sol`: the best solution 
- `best_f_val`: value of the objective at the solution
- `num_fun`: the number of function evaluations
- `num_oracle`: the number of independence oracle evaluations
"""
function repeated_greedy(gnd::Array{<:Integer}, f_diff, ind_add_oracle; num_sol::Integer=0, k::Integer=0, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing}=nothing, monotone::Bool=false, epsilon::Float64=0.0, opt_size_ub::Integer=maximum(gnd), verbose_lvl::Integer=1)


    # check that inputs make sense
    @assert (num_sol > 0) || (k > 0) "At least num_sol or k need to be specified."
    @assert dimension_check(gnd, knapsack_constraints) "There are more elemnents in knapsack constraints than in the ground set."
    @assert ignorable_knapsack(knapsack_constraints) || (epsilon > 0) "Because knapsack constraints are provided, please specify a non-zero value of epsilon for running the density search."
    @assert ignorable_knapsack(knapsack_constraints) || (k > 0) "Because knapsack constraints are provided, please specify the independence system parameter k for running the density search."
    @assert 0.0 <= epsilon <= 1.0 "Epsilon needs to be set in the range [0, 1]"

    # initialize a priority queue, initialize algorithm parameters
    num_sol, run_density_search, beta_scaling  = init_rg_params(num_sol, k, monotone, knapsack_constraints, epsilon)
    pq, num_fun, num_oracle = initialize_pq(gnd, f_diff, ind_add_oracle, 1, knapsack_constraints)

    # verbosity levels
    info_verbose = verbose_lvl >= 1
    alg_verbose = verbose_lvl >= 2

    if info_verbose

        # print problem info
        n = length(gnd)
        println("Running repeated greedys\n============================")
        println("Ground set has $n elements")
        if k > 0
            println("Independence system is $k-system")
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

        if run_density_search
            println("A grid search for the best density ratio will be run with the parameters:")
            println("\tbeta scaling term is $beta_scaling")
            println("\terror term is $epsilon")
            println("\tbound on largest set is $opt_size_ub")
        end
    end

    # run the algorithms
    if run_density_search
        delta = epsilon
        best_sol, best_val, num_f, num_or = repeated_greedy_density_search(pq, num_sol, f_diff, ind_add_oracle, knapsack_constraints, beta_scaling, delta, 0.0, opt_size_ub; verbose=alg_verbose)

    else
        # run the plain repeated greedy with density ratio = 0
        density_ratio = 0.0
        best_sol, best_val, num_f, num_or = repeated_greedy_alg(pq, num_sol, f_diff, ind_add_oracle, knapsack_constraints, density_ratio, epsilon, opt_size_ub; verbose=alg_verbose)
    end

    # update the number of oracle queries
    num_fun += num_f
    num_oracle += num_or

    # print final info
    if info_verbose
        if length(best_sol) <= 10
            print("\n\nObtained solution S = ")
            printlnset(best_sol)
        end
        println("Obtained solution has value $best_val")
        println("Required $num_fun function evaluations and $num_oracle independence queries")
    end

    return best_sol, best_val, num_fun, num_oracle
end