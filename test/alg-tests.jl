# alg-tests.jl
# Chris Harshaw, October 2020
#
# These are tests of the algorithms.
#
# Certain aspects of the produced solutions are tested with @test statements. 
# However, checking the verbose outputs are often the most helpful diagnostic.



# CREATE TEST INSTANCE 1 -- AN ALMOST LINEAR FUNCTION
# Objective is f(S) = c(S) + a * sqrt( |S| ), where c is a modular function.
# This objective is monotone submodular and when a is small, it is basically a modular function.

if verbose
    print_header("CREATING TEST INSTANCE 1: ALMOST LINEAR OBJECTIVE + CARDINALITY CONSTRAINT")
end

# parameters
n = 10 # size of ground set
k = 1 # k=1-extendible
kc = 3 # cardinality limit

# create the objective function
c = sort(collect(1:n), rev=true) / n + ones(n)
alpha = 1e-5
f_diff(elm, sol) = lin_plus_sqrt_card_diff(elm, sol, c, alpha)

# create the independence oracle and ground set list 
ind_add_oracle(elm, sol) = card_add_ind(elm, sol, kc)
gnd = collect(1:n)

# define the knapsack constraints 
knapsack_constraints = (1/kc) * ones(2, n)
knapsack_constraints[1,1] = 1 - 1/(2*kc)
knapsack_constraints[2,2] = 1 - 1/(2*kc)

if verbose
    println("The test instance has $n elements in the ground set.")
    println("The independence system is a cardinality consraint of size $kc, which is $k-extendible.")
    println("\nHere are the two knapsack constraints:")
    myshow(knapsack_constraints)
end


# CREATE TEST INSTANCE 2 -- A GRAPH CUT FUNCTION
# Objective is a graph cut on the complete graph, which is f(S) = (n - |S|) * (|S|)
# This objective is non-monotone submodular but it is also simple because it is a function only of cardinality.

if verbose
    print_header("CREATING TEST INSTANCE 2: GRAPH CUT OBJECTIVE + CARDINALITY CONSTRAINT")
end

# create a graph cut on the complete graph
A = ones(n,n)
for i=1:n
    A[i,i] = 0.0
end
wf_diff(elm, sol) = weighted_cut_diff(elm, sol, A)

if verbose
    println("The test instance has $n elements in the ground set.")
    println("The independence system is a cardinality consraint of size $kc, which is $k-extendible.")
end



# TEST 1 -- simultaneous greedys (no knapsack)
if verbose
    print_header("TEST 1: SIMULTANEOUS GREEDYS (no knapsack)")
end
println("\nTesting algorithm output while changing the number of solutions")

num_sol_vals = [1, 2, 3, 4]
for num_sol in num_sol_vals

    # run the algorithm, specifying number of solutions
    print("\ttesting algorithm with $num_sol solutions... S = ")
    best_sol, best_f_val, num_fun, num_oracle = simultaneous_greedys(gnd, f_diff, ind_add_oracle, num_sol=num_sol, verbose_lvl=0)
    printlnset(best_sol)

    # check that the solution contains exactly 1 of the num_sol consecutive elements
    for i=1:kc
        start_elm = 1 + (i-1)*num_sol
        end_elm = i * num_sol
        num_consecutive_elm = sum([elm in best_sol for elm in start_elm:end_elm])
        @test num_consecutive_elm <= 1
    end
end

# TEST 2 -- simultaneous greedys (knapsack)
if verbose
    print_header("TEST 2: SIMULTANEOUS GREEDYS (knapsack)")
end
best_sol, best_f_val, num_fun, num_oracle = simultaneous_greedys(gnd, f_diff, ind_add_oracle, k=k, knapsack_constraints=knapsack_constraints, epsilon=0.1, extendible=true, verbose_lvl=2)

if verbose
    println("\nTesting that the solution doesn't contain both 1 and 2...")
end
@test length(intersect(best_sol, Set([1,2]))) <= 1


# TEST 3 -- the greedy algorithm (no knapsack)
if verbose
    print_header("TEST 3: GREEDY (without knapsack)")
end
best_sol, best_f_val, num_fun, num_oracle = greedy(gnd, f_diff, ind_add_oracle);

if verbose
    println("Testing that the greedy algorithm returns S = { 1, 2, 3 }...")
end
@test best_sol == Set([1,2,3])


# TEST 4 -- the greedy algorithm (with knapsack)
if verbose
    print_header("TEST 4: GREEDY (knapsack)")
end
best_sol, best_f_val, num_fun, num_oracle = greedy(gnd, f_diff, ind_add_oracle, knapsack_constraints=knapsack_constraints);

if verbose
    println("Testing that the greedy algorithm returns S = { 1 }...")
end
@test best_sol == Set([1])


# TEST 5 -- the sample greedy algorithm
if verbose
    print_header("TEST 5: SAMPLE GREEDY")
end
if verbose
    println("Check that the running looks correct\n--------------------------->\n")
end
sample_prob = 0.5
best_sol, best_f_val, num_fun, num_oracle = sample_greedy(gnd, sample_prob, f_diff, ind_add_oracle, verbose=verbose)
if verbose
    println("Testing that the returned solution is feasible wrt cardinality constraints...")
end
@test length(best_sol) <= kc


# TEST 6 -- the unconstrained maximization algorithm (monotone function)
if verbose
    print_header("TEST 6: UNCONSTRAINED MAXIMIZATION (monotone function)")
end
sol, f_val, num_fun = deterministic_usm(gnd, f_diff)
if verbose
    println("Testing that the returned solution is the entire ground set...")
end
@test sol == Set(gnd)


# TEST 7 -- the unconstrained maximization algorithm (non-monotone function, with print out)
if verbose
    print_header("TEST 7: UNCONSTRAINED MAXIMIZATION (non-monotone function)")
end
sol, f_val, num_fun = deterministic_usm(gnd, wf_diff)
if verbose
    println("Testing that the returned solution has size exactly n/2...")
end
@test length(sol) == div(n,2)


# TEST 8 -- testing the repeated greedy wrapper
if verbose
    print_header("TEST 8: REPEATED GREEDY")
end
best_sol, best_val, num_fun, num_oracle = repeated_greedy(gnd, f_diff, ind_add_oracle, num_sol=4)
if verbose
    println("Testing that the returned solution is { 1, 2, 3}")
end
@test best_sol == Set([1,2,3])