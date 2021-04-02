# helper-funs.jl
# Chris Harshaw, October 2020
#
# This file contains helper functions for creating priority queues and handling knapsack constraints.
#

using Revise
using DataStructures # for priority queue


"""
    initialize_pq(gnd, f_diff, ind_add_oracle, num_sol, knapsack_constraints)

Initialize a priority queue with keys: (elememnt, solution index, staleness, density) and value: marginal gain
"""
function initialize_pq(gnd::Array{<:Integer}, f_diff, ind_add_oracle, num_sol::Integer, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing})

    # initialize a priority queue with initial marginal gain (and record largest maximum gain)
    #   key: element, solution index, staleness, density
    #   value: f(e|S) evaluate for a previous (possibly current) solution S 

    num_oracle = 0
    num_fun = 0

    pq = PriorityQueue{Tuple{Int64, Int64, Int64, Float64}, Float64}(Base.Order.Reverse)
    emptyset = Set{Int64}()
    for elm in gnd

        # if the element is feasible
        num_oracle += 1
        if ind_add_oracle(elm, emptyset) && feasible_knapsack_elm(elm, knapsack_constraints)

            # get the remaining keys (density & staleness) and value (marginal gain)
            density = get_density(elm, knapsack_constraints)
            prev_size = 0
            gain = f_diff(elm, emptyset)
            num_fun += 1

            # add to the priority queue
            for sol_ind = 1:num_sol
                enqueue!(pq, (elm, sol_ind, prev_size, density), gain)
            end
        end
    end

    return pq, num_oracle, num_fun
end


# KNAPSACK FUNCTIONALITY

function feasible_knapsack_elm(elm::Integer, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing})
    if knapsack_constraints === nothing
        return true
    else
        return all(knapsack_constraints[:,elm] .<= 1.0)
    end
end

function get_density(elm::Integer, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing})
    if knapsack_constraints === nothing
        return 0.0
    else
        return sum(knapsack_constraints[:,elm])
    end
end

function knapsack_feasible_to_add(elm::Integer, set_to_update::Integer, sol_costs::Union{Array{<:AbstractFloat}, Nothing}, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing})
    if knapsack_constraints === nothing 
        return true
    else
        return all(sol_costs[:,set_to_update] + knapsack_constraints[:,elm] .<= 1)
    end
end

function ignorable_knapsack(knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing})
    if knapsack_constraints === nothing
        return true
    else
        m, n_k = size(knapsack_constraints)
        return all(reshape(sum(knapsack_constraints, dims=2), m) .<= 1)
    end
end

function num_knapsack(knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing})
    if knapsack_constraints === nothing
        return 0
    else
        return size(knapsack_constraints)[1]
    end
end

function init_knapsack_costs(num_sol::Integer, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing})
    if knapsack_constraints === nothing
        return nothing
    else
        m,n = size(knapsack_constraints)
        sol_costs = zeros(m,num_sol)
        return sol_costs
    end
end

function update_sol_costs!(elm_to_add::Integer, set_to_update::Integer, sol_costs::Union{Array{<:AbstractFloat}, Nothing}, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing})
    if knapsack_constraints !== nothing
        sol_costs[:, set_to_update] += knapsack_constraints[:, elm_to_add]
    end
end  

function dimension_check(gnd::Array{<:Integer}, knapsack_constraints::Union{Array{<:AbstractFloat,2}, Nothing})
    if knapsack_constraints === nothing
        return true
    else
        n = maximum(gnd)
        m, n_k = size(knapsack_constraints)
        return n <= n_k
    end
end


# INDEPENDENCE CONSTRAINTS


function intersection_ind_oracle(elm::Integer, sol::Set{<:Integer}, ind_list...)

    # iterate through all independence oracles in the list
    for ind_add_oracle in ind_list
        if !ind_add_oracle(elm, sol)
            return false
        end
    end
    return true
end

# PRINTING 

"""
    myshow(x)

Displays things in detail, for example matrices correctly.
"""
function myshow(x)
    show(stdout, "text/plain", x)
    println()
end

function printset(sol::Set{<:Integer})

    # handle the empty set
    if length(sol) == 0
        print("{ }")

    else

        # print a sorted larger set
        sol_list = sort(collect(sol))
        n = length(sol)
        print("{ ")
        for i=1:(n-1)
            elm = sol_list[i]
            print("$elm, ")
        end
        elm = sol_list[n]
        print("$elm }")
    end
end

function printlnset(sol::Set{<:Integer})
    printset(sol::Set{<:Integer})
    println("")
end
