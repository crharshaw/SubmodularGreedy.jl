# test-instance.jl
# Chris Harshaw, October 2020
#
# This file contains objective functions and constraints
# for simple test instances.
#

using Revise

function lin_plus_sqrt_card_diff(elm::Integer, sol::Set{<:Integer}, c::Array{<:Real}, alpha::AbstractFloat)
    if elm in sol
        return 0.0
    else
        k = length(sol)
        return c[elm] + alpha * (sqrt(k+1) - sqrt(k))
    end
end


function card_add_ind(elm::Integer, sol::Set{<:Integer}, k::Integer)
    return length(union(sol, elm)) <= k
end

function weighted_cut_diff(elm::Integer, sol::Set{<:Integer}, A::Array{<:Real, 2})

    # return 0 if the element is in the solution
    if elm in sol
        return 0.0
    end

    # the difference = { weight of edges from e to not S } - { weight of edges from e to S } 
    n = size(A)[2]
    cut_diff = sum([ (i in sol ? -1.0 : 1.0) * A[elm,i] for i in 1:n])
    return cut_diff
end