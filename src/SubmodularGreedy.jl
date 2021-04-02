module SubmodularGreedy

using DataStructures
using StatsBase
using Revise

include("helper-funs.jl")
export intersection_ind_oracle

include("simultaneous-greedys.jl")
export simultaneous_greedys

include("sample-greedy.jl")
export greedy, sample_greedy

include("repeated-greedy.jl")
export repeated_greedy, deterministic_usm

include("test-instance.jl")
export card_add_ind, lin_plus_sqrt_card_diff, weighted_cut_diff

end # module
