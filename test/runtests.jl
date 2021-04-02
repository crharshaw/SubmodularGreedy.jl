using Test
using SubmodularGreedy

verbose = true

# PRINTING FUNCTIONS 

function myshow(x)
    show(stdout, "text/plain", x)
    println()
end

function print_header(title_str)
    header_str = "=================================================="
    println("\n" * header_str)
    println(title_str)
    println(header_str * "\n")
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

# run the tests!
include("helper-tests.jl")
include("alg-tests.jl")

println("\n\nAll tests pass!")