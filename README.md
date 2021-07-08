# SubmodularGreedy.jl

`SubmodularGreedy.jl` is a Julia package for constrained submodular maximization using greedy-based algorithm.
In particular, the following algorithms (including linear time and knapsack variants) are included:

1. Greedy: `greedy`
2. Sample Greedy: `sample_greedy`
3. Repeated Greedy: `repeated_greedy`
4. Simultaneous Greedys: `simultaneous_greedys`

The algorithms above are presented and analyzed in these papers:

1. Moran Feldman, Christopher Harshaw, and Amin Karbasi. "Simultaneous Greedys: A Swiss Army Knife for Constrained Submodular Maximization", arXiv:2009.13998. 2020.
2. Moran Feldman, Christopher Harshaw, and Amin Karbasi. "Greed is Good: Near-Optimal Submodular Maximization via Greedy Optimization", COLT 2017.

## Installing this package

To install this package, you must first have installed the Julia programming language.
If you have not done this, visit [julialang.org](https://julialang.org/) for instructions on how to do so.

The best way to install this package is using Julia's builtin package manager, `Pkg`. 
SubmodularGreedy is currently an unregistered package and so the command for adding it is slightly different than for registered packages.
We discuss how to do this for Julia versions v1.0 and higher.
1. At the command line, type `julia` to enter an interactive Julia session.
2. Enter the package manager by pressing `]`.
3. To download our SubmodularGreedy package, enter the command `add https://github.com/crharshaw/SubmodularGreedy.jl.git`
4. Exit the package manager by pressing press backspace or ^C.

Now SubmodularGreedy is installed with your version of Julia and you can use it.

## How to use this package

The main functionality of this package is to run various greedy-based optimization methods given user-defined value and independence oracles. 


Here is an example of how to use the function `sample_gs_walk`.

```julia
# import the package
using SubmodularGreedy

# construct the ground set
n = 20
gnd = collect(1:n)

# construct a random linear function (marginal gain oracle)
coeff = rand(n)
f_diff(elm,sol) = elm in sol ? 0 : coeff[elm] 

# construct cardinality constraint
card_limit = 4
ind_add_oracle(elm,sol) = length(union(sol, elm)) <= card_limit

# run the greedy algorithm
greedy(gnd, f_diff, ind_add_oracle)
```
For more algorithms, see the documentation inside the package or access the docstrings by typing `?` followed by the name of the function and pressing `Enter`.
For example, to see the documentation for Repeated Greedy, type `?repeated_greedy`.

For a more detailed tutorial of how to use this package, please see the Jupyter notebook `SubmodularGreedy.jl Tutorial.ipynb` in the `notebooks` directory.

## Support
The development of this package was supported in part by 
an NSF Graduate Research Fellowship (DGE1122492) awarded to Christopher Harshaw.
Moran Feldman was supported by ISF grants no. 1357/16 and 459/20.
Amin Karbasi was supported in part by NSF (IIS- 1845032), ONR (N00014-19-1-2406), and TATA  Sons Private Limited. 