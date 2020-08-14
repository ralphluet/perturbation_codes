# Run three cases 1) baseline, 2) more DCT coefficients, 3) copula perturbation

## Package Calls

using Plots, Distributions, BenchmarkTools, JLD2, DataFrames, ForwardDiff, NLsolve
using SparseArrays, LinearAlgebra, Random, LaTeXStrings
using Arpack: eigs; using SpecialFunctions: erf; using FFTW: dct
using Parameters, Setfield, MCMCChains, StatsPlots, Optim, CSV, OrderedCollections
using Flatten; import Flatten: flattenable
using FieldMetadata; import FieldMetadata: prior, label

# Plots, Distributions, BenchmarkTools, JLD2, DataFrames, ForwardDiff, NLsolve SparseArrays, LinearAlgebra, Random, LaTeXStrings, Arpack, SpecialFunctions FFTW, Parameters, Setfield, MCMCChains, StatsPlots, Optim, CSV, OrderedCollections, Flatten, FieldMetadata


shock_names = [:A, :Z, :ZI, :μ, :μw, :Gshock, :Rshock]

state_names=["A", "Z", "ZI", "RB" , "μ", "σ", "Ylag", "Blag","Btarget", "Zlag", "Tlag",
"Glag", "Ilag", "wlag", "qlag", "Nlag", "Clag", "πlag", "σlag", "rlag", "RBlag",
 "τlevlag", "Gshock", "Tlevshock", "Rshock", "Sshock"]
control_names=["r", "w", "K", "π" , "Y" ,"C", "q",  "N", "mc",
"LP", "T", "I", "B", "BY","TY","G","τlev",
"GiniW", "GiniC", "GiniX", "GiniI", "sdlgC", "P9010C", "I90share",
"P9010I", "w90share", "P10C", "P50C", "P90C", "Ygrowth", "Bgrowth", "Zgrowth",
"Ggrowth", "Igrowth", "wgrowth", "qgrowth", "Ngrowth", "Cgrowth", "πgrowth", "σgrowth",
"τlevgrowth", "rgrowth", "RBgrowth", "Tgrowth", "profits"]

state_names_ascii=[ "A", "Z", "ZI", "RB", "mu", "sigma", "Ylag", "Blag",
"Zlag","Tlag", "Glag", "Ilag", "wlag", "qlag", "Nlag", "Clag", "pilag", "sigmalag",
"rlag", "RBlag", "taulevlag", "Gshock", "Tlevshock", "Rshock" ,"Sshock"]
control_names_ascii=["r", "w", "K", "pi" , "Y" ,"C", "q",  "N", "mc",
"u", "LP", "T", "I", "B", "BY", "TY", "G", "taulev", "GiniW","GiniC",
"GiniX","GiniI","sdlgC","P9010C","I90share",
"P9010I", "w90share", "P10C", "P50C", "P90C", "Ygrowth", "Bgrowth", "Zgrowth", "Ggrowth",
"Igrowth", "wgrowth", "qgrowth", "Ngrowth", "Cgrowth", "pigrowth", "sigmagrowth",
"taulevgrowth", "rgrowth", "RBgrowth", "Tgrowth", "profits"]

aggr_names = [state_names; control_names]
aggr_names_ascii=[state_names_ascii; control_names_ascii]

distr_names=["GiniW", "GiniC", "GiniX", "GiniI", "sdlgC", "P9010C", "I90share", "P9010I", "w90share", "P10C", "P50C", "P90C"]

n_FD = 5 # number of derivatives to be calculated simultaneously,
# optimal choice depends on CPU and memory of machine


# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
include("Structs_rbc.jl")

e_set = EstimationSettings()
@make_struct IndexStruct state_names control_names
@make_struct_aggr IndexStructAggr aggr_names

include("NumericalBasics.jl")
include("HetAgentsFcns.jl")
include("LinearizationFunctions.jl")

@make_fn produce_indexes state_names control_names
@make_fnaggr produce_indexes_aggr aggr_names
m_par           = ModelParameters( )

n_par           = NumericalParameters( )


## # Baseline

@set! n_par.reduc = 1e-6
@set! n_par.baseline = true
@set! n_par.copula = false

save_file1 = "../results/table6/DHERROR_DCT6.txt"

include("script_rbc.jl")

## # More DCTs

@set! n_par.reduc = 1e-7
@set! n_par.baseline = false
@set! n_par.copula = false

save_file1 = "../results/table6/DHERROR_DCT7.txt"

include("script_rbc.jl")

## # Copula

@set! n_par.reduc = 1e-6
@set! n_par.baseline = false
@set! n_par.copula = true

save_file1 = "../results/table6/DHERROR_DCT6_COP.txt"

include("script_rbc.jl")
