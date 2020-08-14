
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
## STEP 1: Find stationary equilibrium and run compression
# XSS contains all state and control variable values in stationary equilibrium,
# XSSaggr contains the subset of aggregate state and controls,
# indexes / indexes_aggr contains the indexes that link variable names to XSS entries,
# compressionIndexes contain the indexes of the parameters of the DCT that are retained,
# COPULA stores the Copula of the stationary equilibrium distribution,
# CDF_i stores the marginal CDFs for liquid, illiquid assets and income.
XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, Copula, n_par, m_par, d,
CDF_SS, CDF_m, CDF_k, CDF_y, distrSS  = find_SS(12)

# # Convention: profits is the last control in the list of control variables
ntotal                  = indexes.profits
@set! n_par.ntotal      = ntotal
@set! n_par.ncontrols   = length(compressionIndexes[1]) + length(compressionIndexes[2])  + n_par.naggrcontrols
@set! n_par.nstates     = n_par.nstates + length(compressionIndexes[3])
# ------------------------------------------------------------------------------
#   Next we need to define the equilibrium error function F(X,X',XSS)
#   X=[State;Controls] difference from SS (States = distribution, productivity, etc.,
## STEP 2: Linearization/Perturbation
#   Controls = consumption policy or MUs, prices, aggregate capital stock
# ------------------------------------------------------------------------------

A=zeros(ntotal,ntotal)
B=zeros(ntotal,ntotal)

@set! n_par.LOMstate_save = zeros(n_par.nstates, n_par.nstates)
@set! n_par.State2Control_save = zeros(n_par.ncontrols, n_par.nstates)

State2Control, LOMstate, alarm_sgu, nk, A, B = SGU(XSS, copy(A), copy(B), m_par, n_par, indexes, Copula, compressionIndexes, distrSS; estim = false)

# ------------------------------------------------------------------------------
#   Next we define the observed variables, their transformation (level or growth rate),
## STEP 3: Estimation - Mode Finding
#   and whether they are measured with error. We then load the data and construct
#   the observation equation. Finally, everything is passed to the optimizer
# ------------------------------------------------------------------------------
rng = MersenneTwister(1234);
nlag=1000


TFPshock = 0.007.*randn!(rng, zeros(nlag))
MPshock = 0.0.*randn!(rng, zeros(nlag))
UNCshock = 0.0.*randn!(rng, zeros(nlag))


x0      = zeros(size(LOMstate,1), 1)

MX = [I; State2Control]
x = x0*ones(1,nlag+1)
IRF_state_sparse = zeros(indexes.profits, nlag)
for t = 1:nlag
    x[indexes.Sshock,t]=UNCshock[t]
    x[indexes.Rshock,t]=MPshock[t]
    x[indexes.Z,t]=x[indexes.Z,t]+TFPshock[t]

    IRF_state_sparse[:, t]= (MX*x[:,t])'
    x[:,t+1] = LOMstate * x[:,t]
end

HX = [LOMstate; State2Control]


## #

# DenHaan error


Γ  = shuffleMatrix(distrSS, n_par)
# Matrices for discrete cosine transforms
DC = Array{Array{Float64,2},1}(undef,3)
DC[1]  = mydctmx(n_par.nm)
DC[2]  = mydctmx(n_par.nk)
DC[3]  = mydctmx(n_par.ny)
IDC    = [DC[1]', DC[2]', DC[3]']

DCD = Array{Array{Float64,2},1}(undef,3)
DCD[1]  = mydctmx(n_par.nm-1)
DCD[2]  = mydctmx(n_par.nk-1)
DCD[3]  = mydctmx(n_par.ny-1)
IDCD    = [DCD[1]', DCD[2]', DCD[3]']

length_X0 = indexes.profits # Convention is that profits is the last control
X0 = zeros(length_X0) #.+ ForwardDiff.Dual(0.0,tuple(zeros(n_FD)...))

QuantDiff  = denHaan(distrSS, 1.0, X0,X0 ,XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC,DCD,IDCD)[1]

println(QuantDiff)
## #

Krealized_RC, Brealized_RC, RBminus= denHaanLoop(x,distrSS, XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC,DCD,IDCD)

## ###

BLINEAR = exp.(XSS[indexes.BSS] .+ IRF_state_sparse[indexes.B,:])
KLINEAR = exp.(XSS[indexes.KSS] .+ IRF_state_sparse[indexes.K,:])

DENHAANERROR=zeros(2,2)
DENHAANERROR[1,1] = 100*mean(abs.(log.(Brealized_RC./BLINEAR)))
DENHAANERROR[1,2] = 100*maximum(abs.(log.(Brealized_RC./BLINEAR)))
DENHAANERROR[2,1] = 100*mean(abs.(log.(Krealized_RC./KLINEAR)))
DENHAANERROR[2,2] = 100*maximum(abs.(log.(Krealized_RC./KLINEAR)))
## #

using DelimitedFiles

writedlm(save_file1, DENHAANERROR)
