
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
#
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


###
if n_par.baseline == true
    x0      = zeros(size(LOMstate,1), 1)
    x0[indexes.Sshock] = 0.54 # m_par.sigmaS

    MX = [I; State2Control]
    nlag=40
    x = x0*ones(1,nlag+1)
    IRF_state_sparse = zeros(indexes.profits, nlag)
    for t = 1:nlag
            IRF_state_sparse[:, t]= (MX*x[:,t])'
            x[:,t+1] = LOMstate * x[:,t]
    end

    using DataFrames
    df = DataFrame(Y = IRF_state_sparse[indexes.Y,:]*100,
                   C = IRF_state_sparse[indexes.C,:]*100,
                   I = IRF_state_sparse[indexes.I,:]*100,

                   )
    print(df)

    using CSV
    CSV.write("export_unc_irf.csv", df)

end

# ------------------------------------------------------------------------------
#   Next we define the observed variables, their transformation (level or growth rate),
## STEP 3: Estimation - Mode Finding
#   and whether they are measured with error. We then load the data and construct
#   the observation equation. Finally, everything is passed to the optimizer
# ------------------------------------------------------------------------------
rng = MersenneTwister(1234);
nlag=1000


TFPshock = 0.007.*randn!(rng, zeros(nlag))
MPshock  = 0.0015.*randn!(rng, zeros(nlag))
UNCshock = 0.054.*randn!(rng, zeros(nlag))



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




## ###
function hodrick_prescott_filter(y, lambda)
    T = size(y, 1)
    matrix = zeros(T, T)

    matrix[1, 1:3] = [1 + lambda, -2 * lambda, lambda]
    matrix[2, 1:4] = [-2 * lambda, 1 + 5 * lambda, -4 * lambda, lambda]

    for i = 3 : T - 2
        matrix[i, i-2 : i+2] = [lambda, -4*lambda, 1 + 6 * lambda, -4 * lambda, lambda]
    end

    matrix[T-1, T-3:T] = [lambda, -4 * lambda, 1 + 5 * lambda, -2 * lambda]
    matrix[T, T-2:T] = [lambda, -2 * lambda, 1 + lambda]

    trend = matrix \ y
    cycle = y - trend

    return trend, cycle
end;

YT=exp.(XSS[indexes.YSS] .+ IRF_state_sparse[indexes.Y,:])
CT=exp.(XSS[indexes.CSS] .+ IRF_state_sparse[indexes.C,:])
IT=exp.(XSS[indexes.ISS] .+ IRF_state_sparse[indexes.I,:])

YTR, YC = hodrick_prescott_filter(log.(YT), 1600)
CTR, CC = hodrick_prescott_filter(log.(CT), 1600)
ITR, IC = hodrick_prescott_filter(log.(IT), 1600)

BCM=zeros(4,1)
BCM[1]=100*std(YC)
BCM[2]=100*std(CC)
BCM[3]=100*std(IC)

BCMNOHP=zeros(4,1)
BCMNOHP[1]=100*std(IRF_state_sparse[indexes.Y,:])
BCMNOHP[2]=100*std(IRF_state_sparse[indexes.C,:])
BCMNOHP[3]=100*std(IRF_state_sparse[indexes.I,:])

QT=exp.(IRF_state_sparse[indexes.q,:])
PIT=exp.(IRF_state_sparse[indexes.π,:])
rT=exp.(XSS[indexes.rSS] .+ IRF_state_sparse[indexes.r,:]) .- 1.0
RBT=exp.(XSS[indexes.RBSS] .+ IRF_state_sparse[indexes.RB,:])

LP=100*100*PIT[2:end].*(((QT[2:end]+rT[2:end])./QT[1:end-1])-(RBT[2:end]./PIT[2:end]));

BCM[4] = mean(LP)./std(LP)


using DelimitedFiles

writedlm(save_file1, BCM)


## #
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
F  = FSYSNODUAL(X0,X0,XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC,DCD,IDCD)

FT=zeros(eltype(X0),length_X0,nlag)
BT=zeros(eltype(X0),1,nlag)



FT[:,1]  = FSYSNODUAL(X0,IRF_state_sparse[:,1] ,XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC,DCD,IDCD)

for t = 1:nlag

    try

    FT[:,t] = FSYSNODUAL(IRF_state_sparse[:,t],IRF_state_sparse[:,t+1],XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC,DCD,IDCD)

    catch

    end


end

## #
FTERROR=zeros(2,2)

FTERROR[1,1]= mean(abs.(FT[indexes.B,:]))*100
FTERROR[1,2]= maximum(abs.(FT[indexes.B,:]))*100
FTERROR[2,1]= mean(abs.(FT[indexes.K,:]))*100
FTERROR[2,2]= maximum(abs.(FT[indexes.K,:]))*100

writedlm(save_file2, FTERROR)
