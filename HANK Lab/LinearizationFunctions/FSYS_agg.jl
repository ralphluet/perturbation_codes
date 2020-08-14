function Fsys_agg(X::AbstractArray, XPrime::AbstractArray, Xss::Array{Float64,1},distrSS::AbstractArray, m_par::ModelParameters,
              n_par::NumericalParameters, indexes::Union{IndexStructAggr,IndexStruct})
              # The function call with Duals takes
              # Reserve space for error terms
    F = zeros(eltype(X),size(X))
    ############################################################################
    #            I. Read out argument values                                   #
    ############################################################################

    ############################################################################
    # I.1. Generate code that reads aggregate states/controls
    #      from steady state deviations. Equations take the form of:
    # r       = exp.(Xss[indexes.rSS] .+ X[indexes.r])
    # rPrime  = exp.(Xss[indexes.rSS] .+ XPrime[indexes.r])
    ############################################################################

    @generate_equations(aggr_names)

mcw = 1.0
u = 1.0

    # Elasticities and steepness from target markups for Phillips Curves
    η        = μ / (μ - 1.0) # demand elasticity
    κ        = η * (m_par.κ / m_par.μ) * (m_par.μ - 1.0) # implied steepnes of phillips curve
    # ηw       = μw / (μw - 1.0) # demand elasticity wages
    # κw       = ηw * (m_par.κw / m_par.μw) * (m_par.μw - 1.0) # implied steepnes of wage phillips curve

    # Capital Utilization
    MPK_SS   = exp(Xss[indexes.rSS]) - 1.0 + m_par.δ_0
    δ_1      = MPK_SS
    δ_2      = δ_1 .* m_par.δ_s
    # Auxiliary variables
    Kserv    = K * u # Effective capital
    MPKserv  = mc .* Z .* m_par.α .* (Kserv ./ N) .^(m_par.α - 1.0) # marginal product of Capital
    depr     = m_par.δ_0 + δ_1 * (u - 1.0) + δ_2 / 2.0 * (u - 1.0)^2.0 # depreciation

    Wagesum         = N*w       # Total wages in economy t
    WagesumPrime    = NPrime*wPrime # Total wages in economy t+1
    ############################################################################
    #           III. Error term calculations (i.e. model starts here)          #
    ############################################################################

    ############################################################################
    #           III. 1. Aggregate Part                                         #
    ############################################################################

    #-------- States -----------#
    # Error Term on exogeneous States
    # Shock processes
    F[indexes.Gshock] = log.(GshockPrime) - m_par.ρ_Gshock * log.(Gshock)
    F[indexes.Tlevshock] = log.(TlevshockPrime) - m_par.ρ_τshock * log.(Tlevshock)

    F[indexes.Rshock] = log.(RshockPrime) - m_par.ρ_Rshock * log.(Rshock)
    F[indexes.Sshock] = log.(SshockPrime) - m_par.ρ_Sshock * log.(Sshock)
    # Observables only correlated with

    F[indexes.A]      = log.(APrime) - m_par.ρ_A * log.(A)    # (unobserved) Private bond return fed-funds spread (produces goods out of nothing if negative)
    F[indexes.Z]      = log.(ZPrime) - m_par.ρ_Z * log.(Z)    # TFP
    F[indexes.ZI]     = log.(ZIPrime) - m_par.ρ_ZI * log.(ZI) # Invsestment-good Productivity

    F[indexes.μ]      = log.(μPrime) - m_par.ρ_μ * log.(μ)  - (1.0 - m_par.ρ_μ) * log.(m_par.μ) # Process for markup target

    F[indexes.σ]      = log.(σPrime) - (m_par.ρ_s * log.(σ) + (1.0 - m_par.ρ_s) *
                        m_par.Σ_n * (log(NPrime) - Xss[indexes.NSS]) + log(Sshock)) # Idiosyncratic income risk (contemporaneous reaction to employment)

    F[indexes.LP]  = log.(LP) - (log((qPrime + rPrime - 1.0)/q) - log(RBPrime / πPrime))

    # Endogeneous States
    F[indexes.Ylag] = log(YlagPrime) - log(Y)
    F[indexes.Blag] = log(BlagPrime) - log(B)
    F[indexes.Zlag] = log(ZlagPrime) - log(Z)
    F[indexes.Glag] = log(GlagPrime) - log(G)
    F[indexes.Ilag] = log(IlagPrime) - log(I)
    F[indexes.wlag] = log(wlagPrime) - log(w)
    F[indexes.Tlag] = log(TlagPrime) - log(T)

    F[indexes.qlag] = log(qlagPrime) - log(q)
    F[indexes.Nlag] = log(NlagPrime) - log(N)
    F[indexes.Clag] = log(ClagPrime) - log(C)
    F[indexes.πlag] = log(πlagPrime) - log(π)
    F[indexes.σlag] = log(σlagPrime) - log(σ)
    F[indexes.rlag] = log(rlagPrime) - log(r)
    F[indexes.RBlag] = log(RBlagPrime) - log(RB)
    F[indexes.τlevlag] = log(τlevlagPrime) - log(τlev)

    # Growth rates
    F[indexes.Ygrowth] = log(Ygrowth) - log(Y/Ylag)
    F[indexes.Tgrowth] = log(Tgrowth) - log(T/Tlag)
    F[indexes.Bgrowth] = log(Bgrowth) - log(B/Blag)
    F[indexes.Zgrowth] = log(Zgrowth) - log(Z/Zlag)
    F[indexes.Ggrowth] = log(Ggrowth) - log(G/Glag)
    F[indexes.Igrowth] = log(Igrowth) - log(I/Ilag)
    F[indexes.wgrowth] = log(wgrowth) - log(w/wlag)
    # F[indexes.mcwwgrowth] = log(mcwwgrowth) - log(mcww/mcwwlag)
    # F[indexes.mcwwlag] = log(mcwwlagPrime) - log(mcww)

    F[indexes.qgrowth] = log(qgrowth) - log(q/qlag)
    F[indexes.Ngrowth] = log(Ngrowth) - log(N/Nlag)
    F[indexes.Cgrowth] = log(Cgrowth) - log(C/Clag)
    F[indexes.πgrowth] = log(πgrowth) - log(π/πlag)
    F[indexes.σgrowth] = log(σgrowth) - log(σ/σlag)
    F[indexes.τlevgrowth] = log(τlevgrowth) - log(τlev/τlevlag)
    F[indexes.rgrowth] = log(rgrowth) - log(r/rlag)
    F[indexes.RBgrowth] = log(RBgrowth) - log(RB/RBlag)

    N_GAP = employment(K, 1.0 ./ (m_par.μ*m_par.μw), m_par)
    Y_GAP = output(K,Z,N_GAP, m_par)

    #  Taylor rule and interest rates
    F[indexes.RB]   = log(RBPrime) - Xss[indexes.RBSS] -
                      ((1 - m_par.ρ_R) * m_par.θ_π).*log(π) -
                      ((1 - m_par.ρ_R) * m_par.θ_Y) .* log(Y/Y_GAP) -
                      m_par.ρ_R * (log.(RB) - Xss[indexes.RBSS])  - log(Rshock)# Taylor rule

    # Tax rule

    F[indexes.T]    = log(T) - log(τlev .* w .* N .+ τlev .* profits)

    F[indexes.τlev]       = log(τlev) - m_par.ρ_τ * log(τlevlag)  - (1.0 - m_par.ρ_τ) *(Xss[indexes.τlevSS]) -
                        (1.0 - m_par.ρ_τ) * m_par.γ_Yτ * log(Y/Y_GAP) -
                        (1.0 - m_par.ρ_τ) * m_par.γ_Bτ * (log(B)- log(Btarget))  - log(Tlevshock)


    # --------- Controls ------------
    # Deficit rule
    # F[indexes.π]   = log(B) - log(Btarget) -
                     # (m_par.γ_B * (log(B)- log(Btarget)) + m_par.γ_B * (log(RB)  - Xss[indexes.RBSS]) - (m_par.γ_B + m_par.γ_π) * log(π) - m_par.γ_Y * (log(T) - Xss[indexes.TSS]))
    F[indexes.π]   = log(BgrowthPrime) + m_par.γ_B * (log(B)- log(Btarget))  -
                                                              m_par.γ_Y * (log(T) - Xss[indexes.TSS])  -
                                                              m_par.γ_π * log(π) - log(Gshock)



    F[indexes.G] = log(G) - log(BPrime + T - RB/π*B)

    F[indexes.Btarget] = log.(BtargetPrime) - (m_par.ρ_Btarget * log.(Btarget) + (1.0 - m_par.ρ_Btarget) *
                                                       (Xss[indexes.BtargetSS]))

    KG = 0.0

# if n_par.adjust == 1 # spending adjustment
     # F[indexes.τlev]   = log(τlev) - Xss[indexes.τlevSS]
     # # expenditures
     # F[indexes.G] = log(G) - log(BPrime + T - RB/π*B)
     # # Tax Revenues
     # F[indexes.T]    = log(T) - log(τlev .* w .* N .+ τlev .* profits)

 #   elseif n_par.adjust == 2 # tax adjustment
     # F[indexes.G]   = log(G) - Xss[indexes.GSS]
     # F[indexes.T]   = log(G) - log(BPrime + T - RB/π*B)
     # # Tax Revenues
     # F[indexes.τlev]   = log(T) - log(τlev .* w .* N .+ τlev .* profits)
 #   elseif n_par.adjust == 3 # fund, spending
     # F[indexes.τlev]   = log(τlev) - Xss[indexes.τlevSS]
     # # expenditures
     # KGPrime        =  (BPrime - exp(Xss[indexes.BSS]))./q
     # KG             =  (B - exp(Xss[indexes.BSS]))./qlag
     # F[indexes.G]   = log(G) - log(BPrime + T - RB/π*B - q*(KGPrime - KG) + (r - 1.0)*KG )
     # # Tax Revenues
     # F[indexes.T]   = log(T) -  log(τlev .* w .* N .+ τlev .* profits)
 #   elseif n_par.adjust == 4 # fund, taxes
     # F[indexes.G]   = log(G) - Xss[indexes.GSS]
     # # expenditures
     # KGPrime        =  (BPrime - exp(Xss[indexes.BSS]))./q
     # KG             =  (B - exp(Xss[indexes.BSS]))./qlag
     # F[indexes.T]   = log(G) - log(BPrime + T - RB/π*B - q*(KGPrime - KG) + (r - 1.0)*KG )
     # # Tax Revenues
     # F[indexes.τlev]   = log(T) -  log(τlev .* w .* N .+ τlev .* profits)
 # end


    # Phillips Curve to determine equilibrium markup, output, factor incomes (!!add irrelevant shifters of beta!!)
    F[indexes.mc]   = (log.(π)- Xss[indexes.πSS]) -
                     (m_par.β * ((log.(πPrime) - Xss[indexes.πSS]) .* YPrime ./ Y) + κ *(mc - 1 ./ μ ))

    # Wage Phillips Curve (NEW)
    # F[indexes.mcw]  = (log.(πw)- Xss[indexes.πwSS]) - (κw *(mcw - 1 ./ μw) +
    #                     m_par.β * ((log.(πwPrime) - Xss[indexes.πwSS]) .* WagesumPrime ./ Wagesum))
        # worker's wage= mcw * firm's wage
    # Wage Dynamics (NEW)
    # F[indexes.πw]   = log.(w./wlag) - log.(πw./π) # LOM for real wage (what firms pay)

    # Capital utilisation
    # F[indexes.u]    = MPKserv  -  q * (δ_1 + δ_2 * (u - 1.0)) # Optimality condition for utilization

    # Prices
    F[indexes.r]    = log.(r) - log.(1 + MPKserv * u - q * depr ) # rate on capital

    # F[indexes.mcww] =  log.(mcww) - log.(mcw * w)   # wages

    F[indexes.w]    = log.(w) - log.(wage(Kserv, Z * mc, N, m_par))  # wages (NEW)

    F[indexes.profits] = log.(profits) - log.(Y * (1.0 - mc))     # profits: + price setting profits + investment profits missing, but =0 on the margin


    # F[indexes.q]    = 1.0 - ZI * q * (1.0 - m_par.ϕ / 2.0 * (Igrowth - 1.0)^2.0 - # price of capital investment adjustment costs
    #                    m_par.ϕ * (Igrowth - 1.0) * Igrowth)  -
    #                    m_par.β * ZIPrime * qPrime * m_par.ϕ * (IgrowthPrime - 1.0) * (IgrowthPrime)^2.0
    F[indexes.q]    =  q - (m_par.ϕ*(KPrime ./ K .- 1.0) + 1.0)

    # Aggregate Quantities
    F[indexes.I]   = KPrime .-  K .* (1.0 .- depr)  .- ZI .* I .* (1.0 .- m_par.ϕ ./ 2.0 .* (Igrowth -1.0).^2.0) # Capital accumulation equation
    F[indexes.N]   = log.(N) - log.(( (1.0 - τlev) .* (mcw.*w)).^(1.0 / m_par.γ))              # labor supply (NEW)
    F[indexes.Y]   = log.(Y) - log.(Z .* N .^(1.0 .- m_par.α) .* Kserv .^m_par.α)         # production function
    # F[indexes.C]   = log.(Y .- G .- I .+ (A .- 1.0) .* RB .* B ./ π .- (δ_1 * (u - 1.0) + δ_2 / 2.0 * (u - 1.0)^2.0).*K) .- log(C) # Resource constraint
    F[indexes.C]   = log.(Y .- G .- I) .- log(C) # Resource constraint

    # Error Term on prices/aggregate summary vars (logarithmic, controls), here difference to SS value
    # averages
    F[indexes.K]      = log.(K-KG)     - Xss[indexes.KSS]
    F[indexes.B]      = log.(B)     -  Xss[indexes.BSS]
    F[indexes.BY]     = log.(BY)    -  log.(B/Y)
    F[indexes.TY]     = log.(TY)    -  log.(T/Y)

    return F
end

function tot_dual(x::ForwardDiff.Dual)
    a = sum(ForwardDiff.partials(x,:))
    return a
end
function realpart(x::ForwardDiff.Dual)
    a = ForwardDiff.value(x)
    return a
end
