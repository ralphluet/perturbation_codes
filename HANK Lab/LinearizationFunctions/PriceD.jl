function PriceD(QuantGuess,VmPrimeDH,VkPrimeDH,PIPrime,B,K,H,RB,Z,m_par,n_par,JDminus,Xss)
    KPrime=QuantGuess[1];
    PI=QuantGuess[2];


    q=((KPrime/K)-1)*m_par.ϕ + 1.0;

    RBPrime =  exp(Xss[indexes.RBSS] +
    ((1 - m_par.ρ_R) * m_par.θ_π).*log(PI) +
    m_par.ρ_R * (log(RB) - Xss[indexes.RBSS]))

    mc        =  1 ./ m_par.μ

    N       = employment(K, Z ./(m_par.μ*m_par.μw), m_par)

    MPK  = mc .* Z .* m_par.α .* (K ./ N) .^(m_par.α - 1.0) # marginal product of Capital


    w = (wage(K, Z * mc, N, m_par))  # wages (NEW)

    r = (1 + MPK - m_par.δ_0 )

    Y = Z .* N .^(1.0 .- m_par.α) .* K .^m_par.α       # production function


    profits = (Y * (1.0 - mc))     # profits: + price setting profits + investment profits missing, but =0 on the margin

    mcw = 1.0
    u = 1.0

    # Incomes
    eff_int      = ((RB )           .+ (m_par.Rbar .* (n_par.mesh_m.<=0.0))) ./ PI # effective rate (need to check timing below and inflation)
    eff_intPrime = (RBPrime .+ (m_par.Rbar.*(n_par.mesh_m.<=0.0))) ./ PIPrime

    mcw=1.0
    GHHFA=((m_par.γ)/(m_par.γ+1)) # transformation (scaling) for composite good
    inc =[  GHHFA.*(1.0 .- m_par.τ_lev).*((n_par.mesh_y/H) .*mcw.*w.*N).+
    ((1.0 .- mcw).*w.*N).*(1.0 .- m_par.τ_lev),# labor income (NEW)
    (r .- 1.0).* n_par.mesh_k, # rental income
    eff_int .* n_par.mesh_m, # liquid asset Income
    n_par.mesh_k .* q]
    inc[1][:,:,end].= (1.0 .- m_par.τ_lev).*(n_par.mesh_y[:,:,end] .* profits) # profit income net of taxes

    # Calculate optimal policies
    # expected margginal values
    EVkPrime = copy(reshape(VkPrimeDH,(n_par.nm,n_par.nk, n_par.ny)))
    EVmPrime = copy(reshape(VmPrimeDH,(n_par.nm,n_par.nk, n_par.ny)))

    @views @inbounds begin
        for mm = 1:n_par.nm
            EVkPrime[mm,:,:] .= EVkPrime[mm,:,:]*n_par.Π'
            EVmPrime[mm,:,:] .= eff_intPrime[mm,:,:].*(EVmPrime[mm,:,:]*n_par.Π')
        end
    end
    # roughly 20% time
    c_a_star, m_a_star, k_a_star, c_n_star, m_n_star =
    EGM_policyupdate(EVmPrime ,EVkPrime ,q,PI,RB,1.0,inc,n_par,m_par, false) # policy iteration

    # roughly 20% time
    # Error Term on distribution (in levels, states)
    dPrime        = DirectTransition(m_a_star,  m_n_star, k_a_star,JDminus, m_par.λ, n_par.Π, n_par)
    dPrs          = reshape(dPrime,n_par.nm,n_par.nk,n_par.ny)


    T = m_par.τ_lev .* w .* N .+ m_par.τ_lev .* profits
    Btarget = exp(Xss[indexes.BSS])

    BPrime =  B*(B/Btarget)^( - m_par.γ_B) .* (Y./exp(Xss[indexes.YSS]))^m_par.γ_Y .* PI^m_par.γ_π

    BNEW=dot(dPrs,n_par.mesh_m)
    KNEW=dot(dPrs,n_par.mesh_k)

    QuantDiff=[ KNEW-KPrime,BNEW-BPrime];


    return QuantDiff, RBPrime, dPrs

end
