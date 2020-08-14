function denHaan(JDminus::Array{Float64,3}, RBminus::AbstractFloat, X::AbstractArray, XPrime::AbstractArray, Xss::Array{Float64,1}, m_par::ModelParameters,
    n_par::NumericalParameters, indexes::IndexStruct, Γ::Array{Array{Float64,2},1},
    compressionIndexes::Array{Array{Int,1},1}, DC::Array{Array{Float64,2},1},
    IDC::Array{Adjoint{Float64,Array{Float64,2}},1}, DCD::Array{Array{Float64,2},1},
    IDCD::Array{Adjoint{Float64,Array{Float64,2}},1})

    # @generate_equations(aggr_names)


    θm      = uncompress(compressionIndexes[1], XPrime[indexes.Vm], DC,IDC, n_par)
    VmPrimeDH = Xss[indexes.VmSS]+ θm

    VmPrimeDH .= (exp.(VmPrimeDH))

    θk      = uncompress(compressionIndexes[2], XPrime[indexes.Vk], DC,IDC, n_par)
    VkPrimeDH = Xss[indexes.VkSS]+  θk
    VkPrimeDH .= (exp.(VkPrimeDH))


    PIPrime = exp(Xss[indexes.πSS] + XPrime[indexes.π])
    KPrime = exp(Xss[indexes.KSS] + XPrime[indexes.K])

    Z = exp(Xss[indexes.ZSS] + X[indexes.Z])

    B=dot(JDminus,n_par.mesh_m)
    K=dot(JDminus,n_par.mesh_k)
    H=n_par.H


    DenHaanSolve(Qg) = PriceD(Qg,VmPrimeDH,VkPrimeDH,PIPrime,B,K,H,RBminus,Z,m_par,n_par,JDminus,Xss)[1];
    QuantGuess=[exp(Xss[indexes.KSS]),1.0];
    # QuantGuess=[KPrime,PIPrime];


    nlstar=nlsolve(DenHaanSolve, QuantGuess,xtol=1e-6)

    QuantDiff, RBPrime, dPrs = PriceD(nlstar.zero,VmPrimeDH,VkPrimeDH,PIPrime,B,K,H,RBminus,Z,m_par,n_par,JDminus,Xss)

    # xstar, fval, iter, distF = broyden(DenHaanSolve,QuantGuess,1e-8,1e-5,500,0.001,0.0001)
#
    # QuantDiff, RBPrime, dPrs = PriceD(xstar,VmPrimeDH,VkPrimeDH,PIPrime,B,K,H,RBminus,Z,m_par,n_par,JDminus,Xss)

return QuantDiff, RBPrime, dPrs
end
