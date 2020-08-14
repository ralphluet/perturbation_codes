function updateV(EVk::Array,#{ForwardDiff.Dual{Nothing,Float64,5},3},
                 c_a_star::Array,#{ForwardDiff.Dual{Nothing,Float64,5},3},
                 c_n_star::Array,#{ForwardDiff.Dual{Nothing,Float64,5},3},
                 m_n_star::Array,#{ForwardDiff.Dual{Nothing,Float64,5},3},
                 r,#::Union{Float64,DualNumbers.Dual{Float64}},
                 Q,#::Union{Float64,DualNumbers.Dual{Float64}},
                 λ,#::Union{Float64,DualNumbers.Dual{Float64}},
                 Π::Array)#{ForwardDiff.Dual{Nothing,Float64,5},2})#Union{Array{Float64,2},Array{DualNumbers.Dual{Float64},2}})

    β::Float64 = m_par.β
    n = size(c_n_star)
    ## Update Marginal Value Bonds
    mutil_c_n = mutil(c_n_star,m_par.ξ) # marginal utility at consumption policy no adjustment
    mutil_c_a = mutil(c_a_star,m_par.ξ) # marginal utility at consumption policy adjustment
    Vm   = λ.*mutil_c_a .+ (1.0-λ ).*mutil_c_n # Expected marginal utility at consumption policy (w &w/o adjustment)

    ## Update marginal Value of Capital
    # Linear interpolate expected future marginal value of capital using savings policy
    Vk = Array{eltype(EVk),3}(undef,n) # Initialize m_n-container
    #EVk_aux = reshape(EVk, n)
@inbounds @views begin
    for j::Int = 1:n[3]
        for k::Int = 1:n[2]
            Vk[:,k,j] = mylinearinterpolate(n_par.grid_m, EVk[:,k,j], m_n_star[:,k,j])
        end
    end
end
    Vk        = r.*Vm .+ λ.*Q.*mutil_c_a .+ (1.0.-λ).*β.*Vk # Expected marginal utility at consumption policy (w &w/o adjustment)
    return Vk, Vm
end
