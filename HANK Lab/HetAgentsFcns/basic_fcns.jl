#----------------------------------------------------------------------------
# Basic Functions: Return on capital, marginal utility and its inverse
#----------------------------------------------------------------------------
# mutil(c::Array, ξ::Float64)     = 1.0 ./ c.^(ξ)
# invmutil(mu::Array, ξ::Float64)  = 1.0 ./ mu.^(1.0/ξ)
mutil(c::Array,#Union{Array{Float64},Array{DualNumbers.Dual{Float64}}},
    ξ::Float64)      = 1.0 ./ c.^ξ#((c.*c).*(c.*c))#
invmutil(mu::Array,#Union{Array{Float64},Array{DualNumbers.Dual{Float64}}},
    ξ::Float64)  = 1.0 ./ mu.^(1.0./ξ)#(sqrt.(sqrt.(mu)))#

# Incomes (K:capital, A: TFP): Interst rate = MPK.-δ, Wage = MPL, profits = Y-wL-(r+\delta)*K
interest(K::Number, A::Number,N::Number, m_par::ModelParameters) = A.* m_par.α .* (K./N) .^(m_par.α - 1.0) .- m_par.δ_0         # A=TFP * MC
wage(K::Number, A::Number,N::Number, m_par::ModelParameters)     = A.* (1-m_par.α) .* (K./N) .^m_par.α # A=TFP * MC
output(K::Number, A::Number,N::Number, m_par::ModelParameters)  = A.* K .^(m_par.α).*N .^(1-m_par.α)
employment(K::Number, A::Number, m_par::ModelParameters)     = (A.* (1 .- m_par.τ_lev).*(1.0-m_par.α) .* K .^(m_par.α )).^(1.0./(m_par.γ+m_par.α)) # A=TFP*MC
#employment(K::Number, A::Number, m_par::ModelParameters)     = (A.* (1.0-m_par.α) .* K .^(m_par.α )).^(1.0/(m_par.γ+m_par.α)) # A=TFP*MC

# (A.* (1.0-m_par.α) .* K .^(m_par.α )).^(1.0/(m_par.γ+m_par.α)) # A=TFP*MC

#.^(1/(1-par.alpha+par.gamma))
# Ns^1/γ= w * (1-τ)*(1-mcw)
# w = mc*Z *(1-α)*K^α*N^(-α)


##### Distributional summaries ##########################
function distrSummaries(distr::AbstractArray,c_a_star::AbstractArray,
                        c_n_star::AbstractArray, n_par::NumericalParameters,
                        inc::AbstractArray, m_par::ModelParameters)
    ## Distributional summaries
    mplusk = zeros(n_par.nk.*n_par.nm)
    for k = 1:n_par.nk
        for m = 1:n_par.nm
            mplusk[m+(k-1)*n_par.nm] = n_par.grid_m[m]+n_par.grid_k[k];
        end
    end
    IX=sortperm(mplusk)
    mplusk          = mplusk[IX]
    moneycapital_pdf= sum(distr, dims=3)
    moneycapital_pdf= moneycapital_pdf[IX];
    moneycapital_cdf= cumsum(moneycapital_pdf);
    S               = [0; cumsum(moneycapital_pdf.*mplusk)]
    giniwealth      = 1-(sum(moneycapital_pdf.*(S[1:end-1]+S[2:end]))/S[end]);

    distr_m = sum(distr,dims=(2,3))[:]
    distr_k = sum(distr,dims=(1,3))[:]
    distr_y = sum(distr,dims=(1,2))[:]

    share_borrower = sum(distr_m.*(n_par.grid_m.<0))

    p50             = count(moneycapital_cdf.<0.5)+1;
    p90             = count(moneycapital_cdf.<0.9)+1;
    w9050           = mplusk[p90]./mplusk[p50];
    w90share        = sum(mplusk[p90:end].*moneycapital_pdf[p90:end])./sum(mplusk.*moneycapital_pdf);
    #sdlogw=sqrt(moneycapital_pdf[:]'*log.(mplusk[:]).^2-(moneycapital_pdf[:]'*log.(mplusk[:]))^2);
    x               = zeros(eltype(c_a_star),(n_par.nm, n_par.nk,n_par.ny,2))
    c               = zeros(eltype(c_a_star),(n_par.nm, n_par.nk,n_par.ny,2))
    distr_x         = zeros(eltype(c_a_star),(n_par.nm, n_par.nk,n_par.ny,2))
    x[:,:,:,1]      = c_a_star
    x[:,:,:,2]      = c_n_star;
    aux_x           = inc[5]#(1 .- m_par.τ_bar).*(1.0 ./ m_par.μw).*w.*N.*n_par.mesh_y./(m_par.γ+1);
    aux_x[:,:,end]  = zeros(n_par.nm,n_par.nk);
    c[:,:,:,1]      = x[:,:,:,1] +aux_x;
    c[:,:,:,2]      = x[:,:,:,2] +aux_x;
    distr_x[:,:,:,1] = m_par.λ .* distr
    distr_x[:,:,:,2] = (1-m_par.λ) .* distr

    IX              = sortperm(x[:]);
    x               = x[IX];
    x_pdf           = distr_x[IX];
    S               = cumsum(x_pdf.*x);
    S               = [0 S'];
    ginicompconsumption = 1-(sum(x_pdf.*(S[1:end-1]+S[2:end]))/S[end]);
    sdlogx          = sqrt(x_pdf[:]'*log.(x[:]).^2-(x_pdf[:]'*log.(x[:]))^2);

    IX              = sortperm(c[:]);
    c               = c[IX];
    c_pdf           = distr_x[IX];
    S               = cumsum(c_pdf.*c);
    c_cdf           = cumsum(c_pdf);

    p10             = count(c_cdf.<0.1)+1;
    p50             = count(c_cdf.<0.5)+1;
    p90             = count(c_cdf.<0.9)+1;

    c9010           = c[p90]./c[p10];

    p10C            = c[p10];
    p50C            = c[p50];
    p90C            = c[p90];
    # P0020C          = c[1:p20].*(c_pdf[1:p20]./sum(c_pdf[1:p20]))

    S               = [0 S'];
    giniconsumption = 1-(sum(c_pdf.*(S[1:end-1]+S[2:end]))/S[end]);
    sdlogc          = sqrt(c_pdf[:]'*log.(c[:]).^2-(c_pdf[:]'*log.(c[:]))^2);

    Yidio           = inc[6]+inc[2]+inc[3] - n_par.mesh_m
    IX              = sortperm(Yidio[:])
    Yidio           = Yidio[IX]
    Y_pdf           = distr[IX]
    Y_cdf           = cumsum(Y_pdf)
    p10             = count(Y_cdf.<0.1)+1
    p90             = count(Y_cdf.<0.9)+1
    y9010           = Yidio[p90]./Yidio[p10]
    I90share     = sum(Yidio[p90:end].*Y_pdf[p90:end])./sum(Yidio.*Y_pdf);


    S               = cumsum(Y_pdf.*Yidio)
    S               = [0 S']
    giniincome      = 1-(sum(Y_pdf.*(S[1:end-1]+S[2:end]))/S[end])
    sdlogy          = 1.0#sqrt(Y_pdf[:]'*log.(Yidio[:]).^2-(Y_pdf[:]'*log.(Yidio[:]))^2);

    return     distr_m, distr_k, distr_y, share_borrower, giniwealth, I90share, ginicompconsumption,#=
            =# sdlogx, c9010, giniconsumption, sdlogc, y9010, giniincome, sdlogy, w90share, p10C, p50C, p90C
end
