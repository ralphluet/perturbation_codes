@metadata prior nothing
@metadata label ""
@metadata latex_label L""
@flattenable @prior @latex_label @label @with_kw struct ModelParameters{T} # true/false denotes whether estimated
	# Household preference parameters
	ξ::T = 4.0    | "xi" | L"\xi" | _  | false  # risk aversion
	γ::T = 1.0    | "gamma" | L"\gamma" |  _  | false  # inverse Frisch elasticity
	β::T = 0.98 | "beta" | L"\beta" |  _  | false  # discount factor
	λ::T = 0.065 | "lambda"  | L"\lambda" | _  | false  # adjustment probability

	# individual income process
	ρ_h::T = 0.98   | "rho" | L"\rho" |  _  | false # autocorrelation income shock
	σ_h::T = 0.06   | "sigma" | L"\sigma" |  _  | false # std of income shocks (steady state)
	ι::T = 0.0625   | "iota" | L"\iota" |  _  | false # probability to return worker
	ζ::T = 0.0005 | "zeta" | L"\zeta" |  _  | false # probability to become entrepreneur

	# Technological parameters
	α::T = 0.3  | "alpha" | L"\alpha" |  _  | false # capital share
	δ_0::T = 0.0135  | "delta" | L"\delta" |  _  | false # depreciation rate
	δ_s::T = 10000.0 | "delta_s" | L"\delta_s" | Gamma(gamma_pars(5.0, 2.0^2)...) | false # depreciation rate increase (flex utilization)
	ϕ::T = 11.4   | "phi" | L"\phi" | Gamma(gamma_pars(4.0, 2.0^2)...)    | true # Capital adjustment costs
	μ::T = 1.05 | "mu" | L"\mu" |  _  | false # markup
	κ::T = 0.09  | "kappa" | L"\kappa" | Gamma(gamma_pars(0.1, 0.02^2)...) | true # Price adjustment costs (in terms of Calvo probs.)
	μw::T = 1.0 | "mu_w" | L"\mu_w" |  _  | false # markup
	κw::T = 1000.0 | "kappa_w" | L"\kappa_w" | Gamma(gamma_pars(0.1, 0.02^2)...) | true # Price adjustment costs (in terms of Calvo probs.)

	# Further steady-state parameters
	ψ::T  = 0.1     | "psi" | L"\psi" |  _  | false # steady-state bond to capital ratio
	τ_lev::T  = 0.3 | "tau_lev" | L"\tau^L" |  _  | false # steady-state income tax rate level
	# τ_prog::T  = 0.18 | "tau_pro" | L"\tau^P" |  _  | false # steady-state income tax rate progressivity

	R::T  = 1.0 | "R" |  L"R"  |  _  | false
	K::T  = 40.0 | "K"  | L"K" |  _  | false
	π::T = 1.0.^0.25 | "Pi"  |  L"\pi" |  _  | false # Steady State Inflation
	RB::T = π*(1.025.^0.25) | "RB" | L"RB"  |  _  | false # Nominal Interest Rate
	Rbar::T = (π*(1.11.^0.25) .- 1.0) | "Rbar" |  L"\bar R"  |  _  | false # borrowing wedge in interest rate

	# exogeneous aggregate "shocks"
	ρ_A::T = 0.9 | "rho_A" | L"\rho_A" | Beta(beta_pars(0.5, 0.2^2)...)        | true # Pers. of bond-spread
	σ_A::T = 0.0 | "sigma_A" | L"\sigma_A" | InverseGamma(ig_pars(0.001, 0.02^2)...) | true # Std of bond-spread shock

	ρ_Z::T = 0.95 | "rho_Z" | L"\rho_Z" | Beta(beta_pars(0.5, 0.2^2)...)        | true # Pers. of TFP
	σ_Z::T = 0.0 | "sigma_Z" | L"\sigma_Z" | InverseGamma(ig_pars(0.001, 0.02^2)...) | true # Std of TFP

	ρ_ZI::T = 0.9 | "rho_Psi" | L"\rho_\Psi" | Beta(beta_pars(0.5, 0.2^2)...)        | true # Pers. of TFP
	σ_ZI::T = 0.0 | "sigma_Psi" | L"\sigma_\Psi" | InverseGamma(ig_pars(0.001, 0.02^2)...) | true # Std of TFP

	ρ_μ::T = 0.9 | "rho_mu" | L"\rho_\mu" | Beta(beta_pars(0.5, 0.2^2)...)        | true # Pers. of price markup
	σ_μ::T = 0.0 | "sigma_mu" | L"\sigma_\mu" | InverseGamma(ig_pars(0.001, 0.02^2)...) | true # Std of cost push shock

	ρ_μw::T = 0.9 | "rho_muw" | L"\rho_{\mu w}" | Beta(beta_pars(0.5, 0.2^2)...)        | true # Pers. of wage markup
	σ_μw::T = 0.0 | "sigma_muw" | L"\sigma_{\mu w}" | InverseGamma(ig_pars(0.001, 0.02^2)...) | true # Std of cost push shock

	# income risk
	ρ_s::T = 0.84 | "rho_sigma" | L"\rho_s" | Beta(beta_pars(0.5, 0.2^2)...)       | false # Persistence of idiosyncratic income risk
	σ_Sshock::T = 0.0  | "sigma_Sshock" | L"\sigma_s" | InverseGamma(ig_pars(0.01, 0.02^2)...) | false #0.54 # std of idiosyncratic income risk
	Σ_n::T = 0.0  | "Sigma_n" | L"\Sigma_N" | Normal(0.0, 100.0) | false #-0.1 # reaction of risk to employment

	# monetary policy
	ρ_R::T = 0.8 | "rho_R" | L"\rho_R" | Beta(beta_pars(0.5, 0.2^2)...)        | true # Pers. in Taylor rule
	σ_Rshock::T = 0.0  | "sigma_Rshock" | L"\sigma_R" | InverseGamma(ig_pars(0.001, 0.02^2)...) | true # Std R
	θ_π::T = 2.0  | "theta_pi" | L"\theta_\pi" | Normal(1.7, 0.3)                   | true # Reaction to inflation
	θ_Y::T = 0.0  | "theta_Y" | L"\theta_y" | Normal(0.125, 0.05)              | true # Reaction to inflation

	# fiscal policy
	γ_B::T = 0.5 | "gamma_B" | L"\gamma_B" | Gamma(gamma_pars(0.1, 0.075^2)...)      | true # reaction of deficit to debt
	γ_π::T = -1.5 | "gamma_pi" | L"\gamma_{\pi}" | Normal(0.0, 1.0) | true # reaction of deficit to inflation
	γ_Y::T = -0.50 | "gamma_Y" | L"\gamma_Y" | Normal(0.0, 1.0) | true # reaction of deficit to output
	ρ_Gshock::T = 0.98 | "rho_Gshock" | L"\rho_G" | Beta(beta_pars(0.5, 0.2^2)...) | true # Pers. in structural deficit
	σ_Gshock::T = 0.00 | "sigma_G" | L"\sigma_G" | InverseGamma(ig_pars(0.001, 0.02^2)...) | true # Std G

	ρ_τ::T = 0.0  | "rho_tau" | L"\rho_\tau" | Beta(beta_pars(0.5, 0.2^2)...)        | false # Pers. in tax level
	σ_Tlevshock::T = 0.00   | "sigma_taushock" | L"\sigma_\tau" | InverseGamma(ig_pars(0.001, 0.02^2)...) | false # Std tax level
	γ_Bτ::T = 0.0 | "gamma_Btau" | L"\gamma_B^\tau" | Gamma(gamma_pars(1.0, 2.0^2)...)  | false # reaction of tax level to debt
	# γ_Bτ::T = 0.0 | "gamma_Btau" | L"\gamma_B^\tau" | Normal(0.0, 1.0) | true # reaction of tax level to debt
	γ_Yτ::T = 0.0 | "gamma_Ytau" | L"\gamma_Y^\tau" | Normal(0.0, 1.0) | false # reaction of tax level to output



	# auxiliary shock parameters
	ρ_Rshock::T = 0.5 | "rho_Rshock" | L"\rho_{Rshock}" | Beta(beta_pars(0.5, 0.2^2)...) | false
	ρ_Pshock::T = 1e-8 | "rho_Pshock" | L"\rho_{Pshock}" | Beta(beta_pars(0.5, 0.2^2)...) | false
	ρ_τshock::T = 1e-8 | "rho_taushock" | L"\rho_{\tau shock}" | Beta(beta_pars(0.5, 0.2^2)...) | false
	ρ_Sshock::T = 1e-8 | "rho_Sshock" | L"\rho_{Sshock}" | Beta(beta_pars(0.5, 0.2^2)...) | false

	ρ_Btarget::T = 0.99999 | "rho_Btarget" | L"\rho_Btarget" | Beta(beta_pars(0.5, 0.2^2)...) | false # Pers. in structural deficit
	σ_Btarget::T = 0.0 |"sigma_Btarget"| L"\sigma_{\bar{B}}"| InverseGamma(ig_pars(0.001, 0.02^2)...) | false
end


@with_kw struct NumericalParameters
	# adjust::Int      = 0
	save_file1::String = "BCM_DCT6.txt"
	save_file2::String = "FTERROR_DCT6.txt"
	ny::Int         = 2     # ngrid income
	nk::Int         = 100   # ngrid capital
	nm::Int         = 100
	nstates::Int    = ny + nk + nm + 14
	naggrstates::Int = 16
	naggrcontrols::Int=16
	ncontrols::Int=16
	ntotal::Int = 100
	naggr::Int = 14
	aggr_names::Array{String,1}=["AHHH"]
	kmin::Float64   = 0.0   # gridmin capital
	kmax::Float64   = 1500.0  # gridmax capital
	mmin::Float64   = 0.0   # gridmin capital
	mmax::Float64   = 1500.0  # gridmax capital
	ϵ::Float64      = 1e-14 # precision
	copula::Bool    = false
	baseline::Bool    = false
	useMATLAB::Bool = false # use MATLAB eigs for finding stationary distribution
	use_parallel::Bool = false # parallel loop over income states
	sol_algo::Symbol = :schur # other options: :schur and :lit

	Π::Matrix{Float64}       = [0.9 0.1; 0.1 0.9] # transition matrix income
	grid_y::Array{Float64,1} = [0.5; 1.5]         # income grid
	grid_k::Array{Float64,1} = exp.(range(log(kmin+1.0), stop = log(kmax+1.0), length = nk)) .- 1.0
	grid_m::Array{Float64,1} = exp.(range(0, stop=log(mmax-mmin+1.0), length = nm)) .+ mmin.- 1.0
	mesh_y::Array{Float64,3} = repeat(reshape(grid_y,(1,1,ny)),outer=[nm, nk, 1])
	mesh_m::Array{Float64,3} = repeat(reshape(grid_m,(nm,1,1)),outer=[1, nk, ny])
	mesh_k::Array{Float64,3} = repeat(reshape(grid_k,(1,nk,1)),outer=[nm, 1, ny])
	reduc::Float64           = 1e-6
	bounds_y::Array{Float64,1} = [0.5; 1; 1.5]
	H::Float64           = 1.0
	Asel::Array{Bool,2} = falses(10,10)
	Bsel::Array{Bool,2} = falses(10,10)
	LOMstate_save::Array{Float64,2} = zeros(nstates, nstates)
	State2Control_save::Array{Float64,2} = zeros(ncontrols, nstates)

	ndraws::Int      = 50
	burnin::Int      = 20
	mhscale::Float64 = 1.0
	dist_guess::Array{Float64,3} = ones(nm, nk, ny)/(nm*nk*ny)
end


@with_kw struct EstimationSettings
	shock_names::Array{Symbol, 1} = [:A, :Z, :ZI, :μ, :μw, :Gshock, :Rshock]
	observed_vars_input::Array{Symbol, 1} = [:Ygrowth, :Igrowth, :Cgrowth, :N, :wgrowth, :RB,  :π]

	data_rename::Dict{Symbol,Symbol} = Dict(
	:pi=>:π, :sigma2=>:σ
	)

	growth_rate_select::Array{Bool, 1} = repeat([false], length(observed_vars_input))
	me_treatment::Symbol = :unbounded
	me_std_cutoff::Float64 = 0.2

	# meas_error_input::Array{Symbol, 1} = [ :σ, :Tgrowth, :τobs1, :w90share, :I90share]
	# meas_error_distr::Array{InverseGamma{Float64}, 1} = [InverseGamma(ig_pars(0.05, 0.1^2)...),InverseGamma(ig_pars(0.0005, 0.001^2)...),
	# InverseGamma(ig_pars(0.0005, 0.001^2)...),InverseGamma(ig_pars(0.0005, 0.001^2)...), InverseGamma(ig_pars(0.0005, 0.001^2)...)]
	#
	# meas_error_input::Array{Symbol, 1} = [ :Tgrowth]
	# meas_error_distr::Array{InverseGamma{Float64}, 1} = [InverseGamma(ig_pars(0.005, 0.001^2)...)]

	meas_error_input::Array{Symbol, 1} = [ ]
	meas_error_distr::Array{InverseGamma{Float64}, 1} = [ ]


	mode_start_file::String = "../saves/hank_2asset_mcmc_0705_baseline_G_chain_all.jld2"

	data_file::String = "../../Data/julia_ready/bbl_data_inequality_2018_FDMV.csv"
	save_mode_file::String = "../saves/hank_2asset_mcmc_1605_baseline_G.jld2"
	save_posterior_file::String = "../saves/hank_2asset_mcmc_1605_baseline_G.jld2"

	fd_flag::Bool = any(growth_rate_select)
	max_iter_mode::Int = 6000
	ndraws::Int      = 25000
	burnin::Int      = 10000
	mhscale::Float64 = 0.4
	debug_print::Bool = true
	mode_compute::Bool = true

end

macro make_struct_aggr(struct_name, a_names)
	# fields=[:($(entry.args[1])::$(entry.args[2])) for entry in var_names]
	# fieldsSS=[:($(Symbol((entry.args[1]), "SS"))::$(entry.args[2])) for entry in var_names]
	a_names=Symbol.(eval((a_names)))
	n_states = length(a_names)

	fields_states 	  = [:($(a_names[i])::Int) for i = 1:n_states]
	fieldsSS_states   = [:($(Symbol(a_names[i], "SS")) ::Int) for i = 1:n_states]
	# fields_controls   =[:($(control_names[i])::Int) for i = 1:n_controls]
	# fieldsSS_controls =[:($(Symbol(control_names[i], "SS")) ::Int) for i =1:n_controls]

	# fields=[:($(Symbol(i, "::Int"))) for i in var_names]
	# fieldsSS=[:($(Symbol(i, "SS::Int"))) for i in var_names]
	esc(quote struct $struct_name
	$(fieldsSS_states...)
	$(fields_states...)
end
end)
end

macro make_struct(struct_name, s_names, c_names)
	# fields=[:($(entry.args[1])::$(entry.args[2])) for entry in var_names]
	# fieldsSS=[:($(Symbol((entry.args[1]), "SS"))::$(entry.args[2])) for entry in var_names]
	state_names=Symbol.(eval((s_names)))
	n_states = length(state_names)
	control_names=Symbol.(eval((c_names)))
	n_controls = length(control_names)

	fields_states 	  = [:($(state_names[i])::Int) for i = 1:n_states]
	fieldsSS_states   = [:($(Symbol(state_names[i], "SS")) ::Int) for i = 1:n_states]
	fields_controls   =[:($(control_names[i])::Int) for i = 1:n_controls]
	fieldsSS_controls =[:($(Symbol(control_names[i], "SS")) ::Int) for i =1:n_controls]

	# fields=[:($(Symbol(i, "::Int"))) for i in var_names]
	# fieldsSS=[:($(Symbol(i, "SS::Int"))) for i in var_names]
	esc(quote struct $struct_name
	distr_m_SS	::Array{Int,1}
	distr_k_SS	::Array{Int,1}
	distr_y_SS	::Array{Int,1}
	DSS     	::Array{Int,1}
	$(fieldsSS_states...)
	VmSS     	::Array{Int,1}
	VkSS     	::Array{Int,1}
	$(fieldsSS_controls...)
	distr_m		::Array{Int,1}
	distr_k		::Array{Int,1}
	distr_y		::Array{Int,1}
	D     		::Array{Int,1}
	$(fields_states...)
	Vm     		::Array{Int,1}
	Vk     		::Array{Int,1}
	$(fields_controls...)
end
end)
end
# macro make_struct(struct_name, s_names, c_names)
#     # fields=[:($(entry.args[1])::$(entry.args[2])) for entry in var_names]
# 	# fieldsSS=[:($(Symbol((entry.args[1]), "SS"))::$(entry.args[2])) for entry in var_names]
# 	state_names=Symbol.(eval((s_names)))
# 	n_states = length(state_names)
# 	control_names=Symbol.(eval((c_names)))
# 	n_controls = length(control_names)
#
# 	fields_states 	  = [:($(state_names[i])::Int) for i = 1:n_states]
# 	fieldsSS_states   = [:($(Symbol(state_names[i], "SS")) ::Int) for i = 1:n_states]
# 	fields_controls   =[:($(control_names[i])::Int) for i = 1:n_controls]
# 	fieldsSS_controls =[:($(Symbol(control_names[i], "SS")) ::Int) for i =1:n_controls]
#
# 	# fields=[:($(Symbol(i, "::Int"))) for i in var_names]
# 	# fieldsSS=[:($(Symbol(i, "SS::Int"))) for i in var_names]
#      esc(quote struct $struct_name
# 	    distr_m_SS	::Array{Int,1}
# 	    distr_k_SS	::Array{Int,1}
# 	    distr_y_SS	::Array{Int,1}
# 	    VmSS     	::Array{Int,1}
# 	    VkSS     	::Array{Int,1}
# 		$(fieldsSS_states...)
# 		$(fieldsSS_controls...)
# 	    distr_m		::Array{Int,1}
# 	    distr_k		::Array{Int,1}
# 	    distr_y		::Array{Int,1}
# 	    Vm     		::Array{Int,1}
# 	    Vk     		::Array{Int,1}
# 		$(fields_states...)
#         $(fields_controls...)
#         end
#     end)
# end



# @with_kw struct IndexStruct
#     distr_m_SS::Array{Int,1}
#     distr_k_SS::Array{Int,1}
#     distr_y_SS::Array{Int,1}
#     ASS       ::Int
#     ZSS       ::Int
#     RBSS       ::Int
#     BSS       ::Int
#     GSS       ::Int
#     IlagSS ::Int
#     μSS       ::Int
#     YlagSS    ::Int
#     τSS       ::Int
#     σSS       ::Int
#     λSS        ::Int
#     GshockSS  ::Int
#     TshockSS  ::Int
#     RshockSS  ::Int
#     SshockSS  ::Int
#     LshockSS  ::Int
#     # Controls (SS)
#     VmSS     ::Array{Int,1}
#     VkSS     ::Array{Int,1}
#     rSS        ::Int
#     wSS       ::Int
#     KSS        ::Int
#     πSS        ::Int
#     YSS        ::Int
#     CSS        ::Int
#     qSS        ::Int
#     NSS        ::Int
#     mcSS       ::Int
#     uSS        ::Int
#     TSS        ::Int
#     ISS        ::Int
#     profitsSS  ::Int
#     #States (dev)
#     distr_m    ::Array{Int,1}
#     distr_k    ::Array{Int,1}
#     distr_y    ::Array{Int,1}
#     A          ::Int
#     Z          ::Int
#     RB          ::Int
#     G          ::Int
#     Ilag    ::Int
#     μ          ::Int
#     Ylag       ::Int
#     τ          ::Int
#     σ          ::Int
#     λ          ::Int
#     Gshock  ::Int
#     Tshock  ::Int
#     Rshock  ::Int
#     Sshock  ::Int
#     Lshock  ::Int
#
#     #Controls (dev)
#     Vm       ::Array{Int,1}
#     Vk       ::Array{Int,1}
#     r         ::Int
#     w         ::Int
#     K          ::Int
#     π          ::Int
#     Y          ::Int
#     C          ::Int
#     q          ::Int
#     N          ::Int
#     mc         ::Int
#     u         ::Int
#     T          ::Int
#     I          ::Int
#     B          ::Int
#     profits    ::Int
# end
