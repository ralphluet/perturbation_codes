function Tauchen(rho::Float64, N::Int; sigma::Float64 = 1.0, mue::Float64 = 0.0) #, transtype::Symbol = :importance)
#TAUCHEN Generates a discrete approximation to an AR-1 Process, following
#   Tauchen (1987).
#   [grid(),P] = Tauchen(rho,N,sigma, mue, type())
#
#   returns a state vector GRID and a transition probability matrix P.
#
#   Input arguments:
#   rho:    autocorrelation coefficient
#   N:      number of gridpoints
#   sigma:  long-run variance
#   mue:    mean of the ar-1 process
#   type(): {importance | equi | s>imple} - type of grid()-Transition generating
#   algorithm.
#
#   importance: importance sampling. Each bin has probability 1/N to
#               realize
#   equi:       bin-centers are equi-spaced between +-3 std()
#   simple:     like equi + Transition Probabilities are calculated without
#               using integrals.
#   simple_importance: like simple but with grid from importance
#

#   Author: Christian Bayer, Uni Bonn, 03.05.2010

# Note: Remove "Ref()"'s when https://github.com/JuliaStats/Distributions.jl/pull/760 is
#       included in Distributions.jl tagged version

dis = Normal()
pr_ij(x, bound1, bound2, rho, sigma_e) = pdf.(Ref(dis), x) .*
     (cdf.(Ref(dis), (bound2 - rho .* x) ./ sigma_e) -
     cdf.(Ref(dis), (bound1 - rho .* x) ./ sigma_e) )

# if transtype == :importance # Importance Sampling
    grid_probs = range(0.0, stop = 1.0, length = N+1)   # generate equi-likely bins

    bounds = quantile.(Ref(dis), grid_probs[1:end])     # corresponding bin bounds

    # replace [-]Inf bounds, by finite numbers
    bounds[1] = bounds[2] - 1.0e2
    bounds[end] = bounds[end-1] + 1.0e2

    # Calculate grid() - centers
    grid_vec = N * (pdf.(Ref(dis), bounds[1:end-1]) - pdf.(Ref(dis), bounds[2:end]))

    sigma_e = sqrt(1 - rho^2)# Calculate short run variance
    P = fill(0.0, (N, N)) # Initialize Transition Probability Matrix

    for j = 1:N
        p(x) = pr_ij(x,bounds[j], bounds[j+1], rho, sigma_e)
        for i = 1:floor(Int, (N-1)/2)+1 # Exploit Symmetrie to save running time
            P[i, j] = my_integrate(p, bounds[i], bounds[i+1]) # Evaluate Integral
        end
    end
    # Exploit Symmetrie Part II
    P[floor(Int, (N - 1) / 2) + 2:N, :] = P[(ceil(Int, (N - 1) / 2):-1:1), end:-1:1]

# elseif transtype == :equi # use +-3 std equi-spaced grid()
#
#     ## Equi-spaced
#     stepsize = 2.0 / (N - 1)
#     grid_vec = range(-1.0, stop = 1.0, step = stepsize)
#     # display(typeof(grid))
#     bounds = vcat(-1.0e2, grid_vec[1:end-1] .+ stepsize/2.0, 1.0e2) # Bounds of Bins
#     # display(typeof(bounds))
#     sigma_e = sqrt(1-rho^2)# Calculate short run variance
#     P = fill(0.0, (N, N)) # Initialize Transition Probability Matrix
#
#     for j = 1:N
#         p(x) = pr_ij(x,bounds[j],bounds[j+1], rho, sigma_e)
#         for i = 1:floor(Int,(N-1)/2)+1 # Exploit Symmetrie to save running time
#             P[i, j] = my_integrate(p, bounds[i], bounds[i+1]) # Evaluate Integral
#         end
#     end
#     # Exploit Symmetrie Part II
#     P[floor(Int, (N-1)/2) + 2:N, :] = P[(ceil(Int, (N-1)/2):-1:1),end:-1:1]
#
# elseif transtype == :simple # use simple Transition Probabilities
#         # Generate Grid
#         stepsize = 6.0 / (N - 1)
#         grid_vec = range(-3, stop = 3, step = stepsize)
#         bounds = Array{Float64}
#         sigma_e = sqrt(1-rho^2)# Calculate short run STD
#         P = fill(0.0, (N, N)) # Initialize Transition Probability Matrix
#         for i = 1:N
#             P[i, 1] = cdf.(Ref(dis), (grid_vec[1] + step / 2 - rho * grid_vec[i]) / sigma_e)
#             P[i, N] = 1 - cdf.(Ref(dis), (grid_vec[N] - step / 2 - rho * grid_vec[i]) / sigma_e)
#             for j = 2:N-1
#                 P[i, j] = cdf.(Ref(dis), (grid_vec[j] + step / 2 - rho * grid_vec[i]) / sigma_e) -
#                        cdf.(Ref(dis), (grid_vec[j] - step / 2 - rho * grid_vec[i]) / sigma_e)
#             end
#         end
#
# elseif transtype == :simple_importance # use simple Transition Probabilities
#         # Generate Grid
#         grid_probs = range(0, stop = 1, length = N+1)   # generate equi-likely bins
#         bounds     = quantile.(Ref(dis),grid_probs)     # corresponding bin bounds
#
#         # Calculate grid() - centers
#         grid_vec = N * ( pdf.(Ref(dis), bounds[1:end-1]) - pdf.(Ref(dis), bounds[2:end]) )
#
#         # replace [-]Inf bounds, by finite numbers
#         bounds[1] = bounds[2] - 1e2
#         bounds[end] = bounds[end - 1] + 1.0e2
#
#         sigma_e = sqrt(1 - rho^2) # Calculate short run variance
#         P = fill(0.0, (N, N)) # Initialize Transition Probability Matrix
#         for i = 1:floor(Int,(N-1)/2+1) # Exploit Symmetrie to save running time
#             P[i, 1] = cdf.(Ref(dis), (bounds[2] - rho * grid_vec[i]) / sigma_e)
#             P[i, N] = 1 - cdf.(Ref(dis), (bounds[end-1] - rho * grid_vec[i]) / sigma_e)
#             for j = 2:N-1
#                 P[i, j]= cdf.(Ref(dis), (bounds[j+1] - rho * grid_vec[i]) / sigma_e) -
#                         cdf.(Ref(dis), (bounds[j] - rho * grid_vec[i]) / sigma_e)
#             end
#         end
#         # Exploit Symmetrie Part II
#         P[floor(Int, (N-1)/2+2):end,:] = P[ceil(Int,(N-1)/2):-1:1,end:-1:1]
# else
#         error("Transition type not correctly defined")
# end

# Make sure P is a Probability Matrix
P = P ./ sum(P, dims = 2)

grid_vec   = grid_vec .* sigma .+ mue
lmul!(sigma, bounds)
# bounds = bounds .* sigma

return grid_vec, P, bounds

end

function ExTransition(rho::Number,bounds::Array{Float64,1},riskscale::Number)
#similar to TAUCHEN
N = length(bounds)-1
# Assume Importance Sampling
        sigma_e=riskscale*sqrt(1-rho^2)# Calculate short run variance
        P=zeros(typeof(riskscale),N,N) # Initialize Transition Probability Matrix

        for i=1:floor(Int,(N-1)/2)+1
            nodes, weights = my_qnwcheb(500, bounds[i], bounds[i+1])
            for j=1:N
            p(x) = pr_ij(x,bounds[j],bounds[j+1],rho,sigma_e)
             # Exploit Symmetrie to save running time
                P[i,j]=dot(weights, p.(nodes))# my_integrate(p,bounds[i],bounds[i+1]) # Evaluate Integral
            end
        end

       # Exploit Symmetrie Part II
        P[floor(Int,(N-1)/2)+2:N,:]=P[(ceil(Int,(N-1)/2):-1:1),end:-1:1]

# Make sure P is a Probability Matrix
P = P./ sum(P, dims = 2)

return P

end

# function pr_ij(x,bound1,bound2,rho,sigma_e)
#     dis = Normal()
#     p = pdf.(dis,x) .* (cdf.(dis,(bound2 - rho.*x)./sigma_e) - cdf.(dis,(bound1 - rho.*x)./sigma_e) )
# return p
# end
function pr_ij(x,bound1,bound2,rho,sigma_e)
    mycdf(x) = 0.5 + 0.5 * erf.(x / sqrt(2.0))
    mypdf(x) =  1/sqrt(2*π).*exp.(-x.^2/2.0)
    p = mypdf.(x) .* (mycdf.((bound2 - rho.*x)./sigma_e) - mycdf.((bound1 - rho.*x)./sigma_e) )
return p
end
