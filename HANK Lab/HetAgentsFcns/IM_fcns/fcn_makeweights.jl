function MakeWeights(xpol,grid)
    idx = zeros(Int64, size(xpol))
    weightright = zeros(typeof(xpol[1]), size(xpol))
    weightleft  = zeros(typeof(xpol[1]), size(xpol))
    dx = diff(grid)
    for i= 1:length(xpol)
        if xpol[i] .<= grid[1]
            idx[i] = 1
        elseif xpol[i] .>= grid[end]
            idx[i] = length(grid) .- 1
        else
            idx[i] = locate(xpol[i], grid)
        end
        weightright[i] = (xpol[i] - grid[idx[i]]) / dx[idx[i]]
        if weightright[i] >= 1.0
            weightright[i] = 1.0 - 1.0e-14
        elseif weightright[i] <= 0.0
            weightright[i] = 1.0e-14
        end
        weightleft[i] = 1.0 - weightright[i]
    end
    return idx, weightright, weightleft
end

function MakeWeightsLight(xpol,grid)
    idx = zeros(Int64, size(xpol))
    weightright = zeros(eltype(xpol), size(xpol))
    #weightleft  = zeros(typeof(xpol[1]), size(xpol))
    dx = diff(grid)
    @inbounds begin
        for i in eachindex(xpol)
            if xpol[i] .<= grid[1]
                idx[i] = 1
            elseif xpol[i] .>= grid[end]
                idx[i] = length(grid) .- 1
            else
                idx[i] = locate(xpol[i], grid)
            end
            weightright[i] = (xpol[i] .- grid[idx[i]]) ./ dx[idx[i]]
            if weightright[i] >= 1.0
                weightright[i] = 1.0 - 1.0e-14
            elseif weightright[i] <= 0.0
                weightright[i] = 1.0e-14
            end
        #    weightleft[i] = 1.0 - weightright[i]
        end
    end
    return idx, weightright#, weightleft
end
# function MakeWeightsScalar(xpol,grid)
#     idx = 0
#     weightright = zeros(eltype(xpol))
#     #weightleft  = zeros(typeof(xpol[1]), size(xpol))
#     dx = diff(grid)
#     #for i in eachindex(xpol)
#         if xpol .<= grid[1]
#             idx = 1
#         elseif xpol .>= grid[end]
#             idx = length(grid) .- 1
#         else
#             idx = locate(xpol, grid)
#         end
#         weightright = (xpol - grid[idx]) / dx[idx]
#         if weightright >= 1.0
#             weightright = 1.0 - 1.0e-14
#         elseif weightright <= 0.0
#             weightright = 1.0e-14
#         end
#     #    weightleft[i] = 1.0 - weightright[i]
#     # end
#     return idx, weightright#, weightleft
# end
