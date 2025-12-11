#using CairoMakie
using GLMakie

include("context.jl")

function plot_iso_bce(BECContext::IsoBECContext)
    fig = Figure()
    ax = Axis3(fig[1, 1], xlabel="X", ylabel="Y", zlabel="Z", title="Isosurface of BEC Structure", aspect = :data)

    rho = abs.(BECContext.psi)
    
    x = BECContext.x
    y = BECContext.y
    z = BECContext.z

    max_val = maximum(rho)
    alpha = 0.8
    colormap = :plasma
    transparency = true
    algorithm = :iso
    isovals = [0.2 * max_val, 0.5 * max_val, 0.8 * max_val]
    alphas = [0.3, 0.5, 0.9]

    for (i, isovalue) in enumerate(isovals)
        volume!(ax,
                (x[begin], x[end]),
                (y[begin], y[end]),
                (z[begin], z[end]),
                rho,
                algorithm = :iso,
                isovalue = isovalue,
                alpha = alphas[i],
                colormap = :YlGn,
                transparency = true,
                isorange = 0.1 * max_val)
    end

    fig
end

function plot_heatmap_cut(data)
    fig = Figure()
    ax  = Axis(fig[1, 1], xlabel="X", ylabel="Y", title="Density Cut at Z=0")

    heatmap!(ax, data, colormap = :viridis)

    fig
end
