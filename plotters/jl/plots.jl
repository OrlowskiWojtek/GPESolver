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

function plot_energies(context::EnergiesContext)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Iteration", ylabel="Energy (meV)")

    # Energy conversion factor from Hartree to meV
    conversion_factor = 27.2114 * 1000

    lines!(ax, context.iter, context.e_kin * conversion_factor, label="Kinetic Energy")
    lines!(ax, context.iter, context.e_pot * conversion_factor, label="Potential Energy")
    lines!(ax, context.iter, context.e_int * conversion_factor, label="Internal Energy")
    lines!(ax, context.iter, context.e_ext * conversion_factor, label="External Energy")
    lines!(ax, context.iter, context.e_bmf * conversion_factor, label="BMF Energy")
    lines!(ax, context.iter, context.e_tot * conversion_factor, label="Total Energy")

    axislegend(ax)
    fig
end

function plot_slice(slice::IsoBECSlice)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "x [nm]", ylabel = "y [nm]")

    rho = abs.(slice.psi)
    heatmap!(ax, slice.x, slice.y, rho, colormap = :plasma)

    return fig;
end


