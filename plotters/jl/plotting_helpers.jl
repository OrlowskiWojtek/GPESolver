using CairoMakie
using GLMakie

if !isdefined(@__MODULE__, :_bec_density_colormap)
    function _bec_density_colormap()
        positions = [0.0, 3/24, 3/12, 6/12, 10/12, 1.0]
        colors = [
            CairoMakie.RGB(0.0, 0.0, 1.0),
            CairoMakie.RGB(1.0, 1.0, 1.0),
            CairoMakie.RGB(0.0, 1.0, 0.0),
            CairoMakie.RGB(1.0, 1.0, 0.0),
            CairoMakie.RGB(1.0, 0.0, 0.0),
            CairoMakie.RGB(0.0, 0.0, 0.0),
        ]

        return CairoMakie.cgrad(colors, positions)
    end
end
