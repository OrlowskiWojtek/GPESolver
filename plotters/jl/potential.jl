using CairoMakie
using GLMakie

## Finding correct potential profile for ball hitting (only x dimension right now)

##

x  = LinRange(-14520 / 2, 14520 /2 , 1000) * 18.897261260649092 # nm to au
Cdd = 0.0165267

dd = 1500 * 18.897261260649092
wzl = 120 * 4.1356e-12 / 27211.6; # angular frequency of harmonic potential - z direction
wrl = 60. * 4.1356e-12 / 27211.6; # angular frequency of harmonic potential - y direction

m   = 163.929 / 0.000548579909;
aa  = m * wrl^2 / (4. * dd^2) # aa depends on dd / 2 -> kinda weird, but ok
vx  = aa * x .^ 4

##

b  = 0.5 * m * wrl^2;
vx_barrier  = @. -b * (x) ^ 2 + aa * x ^ 4
vx_normal   = @. aa * (x) ^ 4

vx_offset   = Float64[]
for _x in x
    val = 0.

    if(_x < 0)
        val = (aa * _x^4)
    end

    push!(vx_offset, val)
end
##

fig  = Figure();
ax   = Axis(fig[1,1]);

x_um = x ./ 18.897261260649092 ./ 1000
 
lines!(ax, x_um, vx_barrier, label = "with barrier", color = :blue)
lines!(ax, x_um, vx_normal, label = "normal", color = :black)
lines!(ax, x_um, vx_offset, label = "huśtawka", color = :orange, linestyle = :dash)

xlims!(ax, (-5, 5))
ylims!(ax, (-1e-14, 1e-13))

axislegend()

fig
