using CairoMakie
using GLMakie

## Finding correct potential profile for ball hitting (only x dimension right now)

##

x  = LinRange(-20000 / 2, 20000 /2 , 1000) * 18.897261260649092 # nm to au

dd    = 6000 * 18.897261260649092
dd_um = 6000 / 1000

wzl = 120 * 4.1356e-12 / 27211.6; # angular frequency of harmonic potential - z direction
wrl = 60. * 4.1356e-12 / 27211.6; # angular frequency of harmonic potential - y direction

m   = 163.929 / 0.000548579909;
aa  = m * wrl^2 / (4. * dd^2) # aa depends on dd / 2 -> kinda weird, but ok
vx  = aa * x .^ 4

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

##

n = 3
x0 = -dd
step = 2 * dd / n
V0 = b * x0^2
is_left = [_x < x0 + step for _x in x]
k = [left ? π * n / (2 * x0) : π * n / x0 for left in is_left]

step_offset = [left ? step : 0 for left in is_left]
periodic = [V0 * cos(_k * (_x - x0 - _stp_off))  for (_k, _x, _stp_off) in zip(k,x, step_offset)]

vx = @. aa * x^4

fig = Figure();
ax  = Axis(fig[1,1]);

lines!(ax, x_um, vx, color = :blue)
lines!(ax, x_um, periodic, color = :red)
lines!(ax, x_um, vx .+ periodic, color = :orange)
vlines!(ax, [-dd_um, dd_um], color = :black, linestyle = :dash)

fig

##

n = 3
x0 = -dd
step = 2 * dd / n
V0 = aa * x0^4
vx = @. aa * x^4

sigma = step / 3
movement = @. V0 * exp(-(x - x0)^2 / sigma^2)

fig = Figure();
ax  = Axis(fig[1,1]);

lines!(ax, x_um, vx, color = :blue)
lines!(ax, x_um, vx - movement, color = :red)
vlines!(ax, [-dd_um, dd_um], color = :black, linestyle = :dash)

fig
