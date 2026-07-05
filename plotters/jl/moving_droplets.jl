using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")

##

#plot_iso_bce(psi)


##

psi = load_from_text(joinpath("../../build", "initial_state.gpe.dat"))

function move_droplet(psi::IsoBECContext)
    slice = get_BEC_slice(psi)
    maxs = find_local_maxima(slice)

    if(length(maxs) < 2)
        error("Can't find enough droplets")
    end

    middle_idx = floor(Int, maxs[1].i + (maxs[2].i - maxs[1].i) / 2)
    move_idxs = 10
    psi_moved = Array{Float64, 3}(undef, psi.nx, psi.ny, psi.nz)

    for i in 1:psi.nx
        for j in 1:psi.ny
            for k in 1:psi.nz
                if(i - move_idxs > 0 && i < middle_idx)
                    psi.psi[i - move_idxs,j,k] = psi.psi[i,j,k]
                    psi.psi[i,j,k] = 0
                end
            end
        end
    end
    
    return psi
end


psi_2 = move_droplet(psi)

##

plot_iso_bce(psi_2)
save_to_text(psi_2, "../../build/moved.gpe.dat")

##

psi_test = load_from_text("../../build/moved.gpe.dat")

plot_iso_bce(psi_test)
