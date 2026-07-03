#!/usr/bin/julia

using HDF5

using TrapVoltages: Potentials, Solutions, PolyFit

# This should be pointing to the RS1394_coarse.va file (without extension).
const potential_prefix = ARGS[1]
const data_prefix = joinpath(@__DIR__, "../data/phoenix/rf_trap")

function compute_trap(file, trap, region=1)
    centers = Solutions.CenterTracker(trap, region)
    potential = Potentials.import_pillbox_64(file, trap=trap)
    # RF is electrode 2 (ground is 1)
    rf_data = potential.data[:, :, :, 2]
    yz_fitter = PolyFit.Fitter(4, 4, sizes=(5, 5))
    rf_y2s = Vector{Float64}(undef, potential.nx)
    rf_yzs = Vector{Float64}(undef, potential.nx)
    rf_z2s = Vector{Float64}(undef, potential.nx)
    ystride_m = potential.stride[2] / 1000
    zstride_m = potential.stride[3] / 1000
    for xidx in 1:potential.nx
        yidx, zidx = get(centers, xidx)
        fit_cache = PolyFit.FitCache(yz_fitter, @view rf_data[xidx, :, :])
        fit = get(fit_cache, (yidx, zidx))
        y2 = fit[2, 0] / ystride_m^2
        yz = fit[1, 1] / zstride_m / ystride_m / 2
        z2 = fit[0, 2] / zstride_m^2
        rf_y2s[xidx] = y2 # 1/m^2
        rf_yzs[xidx] = yz # 1/m^2
        rf_z2s[xidx] = z2 # 1/m^2
    end
    xs_um = Potentials.x_index_to_axis.(Ref(potential), 1:potential.nx) .* 1000
    return xs_um, rf_y2s, rf_yzs, rf_z2s
end

mkpath(dirname(data_prefix))

h5open("$(data_prefix).h5", "w") do fh
    g = create_group(fh, "1") # region 1
    xs_um, rf_y2s, rf_yzs, rf_z2s = compute_trap("$(potential_prefix).va", "phoenix")
    g["xs_um"] = xs_um
    fields = create_group(g, "field")
    fields["y2"] = rf_y2s
    fields["yz"] = rf_yzs
    fields["z2"] = rf_z2s
end
