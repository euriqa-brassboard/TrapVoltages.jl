#!/usr/bin/julia

using TrapVoltages: Potentials, Solutions

using HDF5

# This should be pointing to the RS1394_Peregrine_coarse.va file (without extension).
const potential_prefix = ARGS[1]
const data_prefix = joinpath(@__DIR__, "../data/peregrine/rf_center")

function compute_center(file, trap)
    potential = Potentials.import_pillbox_64(file, trap=trap)
    # RF is electrode 2 (ground is 1)
    centers = Solutions.find_all_flat_points(@view(potential.data[:, :, :, 2]))
    centers_um = similar(centers)
    for i in 1:potential.nx
        centers_um[i, 1] = Potentials.z_index_to_axis(potential, centers[i, 1]) * 1000
        centers_um[i, 2] = Potentials.y_index_to_axis(potential, centers[i, 2]) * 1000
    end
    return centers, centers_um
end

mkpath(dirname(data_prefix))

h5open("$(data_prefix).h5", "w") do fh
    g = create_group(fh, "1") # region 1
    centers, centers_um = compute_center("$(potential_prefix).va", "peregrine")
    g["zy_index"] = centers
    g["zy_um"] = centers_um
end
