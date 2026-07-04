module TrapVoltagesHDF5Ext

using TrapVoltages.Solutions

using HDF5

function Solutions._center_tracker_hdf5(path, region)
    return h5open(path) do fh
        return Solutions.CenterTracker(read(fh[region], "yz_index"))
    end
end

end
