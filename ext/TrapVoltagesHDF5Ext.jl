module TrapVoltagesHDF5Ext

using TrapVoltages.Solutions

using HDF5

function _get_path(name)
    if !('/' in name)
        p = joinpath(@__DIR__, "../data", name, "rf_center.h5")
        if isfile(p)
            return p
        end
    end
    return name
end

function Solutions.CenterTracker(name::AbstractString, region=1)
    return h5open(_get_path(name)) do fh
        return Solutions.CenterTracker(read(fh["$region"], "zy_index"))
    end
end

end
