#!/usr/bin/julia

using Distributed

env_nproc = get(ENV, "TRAP_VOLTAGES_NPROC", nothing)
if env_nproc === nothing
    addprocs(Sys.CPU_THREADS)
else
    env_nproc = parse(Int, env_nproc)
    if env_nproc > 0
        addprocs(env_nproc)
    end
end

pmap(["optimizers",
      "poly_fit",
      "potentials",
      "trap",
      "units"]) do file
    println("Start testing $file")
    @eval module $(Symbol("Test_$(file)_mod"))
    include($(joinpath(@__DIR__, "$(file).jl")))
    end
    println("Done testing $file")
    # So that we do not try to bring the worker-only module back to the main process
    return
end
