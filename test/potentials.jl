#

using TrapVoltages.Potentials

using Test

const data_dir = joinpath(@__DIR__, "test_data")

@testset "Raw Potentials" begin
    p1_v0 = Potentials.import_pillbox_v0_raw(joinpath(data_dir, "dummy_v0_e92.bin"))
    @test p1_v0.electrodes == 92
    @test p1_v0.nx == 5
    @test p1_v0.ny == 4
    @test p1_v0.nz == 3
    @test all(p1_v0.stride .≈ (1e-3, 1e-3, 1e-3))
    @test all(p1_v0.origin .≈ (-3e-3, -1e-3, -2e-3))
    @test size(p1_v0.data) == (3, 4, 5, 92)

    p2_v0 = Potentials.import_pillbox_v0_raw(joinpath(data_dir, "dummy_v0_e96.bin"))
    @test p2_v0.electrodes == 96
    @test p2_v0.nx == 2
    @test p2_v0.ny == 4
    @test p2_v0.nz == 4
    @test all(p2_v0.stride .≈ (1e-3, 1e-3, 1e-3))
    @test all(p2_v0.origin .≈ (-2e-3, -2e-3, 1e-3))
    @test size(p2_v0.data) == (4, 4, 2, 96)

    @test_throws ArgumentError Potentials.import_pillbox_v0_raw(joinpath(data_dir, "dummy_v0_extra.bin"))
    @test_throws ArgumentError Potentials.import_pillbox_v0_raw(joinpath(data_dir, "dummy_v0_short.bin"))
end
