#

using Test
using TrapVoltages

@testset "TrapDesc" begin
    px = TrapDesc("phoenix")
    @test px.name == "phoenix"
    @test length(px.ele_names) == length(px.ele_indices) == 92
    pr = TrapDesc("peregrine")
    @test pr.name == "peregrine"
    @test length(pr.ele_names) == length(pr.ele_indices) == 92
    hoa = TrapDesc("hoa")
    @test hoa.name == "hoa"
    @test length(hoa.ele_names) == length(hoa.ele_indices) == 96

    @test TrapDesc(px) === px
    @test TrapDesc(pr) === pr
    @test TrapDesc(hoa) === hoa

    @test_throws ArgumentError TrapDesc("unknown")

    # Make coverage report happy...
    pabc = TrapDesc("abc", copy(px.ele_names))
    @test px.ele_indices == pabc.ele_indices
end
