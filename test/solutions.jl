#

using TrapVoltages: Solutions, CenterTracker

using Test
using HDF5

@testset "find_flat_point" begin
    xs = 1:10
    ys = 1:10
    function test_find_flat(x0, y0, scale)
        data = [(x - x0)^2 + (y - y0)^2 * scale
                for x in xs, y in ys]
        @test collect(Solutions.find_flat_point(data)) ≈ [x0, y0] atol=1e-8
    end
    for x0 in range(0, 8, 10)
        for y0 in range(0, 8, 10)
            for scale in (-1.0, -0.1, 0.1, 1)
                test_find_flat(x0, y0, scale)
            end
        end
    end

    zs = 1:12
    function x0_z(z)
        z = (z - 500) / 500
        return z * 0.5
    end
    function y0_z(z)
        z = (z - 500) / 500
        return z^3 * 0.3
    end
    for scale in (-1.0, -0.1, 0.1, 1)
        data = [(x - x0_z(z))^2 + (y - y0_z(z))^2 * scale
                for x in xs, y in ys, z in zs]
        all_flats = Solutions.find_all_flat_points(data)
        for (zi, zv) in enumerate(zs)
            @test all_flats[zi, :] ≈ [x0_z(zv), y0_z(zv)] atol=1e-4
        end
    end
end

@testset "CenterTracker" begin
    cpx = CenterTracker("phoenix")
    cpr = CenterTracker("peregrine")

    yz_px1 = get(cpx, 3220)
    yz_px2 = get(cpx, 3221)
    yz_px3 = get(cpx, 3220.5)
    @test yz_px1[1] ≈ 4.1 atol=0.1
    @test yz_px1[2] ≈ 3.5 atol=0.1
    @test yz_px2[1] ≈ 4.1 atol=0.1
    @test yz_px2[2] ≈ 3.5 atol=0.1
    @test yz_px3[1] ≈ (yz_px1[1] + yz_px2[1]) / 2
    @test yz_px3[2] ≈ (yz_px1[2] + yz_px2[2]) / 2

    yz_pr1 = get(cpr, 3220)
    yz_pr2 = get(cpr, 3221)
    yz_pr3 = get(cpr, 3220.5)
    @test yz_pr1[1] ≈ 3.0 atol=0.1
    @test yz_pr1[2] ≈ 11.3 atol=0.1
    @test yz_pr2[1] ≈ 3.0 atol=0.1
    @test yz_pr2[2] ≈ 11.3 atol=0.1
    @test yz_pr3[1] ≈ (yz_pr1[1] + yz_pr2[1]) / 2
    @test yz_pr3[2] ≈ (yz_pr1[2] + yz_pr2[2]) / 2

    @test CenterTracker("hoa").zy_index == CenterTracker("hoa", 1).zy_index
    CenterTracker("hoa", 2)
    CenterTracker("hoa", 3)
    CenterTracker("hoa", 4)
    CenterTracker("hoa", 5)
end
