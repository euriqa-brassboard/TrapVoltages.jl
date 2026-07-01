#

using TrapVoltages: Solutions

using Test

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
