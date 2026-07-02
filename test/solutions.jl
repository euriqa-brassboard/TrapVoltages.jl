#

using TrapVoltages: Solutions, CenterTracker, PolyFit as PF
using TrapVoltages.Units: TrapUnits

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
    @test CenterTracker("hoa", 2).zy_index == CenterTracker(joinpath(@__DIR__, "../data/hoa/rf_center.h5"), 2).zy_index
    CenterTracker("hoa", 3)
    CenterTracker("hoa", 4)
    CenterTracker("hoa", 5)
end

macro test_term(term)
    term_v = Symbol("$(term)_v")
    quote
        if $(esc(term))
            @test $(esc(:(terms.$term))) ≈ $(esc(term_v))
        else
            @test !hasproperty($(esc(:terms)), $(QuoteNode(term)))
        end
    end
end

@testset "Compensation Terms" begin
    fit = PF.Result{3}((2, 2, 4), zeros(3 * 3 * 5))
    for (dx, dy, dz, xy, yz, zx, z2,
         x2, x3, x4, x2z) in Iterators.product((false, true), (false, true),
                                               (false, true), (false, true),
                                               (false, true), (false, true),
                                               (false, true), (false, true),
                                               (false, true), (false, true),
                                               (false, true))
        mask = Solutions.TermMask(dx=dx, dy=dy, dz=dz, xy=xy, yz=yz, zx=zx, z2=z2,
                                  x2=x2, x3=x3, x4=x4, x2z=x2z)

        nterms = dx + dy + dz + xy + yz + zx + z2 + x2 + x3 + x4 + x2z
        if 3 < nterms < 9
            # Limit compilation
            continue
        end

        for _ in 1:10
            unit = TrapUnits(rand() + 1, rand() + 1)
            stride = (rand() + 0.5, rand() + 0.5, rand() + 0.5)

            fit.coefficient .= rand.()
            dx_v = rand()
            dy_v = rand()
            dz_v = rand()
            xy_v = rand()
            yz_v = rand()
            zx_v = rand()
            z2_v = rand()
            x2_v = rand()
            x3_v = rand()
            x4_v = rand()
            x2z_v = rand()

            x2y2z2_v = rand()
            real_x2 = x2y2z2_v + x2_v / 2
            real_y2 = x2y2z2_v - x2_v / 4 - z2_v / 2
            real_z2 = x2y2z2_v - x2_v / 4 + z2_v / 2

            fit[0, 0, 1] = dx_v * stride[1] / 1e6
            fit[0, 1, 0] = dy_v * stride[2] / 1e6
            fit[1, 0, 0] = dz_v * stride[3] / 1e6

            fit[0, 1, 1] = xy_v * stride[1] * stride[2] * unit.V_unit / unit.l_unit_um^2
            fit[1, 1, 0] = yz_v * stride[2] * stride[3] * unit.V_unit / unit.l_unit_um^2
            fit[1, 0, 1] = zx_v * stride[3] * stride[1] * unit.V_unit / unit.l_unit_um^2

            fit[0, 0, 2] = real_x2 * stride[1]^2 * unit.V_unit / unit.l_unit_um^2
            fit[0, 2, 0] = real_y2 * stride[2]^2 * unit.V_unit / unit.l_unit_um^2
            fit[2, 0, 0] = real_z2 * stride[3]^2 * unit.V_unit / unit.l_unit_um^2

            fit[0, 0, 3] = x3_v / 6 * stride[1]^3 * unit.V_unit / unit.l_unit_um^3
            fit[0, 0, 4] = x4_v / 24 * stride[1]^4 * unit.V_unit / unit.l_unit_um^4

            fit[1, 0, 2] = x2z_v / 2 * stride[1]^2 * stride[3] * unit.V_unit / unit.l_unit_um^3

            terms = Solutions.compensate_terms(fit, stride, unit=unit, mask=mask)

            @test_term dx
            @test_term dy
            @test_term dz

            @test_term xy
            @test_term yz
            @test_term zx

            @test_term z2
            @test_term x2
            @test_term x3
            @test_term x4

            @test_term x2z
        end
    end
end
