#

using TrapVoltages: Solutions, CenterTracker, PolyFit as PF,
TrapDesc, Potentials, Traps, TrapUnits

using Test
using LinearAlgebra

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
                for z in zs, y in ys, x in xs]
        all_flats = Solutions.find_all_flat_points(data)
        for (zi, zv) in enumerate(zs)
            @test all_flats[zi, :] ≈ [y0_z(zv), x0_z(zv)] atol=1e-4
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

    @test CenterTracker("hoa").yz_index == CenterTracker("hoa", 1).yz_index
    @test CenterTracker("hoa", 2).yz_index == CenterTracker(joinpath(@__DIR__, "../data/hoa/rf_center.h5"), 2).yz_index
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

unwrap(::Val{V}) where V = V

@testset "Compensation Terms" begin
    fit = PF.Result{3}((4, 2, 2), zeros(5 * 3 * 3))
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

        @test count(unwrap(mask)) == nterms
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

            fit[1, 0, 0] = dx_v * stride[1] / 1e6
            fit[0, 1, 0] = dy_v * stride[2] / 1e6
            fit[0, 0, 1] = dz_v * stride[3] / 1e6

            fit[1, 1, 0] = xy_v * stride[1] * stride[2] * unit.V_unit / unit.l_unit_um^2
            fit[0, 1, 1] = yz_v * stride[2] * stride[3] * unit.V_unit / unit.l_unit_um^2
            fit[1, 0, 1] = zx_v * stride[3] * stride[1] * unit.V_unit / unit.l_unit_um^2

            fit[2, 0, 0] = real_x2 * stride[1]^2 * unit.V_unit / unit.l_unit_um^2
            fit[0, 2, 0] = real_y2 * stride[2]^2 * unit.V_unit / unit.l_unit_um^2
            fit[0, 0, 2] = real_z2 * stride[3]^2 * unit.V_unit / unit.l_unit_um^2

            fit[3, 0, 0] = x3_v / 6 * stride[1]^3 * unit.V_unit / unit.l_unit_um^3
            fit[4, 0, 0] = x4_v / 24 * stride[1]^4 * unit.V_unit / unit.l_unit_um^4

            fit[2, 0, 1] = x2z_v / 2 * stride[1]^2 * stride[3] * unit.V_unit / unit.l_unit_um^3

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

@testset "Compensation Terms Fit" begin
    trap = TrapDesc("phoenix")
    origin = (-500, -4, -4)
    p = Potentials.Potential(length(trap.ele_names), 1001, 9, 9, (1, 1, 1), origin,
                             PermutedDimsArray(zeros(9, 9, 1001, length(trap.ele_names)),
                                               (3, 2, 1, 4)),
                             trap.ele_indices, [[n] for n in trap.ele_names], trap)
    unit = TrapUnits((rand() + 0.5) * 1e-6, rand() + 0.5)
    eles = sort!(collect(Traps.find_electrodes(trap, trap.ele_indices, 0, min_num=20)))
    neles = length(eles)
    @test neles >= 20

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

        @test count(unwrap(mask)) == nterms
        if 2 < nterms < 10 || nterms == 0
            # Limit compilation
            continue
        end

        function fill_potential!(vs, eid)
            dx_c, dy_c, dz_c, xy_c, yz_c, zx_c, z2_c, x2_c, x3_c, x4_c, x2z_c = vs
            x2y2z2_c = rand()
            real_x2 = x2y2z2_c + x2_c / 2
            real_y2 = x2y2z2_c - x2_c / 4 - z2_c / 2
            real_z2 = x2y2z2_c - x2_c / 4 + z2_c / 2

            x1_c = dx_c / 1e6
            y1_c = dy_c / 1e6
            z1_c = dz_c / 1e6

            xy_c = xy_c * unit.V_unit / unit.l_unit_um^2
            yz_c = yz_c * unit.V_unit / unit.l_unit_um^2
            zx_c = zx_c * unit.V_unit / unit.l_unit_um^2

            x2_c = real_x2 * unit.V_unit / unit.l_unit_um^2
            y2_c = real_y2 * unit.V_unit / unit.l_unit_um^2
            z2_c = real_z2 * unit.V_unit / unit.l_unit_um^2

            x3_c = x3_c / 6 * unit.V_unit / unit.l_unit_um^3
            x4_c = x4_c / 24 * unit.V_unit / unit.l_unit_um^4
            x2z_c = x2z_c / 2 * unit.V_unit / unit.l_unit_um^3

            for x in 1:p.nx
                xv = Potentials.x_index_to_axis(p, x)
                x2v = xv^2
                x3v = xv^3
                x4v = xv^4
                for y in 1:p.ny
                    yv = Potentials.y_index_to_axis(p, y)
                    y2v = yv^2
                    xyv = xv * yv
                    for z in 1:p.nz
                        zv = Potentials.z_index_to_axis(p, z)
                        z2v = zv^2
                        yzv = yv * zv
                        zxv = zv * xv
                        x2zv = x2v * zv
                        p.data[x, y, z, eid] = (x1_c * xv + y1_c * yv + z1_c * zv +
                            xy_c * xyv + yz_c * yzv + zx_c * zxv + x2_c * x2v +
                            y2_c * y2v + z2_c * z2v + x3_c * x3v + x4_c * x4v +
                            x2z_c * x2zv)
                    end
                end
            end
        end

        center_idx = (Potentials.x_axis_to_index(p, 0),
                      Potentials.y_axis_to_index(p, 0),
                      Potentials.z_axis_to_index(p, 0))

        for _ in 1:5
            coeffs = rand(11, neles)
            for (i, eid) in enumerate(eles)
                fill_potential!(coeffs[:, i], eid)
            end
            fitting = Potentials.Fitting(p, orders=(4, 2, 2), sizes=(41, 5, 5))
            eles2, coeffs2 = Solutions.compensate_terms(fitting, center_idx;
                                                        unit=unit, mask=mask, min_num=20)
            @test eles2 == eles

            @test coeffs2 ≈ coeffs[[dx, dy, dz, xy, yz, zx, z2, x2, x3, x4, x2z], :] atol=1e-4 rtol=1e-4

            eles3, terms3 = Solutions.solve_compensate(fitting, center_idx;
                                                       unit=unit, mask=mask,
                                                       min_num=20, minmax=true)
            @test eles3 == eles
            @test coeffs2 * hcat(terms3...) ≈ Matrix(I, nterms, nterms) atol=1e-6 rtol=1e-6

            eles4, terms4 = Solutions.solve_compensate(fitting, center_idx;
                                                       unit=unit, mask=mask,
                                                       min_num=20, minmax=false)
            @test eles4 == eles
            @test coeffs2 * hcat(terms4...) ≈ Matrix(I, nterms, nterms) atol=1e-6 rtol=1e-6

            tgt5 = rand(nterms)
            eles5, terms5 = Solutions.solve_target(fitting, center_idx, tgt5;
                                                  unit=unit, mask=mask,
                                                  min_num=20, minmax=true)
            @test eles5 == eles
            @test coeffs2 * terms5 ≈ tgt5 atol=1e-6 rtol=1e-6

            tgt6 = rand(nterms)
            eles6, terms6 = Solutions.solve_target(fitting, center_idx, tgt6;
                                                  unit=unit, mask=mask,
                                                  min_num=20, minmax=false)
            @test eles6 == eles
            @test coeffs2 * terms6 ≈ tgt6 atol=1e-6 rtol=1e-6

            @test_throws ArgumentError Solutions.solve_target(fitting, center_idx,
                                                              zeros(nterms + 1);
                                                              unit=unit, mask=mask,
                                                              min_num=20)
        end
    end
end
