#

import TrapVoltages as TV
import TrapVoltages: PolyFit as PF

using Test

@testset "Fitter construction" begin
    @testset "1D" begin
        f = PF.Fitter(3)
        @test f.orders == (3,)
        @test f.sizes == (4,)
        @test size(f.coefficient) == (4, 4)
    end
    @testset "2D" begin
        f = PF.Fitter(2, 3)
        @test f.orders == (2, 3)
        @test f.sizes == (3, 4)
        @test size(f.coefficient) == (12, 12)
    end
    @testset "custom sizes" begin
        f = PF.Fitter(2; sizes=(5,))
        @test f.orders == (2,)
        @test f.sizes == (5,)
        @test size(f.coefficient) == (5, 3)
    end
    @testset "custom center" begin
        f = PF.Fitter(2; sizes=(5,), center=(1.0,))
        @test f.orders == (2,)
        @test f.sizes == (5,)
    end
    @testset "sizes must exceed orders" begin
        @test_throws ArgumentError PF.Fitter(3; sizes=(3,))
        @test_throws ArgumentError PF.Fitter(3, 2; sizes=(3, 5))
    end
end

@testset "1D polynomial fitting and evaluation" begin
    # Fit a known polynomial: p(x) = 2 + 3x + 0.5x^2
    p = [2, 3, 0.5]
    f = PF.Fitter(2)
    data = [evalpoly(-1, p),
            evalpoly(0, p),
            evalpoly(1, p)]
    res = f \ data
    @test res isa PF.Result{1}
    @test res.orders == (2,)
    # Check coefficients via indexing
    @test res[0] ≈ p[1]
    @test res[1] ≈ p[2]
    @test res[2] ≈ p[3]
    # Evaluate at sample points
    for x in (-1.0, 0.0, 0.5, 1.0)
        @test res(x) ≈ evalpoly(x, p)
    end
end

@testset "2D polynomial fitting" begin
    p(x, y) = 1 + 2x + 3y + 4 * x * y
    f = PF.Fitter(1, 1)
    # sizes = (2, 2), center = (1.5, 1.5), so positions are ±0.5
    # (x, y) at indices: (1,1) -> (-0.5, -0.5), (1,2) -> (-0.5, 0.5),
    #                    (2,1) ->  (0.5, -0.5), (2,2) ->  (0.5, 0.5)
    data = [p(-0.5, -0.5) p(-0.5, 0.5);
            p(0.5, -0.5) p(0.5, 0.5)]
    res = f \ data
    @test res[0, 0] ≈ 1.0
    @test res[1, 0] ≈ 2.0
    @test res[0, 1] ≈ 3.0
    @test res[1, 1] ≈ 4.0
    @test res(0.0, 0.0) ≈ p(0, 0)
    @test res(1.0, 1.0) ≈ p(1, 1)
end

@testset "Result arithmetic" begin
    res1 = PF.Result{1}((1,), zeros(2))
    @test all(res1.coefficient .== 0)
    res2 = PF.Result{1}((1,), [1.0, 0.0])
    @test res2.coefficient == [1, 0]
    res3 = PF.Result{1}((1,), [0.0, 2.0])
    @test res3.coefficient == [0, 2]

    r_n1 = -res1
    r_n2 = -res2
    r_n3 = -res3
    @test (+res1).coefficient == res1.coefficient
    @test r_n1.coefficient == -(res1.coefficient)
    @test r_n2.coefficient == -(res2.coefficient)
    @test r_n3.coefficient == -(res3.coefficient)

    r_add12 = res1 + res2
    r_add23 = res2 + res3
    @test r_add12.coefficient == res2.coefficient
    @test r_add23.coefficient == [1, 2]

    r_sub12 = res1 - res2
    r_sub31 = res3 - res1
    r_sub23 = res2 - res3
    @test r_sub12.coefficient == -(res2.coefficient)
    @test r_sub31.coefficient == res3.coefficient
    @test r_sub23.coefficient == [1, -2]

    r_2mul2 = res2 * 2
    r_3mul5 = 5 * res3
    @test r_2mul2.coefficient == [2, 0]
    @test r_3mul5.coefficient == [0, 10]

    r_2div2 = res2 / 2
    r_2ldiv2 = 2 \ res2
    @test r_2div2.coefficient == [0.5, 0]
    @test r_2ldiv2.coefficient == [0.5, 0]

    for x in -1.0:0.5:1.0
        @test r_n1(x) ≈ -res1(x)
        @test r_n2(x) ≈ -res2(x)
        @test r_n3(x) ≈ -res3(x)

        @test r_add12(x) ≈ res1(x) + res2(x)
        @test r_add23(x) ≈ res2(x) + res3(x)

        @test r_sub12(x) ≈ res1(x) - res2(x)
        @test r_sub31(x) ≈ res3(x) - res1(x)
        @test r_sub23(x) ≈ res2(x) - res3(x)

        @test r_2mul2(x) ≈ res2(x) * 2
        @test r_3mul5(x) ≈ 5 * res3(x)

        @test r_2div2(x) ≈ res2(x) / 2
        @test r_2ldiv2(x) ≈ 2 \ res2(x)
    end

    # Mismatched orders
    res4 = PF.Result{1}((2,), zeros(3))
    @test_throws ArgumentError res1 + res4
    @test_throws ArgumentError res1 - res4
end

@testset "Result indexing (getindex/setindex!)" begin
    f = PF.Fitter(2)
    res = f \ [1.0, 0.0, 1.0]
    original_val = res[1]
    res[1] = 42.0
    @test res[1] == 42.0
    @test res[1] != original_val
end

@testset "gradient (Result)" begin
    @testset "1D gradient" begin
        p(x) = 2 + 3x + 0.5x^2
        pdx(x) = 3 + x
        f = PF.Fitter(2)
        data = [p(-1), p(0), p(1)]
        res = f \ data
        for x in -1.0:0.5:1.0
            @test TV.gradient(res, 1, x) ≈ pdx(x) atol=1e-12
        end
    end
    @testset "2D gradient" begin
        p(x, y) = 1 + 2x + 3y + 4 * x * y
        pdx(x, y) = 2 + 4y
        pdy(x, y) = 3 + 4x
        f = PF.Fitter(1, 1)
        data = [p(-0.5, -0.5) p(-0.5, 0.5)
                p(0.5, -0.5) p(0.5, 0.5)]
        res = f \ data
        for x in -2.0:0.25:2.0
            for y in -2.0:0.25:2.0
                @test TV.gradient(res, 1, x, y) ≈ pdx(x, y) atol=1e-12
                @test TV.gradient(res, 2, x, y) ≈ pdy(x, y) atol=1e-12
            end
        end
    end
end

@testset "shift" begin
    # p(x) = 2 + 3x + 0.5x^2
    # shifted by 1: p(x+1) = 2 + 3(x+1) + 0.5(x+1)^2
    #             = 2 + 3x + 3 + 0.5x^2 + x + 0.5
    #             = 5.5 + 4x + 0.5x^2
    f = PF.Fitter(2)
    data = [-0.5, 2.0, 5.5]
    res1 = f \ data
    shifted1 = PF.shift(res1, (1.0,))
    @test shifted1[0] ≈ 5.5 atol=1e-12
    @test shifted1[1] ≈ 4.0 atol=1e-12
    @test shifted1[2] ≈ 0.5 atol=1e-12
    # Verify shifted polynomial evaluates the same as original at shifted point
    for x in -2.0:0.5:2.0
        @test shifted1(x) ≈ res1(x + 1.0) atol=1e-10
    end

    res2 = PF.Result{2}((2, 2), [0.0, 0.0, 1.0,
                                 -1.0, 1.0, 0.0,
                                 0.0, 0.0, -3.0])
    for x in -2.0:0.25:2.0
        for y in -2.0:0.25:2.0
            @test res2(x, y) ≈ x^2 + x * y - 3x^2 * y^2 - y
        end
    end
    res2_2 = PF.shift(res2, (1.5, -0.5))
    for x in -2.0:0.25:2.0
        for y in -2.0:0.25:2.0
            @test res2_2(x - 1.5, y + 0.5) ≈ x^2 + x * y - 3x^2 * y^2 - y atol=1e-10
        end
    end
    res2_3 = PF.shift(res2, (1.5, 0.0))
    for x in -2.0:0.25:2.0
        for y in -2.0:0.25:2.0
            @test res2_3(x - 1.5, y) ≈ x^2 + x * y - 3x^2 * y^2 - y atol=1e-10
        end
    end
    res2_4 = PF.shift(res2_3, (1.5, -1.25))
    for x in -2.0:0.25:2.0
        for y in -2.0:0.25:2.0
            @test res2_4(x - 3.0, y + 1.25) ≈ x^2 + x * y - 3x^2 * y^2 - y atol=1e-10
        end
    end
end

@testset "shifted_coefficient" begin
    f = PF.Fitter(2)
    data = [-0.5, 2.0, 5.5]
    res = f \ data
    # shifted_coefficient gives the coefficient of the shifted polynomial
    # at a given order
    @test PF.shifted_coefficient(res, (1.0,), 0) ≈ 5.5
    @test PF.shifted_coefficient(res, (1.0,), 1) ≈ 4.0
    @test PF.shifted_coefficient(res, (1.0,), 2) ≈ 0.5
end

@testset "Fitter" begin
    # 1D simple fit
    fitter1 = PF.Fitter(2)
    x1 = [-1.0, 0.0, 1.0]
    y1 = 1.25 .+ x1 ./ 2 .+ x1.^2 .* 3
    res1 = fitter1 \ y1
    @test res1.coefficient ≈ [1.25, 0.5, 3.0] atol=1e-10

    # 1D with custom sizes (overdetermined)
    fitter1 = PF.Fitter(4; sizes=(11,))
    x1 = -5:5
    y1 = 1 .- x1 .* 2 .+ x1.^2 ./ 3 .+ 0.4 .* x1.^4
    res1 = fitter1 \ y1
    @test res1.coefficient ≈ [1, -2, 1/3, 0, 0.4] atol=1e-10

    # 1D with custom center
    fitter1 = PF.Fitter(4; sizes=(11,), center=(3.5,))
    x1 = (1.0:11.0) .- 3.5
    y1 = 1 .- x1 .* 2 .+ x1.^2 ./ 3 .+ 0.4 .* x1.^4
    res1 = fitter1 \ y1
    @test res1.coefficient ≈ [1, -2, 1/3, 0, 0.4] atol=1e-10

    # 3D fit: -y^2 - 3z + y^2*z^3 - 3x^2*z^3 + y^2*z^4 - x^2*y^2*z^4
    fitter2 = PF.Fitter(2, 2, 4)
    # sizes = (3, 3, 5), centers = (2, 2, 3)
    data = Array{Float64}(undef, 3, 3, 5)
    for i in 1:3, j in 1:3, k in 1:5
        x, y, z = i - 2.0, j - 2.0, k - 3.0
        data[i, j, k] = -y^2 - 3z + y^2 * z^3 - 3x^2 * z^3 + y^2 * z^4 - x^2 * y^2 * z^4
    end
    res2 = fitter2 \ data
    for i in 0:2, j in 0:2, k in 0:4
        c = res2[i, j, k]
        if (i, j, k) == (0, 2, 0)
            @test c ≈ -1 atol=1e-10
        elseif (i, j, k) == (0, 0, 1)
            @test c ≈ -3 atol=1e-10
        elseif (i, j, k) == (0, 2, 3)
            @test c ≈ 1 atol=1e-10
        elseif (i, j, k) == (2, 0, 3)
            @test c ≈ -3 atol=1e-10
        elseif (i, j, k) == (0, 2, 4)
            @test c ≈ 1 atol=1e-10
        elseif (i, j, k) == (2, 2, 4)
            @test c ≈ -1 atol=1e-10
        else
            @test c ≈ 0 atol=1e-10
        end
    end

    # 3D overdetermined fit
    fitter2 = PF.Fitter(2, 2, 4; sizes=(10, 21, 51))
    data = Array{Float64}(undef, 10, 21, 51)
    cx, cy, cz = (10 + 1) / 2, (21 + 1) / 2, (51 + 1) / 2
    for i in 1:10, j in 1:21, k in 1:51
        x, y, z = i - cx, j - cy, k - cz
        data[i, j, k] = -y^2 - 3z + y^2 * z^3 - 3x^2 * z^3 + y^2 * z^4 - x^2 * y^2 * z^4
    end
    res2 = fitter2 \ data
    for i in 0:2, j in 0:2, k in 0:4
        c = res2[i, j, k]
        if (i, j, k) == (0, 2, 0)
            @test c ≈ -1 atol=1e-5
        elseif (i, j, k) == (0, 0, 1)
            @test c ≈ -3 atol=1e-5
        elseif (i, j, k) == (0, 2, 3)
            @test c ≈ 1 atol=1e-5
        elseif (i, j, k) == (2, 0, 3)
            @test c ≈ -3 atol=1e-5
        elseif (i, j, k) == (0, 2, 4)
            @test c ≈ 1 atol=1e-5
        elseif (i, j, k) == (2, 2, 4)
            @test c ≈ -1 atol=1e-5
        else
            @test c ≈ 0 atol=1e-5
        end
    end

    # 3D overdetermined fit with custom center
    fitter2 = PF.Fitter(2, 2, 4; sizes=(10, 21, 51), center=(4.4, 13.6, 31.0))
    data = Array{Float64}(undef, 10, 21, 51)
    for i in 1:10, j in 1:21, k in 1:51
        x, y, z = i - 4.4, j - 13.6, k - 31.0
        data[i, j, k] = -y^2 - 3z + y^2 * z^3 - 3x^2 * z^3 + y^2 * z^4 - x^2 * y^2 * z^4
    end
    res2 = fitter2 \ data
    for i in 0:2, j in 0:2, k in 0:4
        c = res2[i, j, k]
        if (i, j, k) == (0, 2, 0)
            @test c ≈ -1 atol=1e-5
        elseif (i, j, k) == (0, 0, 1)
            @test c ≈ -3 atol=1e-5
        elseif (i, j, k) == (0, 2, 3)
            @test c ≈ 1 atol=1e-5
        elseif (i, j, k) == (2, 0, 3)
            @test c ≈ -3 atol=1e-5
        elseif (i, j, k) == (0, 2, 4)
            @test c ≈ 1 atol=1e-5
        elseif (i, j, k) == (2, 2, 4)
            @test c ≈ -1 atol=1e-5
        else
            @test c ≈ 0 atol=1e-5
        end
    end
end

@testset "FitCache" begin
    function check_fit(fit_cache, pos; kwargs...)
        fit = get(fit_cache, pos; kwargs...)
        orders = fit_cache.fitter.orders
        for cidx in CartesianIndices(orders .+ 1)
            order = Tuple(cidx) .- 1
            fit_single = TV.get_single(fit_cache, pos, order; kwargs...)
            fit_full = fit[order...]
            @test fit_single ≈ fit_full atol=1e-5
        end
        return fit
    end

    @testset "1D" begin
        fitter1 = PF.Fitter(4; sizes=(11,))
        f1(x) = 1 - x - x * 2 + x^2 / 3 + 0.4 * x^4
        v1 = f1.(1:101)
        fit_cache1 = PF.FitCache(fitter1, v1)

        function check_fit1(pos; kwargs...)
            fit = check_fit(fit_cache1, (pos,); kwargs...)
            @test fit[4] ≈ 0.4 atol=1e-6
            for x in range(-10, 10, length=50)
                @test fit(x) ≈ f1(x + pos) atol=0.1
            end
        end

        for pos in range(-20, 120, length=20)
            check_fit1(pos)
            for fit_center in range(-20, 120, length=20)
                check_fit1(pos; fit_center=(fit_center,))
            end
        end
    end
    @testset "2D" begin
        fitter2 = PF.Fitter(4, 5; sizes=(11, 20))
        function f2(x, y)
            x = x - 15
            y = y - 50
            return (1 + x + y / 2 + x * 2 + x * y + x^2 / 3 +
                (x / 5)^2 * (y / 10)^3 + (x / 5)^4 * (y / 10)^5)
        end
        v2 = [f2(i, j) for i in 1:30, j in 1:100]
        fit_cache2 = PF.FitCache(fitter2, v2)

        function check_fit2(xpos, ypos; kwargs...)
            fit = check_fit(fit_cache2, (xpos, ypos); kwargs...)
            @test fit[4, 5] ≈ 1 / 5^4 / 10^5 atol=1e-8
            for x in range(-2, 2, length=10)
                for y in range(-2, 2, length=10)
                    expected = f2(x + xpos, y + ypos)
                    @test fit(x, y) ≈ expected atol=1.5 rtol=2e-3
                end
            end
        end

        for xpos in range(-5, 35, length=7)
            for ypos in range(-10, 110, length=7)
                check_fit2(xpos, ypos)
                for xfit_center in range(-5, 35, length=7)
                    for yfit_center in range(-10, 110, length=7)
                        check_fit2(xpos, ypos;
                                   fit_center=(xfit_center, yfit_center))
                    end
                end
            end
        end
    end
end

@testset "get_single" begin
    # p(x) = 1 + 0.5x + 0.1x^2
    npoints = 20
    f = PF.Fitter(2; sizes=(5,))
    center = (npoints + 1) / 2
    data = [1.0 + 0.5*(i - center) + 0.1*(i - center)^2 for i in 1:npoints]
    cache = PF.FitCache(f, data)
    # get_single returns a specific coefficient at a position
    @test TV.get_single(cache, (center,), (0,)) ≈ 1.0 atol=1e-10
    @test TV.get_single(cache, (center,), (1,)) ≈ 0.5 atol=1e-10
    @test TV.get_single(cache, (center,), (2,)) ≈ 0.1 atol=1e-10
end

@testset "FitCache gradient 1D" begin
    fitter1 = PF.Fitter(4; sizes=(11,))
    f1(x) = 1 - x - x * 2 + x^2 / 3 + 0.4 * x^4
    f1_dx(x) = -1 - 2 + 2x / 3 + 0.4 * 4 * x^3
    v1 = [f1(Float64(i)) for i in 1:101]
    fit_cache1 = PF.FitCache(fitter1, v1)
    for pos in range(-20, 120, length=50)
        @test TV.gradient(fit_cache1, 1, pos) ≈ f1_dx(pos) rtol=1e-5
    end
end

@testset "FitCache gradient 2D" begin
    fitter2 = PF.Fitter(4, 5; sizes=(11, 20))
    function f2(x, y)
        x = x - 15
        y = y - 50
        return (1 + x + y / 2 + x * 2 + x * y + x^2 / 3 +
            (x / 5)^2 * (y / 10)^3 + (x / 5)^4 * (y / 10)^5)
    end
    function f2_dx(x, y)
        x = x - 15
        y = y - 50
        return (1 + 2 + y + 2x / 3 + 2 / 5 * (x / 5) * (y / 10)^3 +
            4 / 5 * (x / 5)^3 * (y / 10)^5)
    end
    function f2_dy(x, y)
        x = x - 15
        y = y - 50
        return (1 / 2 + x + 3 / 10 * (x / 5)^2 * (y / 10)^2 +
            5 / 10 * (x / 5)^4 * (y / 10)^4)
    end
    v2 = [f2(Float64(i), Float64(j)) for i in 1:30, j in 1:100]
    fit_cache2 = PF.FitCache(fitter2, v2)
    for xpos in range(-5, 35, length=50)
        for ypos in range(-10, 110, length=50)
            @test TV.gradient(fit_cache2, 1, xpos, ypos) ≈ f2_dx(xpos, ypos) rtol=1e-5
            @test TV.gradient(fit_cache2, 2, xpos, ypos) ≈ f2_dy(xpos, ypos) rtol=1e-5
        end
    end
end

@testset "FitCache caching" begin
    # Verify that repeated access returns the same result (cached)
    npoints = 20
    f = PF.Fitter(2; sizes=(5,))
    center = (npoints + 1) / 2
    data = [1.0 + 0.5*(i - center) + 0.1*(i - center)^2 for i in 1:npoints]
    cache = PF.FitCache(f, data)
    fit1 = get(cache, (center,))
    fit2 = get(cache, (center,))
    @test fit1[0] == fit2[0]
    @test fit1[1] == fit2[1]
    @test fit1[2] == fit2[2]
end
