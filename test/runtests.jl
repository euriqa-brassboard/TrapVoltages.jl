using TrapVoltages
using TrapVoltages: PolyFit
using Test

@testset "PolyFit" begin
    @testset "fitresult_math" begin
        res1 = PolyFit.Result{1}((1,), zeros(2))
        @test all(res1.coefficient .== 0)
        res2 = PolyFit.Result{1}((1,), [1.0, 0.0])
        @test res2.coefficient == [1, 0]
        res3 = PolyFit.Result{1}((1,), [0.0, 2.0])
        @test res3.coefficient == [0, 2]

        @test (+res1).coefficient == res1.coefficient
        @test (-res1).coefficient == -(res1.coefficient)
        @test (-res2).coefficient == -(res2.coefficient)
        @test (-res3).coefficient == -(res3.coefficient)

        @test (res1 + res2).coefficient == res2.coefficient
        @test (res2 + res3).coefficient == [1, 2]

        @test (res1 - res2).coefficient == -(res2.coefficient)
        @test (res3 - res1).coefficient == res3.coefficient
        @test (res2 - res3).coefficient == [1, -2]

        @test (res2 * 2).coefficient == [2, 0]
        @test (5 * res3).coefficient == [0, 10]

        @test (res2 / 2).coefficient == [0.5, 0]
        @test (2 \ res2).coefficient == [0.5, 0]

        # Mismatched orders
        res4 = PolyFit.Result{1}((2,), zeros(3))
        @test_throws AssertionError res1 + res4
        @test_throws AssertionError res1 - res4
    end

    @testset "fitresult_eval" begin
        # 1 + 2x + x^3
        res1 = PolyFit.Result{1}((3,), [1.0, 2.0, 0.0, 1.0])
        for x in -2.0:0.25:2.0
            @test res1(x) ≈ 1 + 2x + x^3
        end
        @test res1[0] == 1
        @test res1[1] == 2
        @test res1[2] == 0
        @test res1[3] == 1
        res1[0] = 1.5
        res1[1] = 0
        res1[2] = -0.5
        res1[3] = 0.25
        for x in -2.0:0.25:2.0
            @test res1(x) ≈ 1.5 - 0.5x^2 + 0.25x^3
        end

        # x^2 + xy - 3x^2y^2 - y
        # Julia column-major: CartesianIndices((3,3)) iterates as
        # (1,1),(2,1),(3,1),(1,2),(2,2),(3,2),(1,3),(2,3),(3,3)
        # corresponding to order pairs (x,y):
        # (0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)
        res2 = PolyFit.Result{2}((2, 2), [0.0, 0.0, 1.0,  # x^n * y^0
                                          -1.0, 1.0, 0.0,  # x^n * y^1
                                          0.0, 0.0, -3.0]) # x^n * y^2
        for x in -2.0:0.25:2.0
            for y in -2.0:0.25:2.0
                @test res2(x, y) ≈ x^2 + x * y - 3x^2 * y^2 - y
            end
        end
        @test res2[0, 0] == 0
        @test res2[0, 1] == -1
        @test res2[0, 2] == 0
        @test res2[1, 0] == 0
        @test res2[1, 1] == 1
        @test res2[1, 2] == 0
        @test res2[2, 0] == 1
        @test res2[2, 1] == 0
        @test res2[2, 2] == -3

        res2[0, 0] = 1
        res2[0, 1] = 0
        res2[0, 2] = 2
        res2[1, 0] = -1
        res2[1, 1] = 0
        res2[1, 2] = 3
        res2[2, 0] = 0
        res2[2, 1] = 0.5
        res2[2, 2] = 0
        for x in -2.0:0.25:2.0
            for y in -2.0:0.25:2.0
                @test res2(x, y) ≈ 1 + 2y^2 - x + 3x * y^2 + 0.5x^2 * y
            end
        end
    end

    @testset "fitresult_shift" begin
        res1 = PolyFit.Result{1}((3,), [1.0, -1.0, 0.0, 1.0])
        for x in -2.0:0.25:2.0
            @test res1(x) ≈ 1 - x + x^3
        end
        res1_2 = PolyFit.shift(res1, (1.5,))
        for x in -2.0:0.25:2.0
            @test res1_2(x - 1.5) ≈ 1 - x + x^3 atol=1e-10
        end

        res2 = PolyFit.Result{2}((2, 2), [0.0, 0.0, 1.0,
                                          -1.0, 1.0, 0.0,
                                          0.0, 0.0, -3.0])
        for x in -2.0:0.25:2.0
            for y in -2.0:0.25:2.0
                @test res2(x, y) ≈ x^2 + x * y - 3x^2 * y^2 - y
            end
        end
        res2_2 = PolyFit.shift(res2, (1.5, -0.5))
        for x in -2.0:0.25:2.0
            for y in -2.0:0.25:2.0
                @test res2_2(x - 1.5, y + 0.5) ≈ x^2 + x * y - 3x^2 * y^2 - y atol=1e-10
            end
        end
        res2_3 = PolyFit.shift(res2, (1.5, 0.0))
        for x in -2.0:0.25:2.0
            for y in -2.0:0.25:2.0
                @test res2_3(x - 1.5, y) ≈ x^2 + x * y - 3x^2 * y^2 - y atol=1e-10
            end
        end
        res2_4 = PolyFit.shift(res2_3, (1.5, -1.25))
        for x in -2.0:0.25:2.0
            for y in -2.0:0.25:2.0
                @test res2_4(x - 3.0, y + 1.25) ≈ x^2 + x * y - 3x^2 * y^2 - y atol=1e-10
            end
        end
    end

    @testset "fitter" begin
        # 1D simple fit
        fitter1 = PolyFit.Fitter(2)
        # center = (3+1)/2 = 2, positions at indices 1,2,3 → x = -1, 0, 1
        x1 = [-1.0, 0.0, 1.0]
        y1 = @. 1.25 + x1 / 2 + x1^2 * 3
        res1 = fitter1 \ y1
        @test res1.coefficient ≈ [1.25, 0.5, 3.0] atol=1e-10

        # 1D with custom sizes (overdetermined)
        fitter1 = PolyFit.Fitter(4; sizes=(11,))
        x1 = collect(range(-5, 5, length=11))
        y1 = @. 1 - x1 * 2 + x1^2 / 3 + 0.4 * x1^4
        res1 = fitter1 \ y1
        @test res1.coefficient ≈ [1, -2, 1/3, 0, 0.4] atol=1e-10

        # 1D with custom center
        fitter1 = PolyFit.Fitter(4; sizes=(11,), center=(3.5,))
        x1 = collect(1.0:11.0) .- 3.5
        y1 = @. 1 - x1 * 2 + x1^2 / 3 + 0.4 * x1^4
        res1 = fitter1 \ y1
        @test res1.coefficient ≈ [1, -2, 1/3, 0, 0.4] atol=1e-10

        # 3D fit: -y^2 - 3z + y^2*z^3 - 3x^2*z^3 + y^2*z^4 - x^2*y^2*z^4
        fitter2 = PolyFit.Fitter(2, 2, 4)
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
        fitter2 = PolyFit.Fitter(2, 2, 4; sizes=(10, 21, 51))
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
        fitter2 = PolyFit.Fitter(2, 2, 4; sizes=(10, 21, 51), center=(4.4, 13.6, 31.0))
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

    @testset "fit_cache" begin
        function check_fit(fit_cache, pos; kwargs...)
            fit = get(fit_cache, pos; kwargs...)
            orders = fit_cache.fitter.orders
            for cidx in CartesianIndices(orders .+ 1)
                order = Tuple(cidx) .- 1
                fit_single = PolyFit.get_single(fit_cache, pos, order; kwargs...)
                fit_full = fit[order...]
                @test fit_single ≈ fit_full atol=1e-5
            end
            return fit
        end

        @testset "1D" begin
            fitter1 = PolyFit.Fitter(4; sizes=(11,))
            f1(x) = 1 - x - x * 2 + x^2 / 3 + 0.4 * x^4
            v1 = [f1(Float64(i)) for i in 1:101]
            fit_cache1 = PolyFit.FitCache(fitter1, v1)

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
            fitter2 = PolyFit.Fitter(4, 5; sizes=(11, 20))
            function f2(x, y)
                x = x - 15
                y = y - 50
                return (1 + x + y / 2 + x * 2 + x * y + x^2 / 3 +
                        (x / 5)^2 * (y / 10)^3 + (x / 5)^4 * (y / 10)^5)
            end
            v2 = [f2(Float64(i), Float64(j)) for i in 1:30, j in 1:100]
            fit_cache2 = PolyFit.FitCache(fitter2, v2)

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

    @testset "gradient" begin
        @testset "Result gradient" begin
            # 1 + 2x + x^3
            res1 = PolyFit.Result{1}((3,), [1.0, 2.0, 0.0, 1.0])
            for x in -2.0:0.25:2.0
                @test PolyFit.gradient(res1, 1, x) ≈ 2 + 3x^2
            end

            # x^2 + xy - 3x^2y^2 - y
            res2 = PolyFit.Result{2}((2, 2), [0.0, 0.0, 1.0,
                                              -1.0, 1.0, 0.0,
                                              0.0, 0.0, -3.0])
            for x in -2.0:0.25:2.0
                for y in -2.0:0.25:2.0
                    @test PolyFit.gradient(res2, 1, x, y) ≈ 2x + y - 6x * y^2
                    @test PolyFit.gradient(res2, 2, x, y) ≈ x - 6x^2 * y - 1
                end
            end
        end

        @testset "FitCache gradient 1D" begin
            fitter1 = PolyFit.Fitter(4; sizes=(11,))
            f1(x) = 1 - x - x * 2 + x^2 / 3 + 0.4 * x^4
            f1_dx(x) = -1 - 2 + 2x / 3 + 0.4 * 4 * x^3
            v1 = [f1(Float64(i)) for i in 1:101]
            fit_cache1 = PolyFit.FitCache(fitter1, v1)
            for pos in range(-20, 120, length=50)
                @test PolyFit.gradient(fit_cache1, 1, pos) ≈ f1_dx(pos) rtol=1e-5
            end
        end

        @testset "FitCache gradient 2D" begin
            fitter2 = PolyFit.Fitter(4, 5; sizes=(11, 20))
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
            fit_cache2 = PolyFit.FitCache(fitter2, v2)
            for xpos in range(-5, 35, length=50)
                for ypos in range(-10, 110, length=50)
                    @test PolyFit.gradient(fit_cache2, 1, xpos, ypos) ≈ f2_dx(xpos, ypos) rtol=1e-5
                    @test PolyFit.gradient(fit_cache2, 2, xpos, ypos) ≈ f2_dy(xpos, ypos) rtol=1e-5
                end
            end
        end
    end
end
