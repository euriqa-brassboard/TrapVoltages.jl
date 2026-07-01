#

using TrapVoltages.Optimizers

using Test

@testset "minmax span" begin
    @test Optimizers.optimize_minmax_span([1 -1]', [0, 1]) ≈ [0.5, 0.5]
    @test Optimizers.optimize_minmax_span([1 2
                                           0 1
                                           1 0], [3, 0, 0]) ≈ [0.75, -0.75, -0.75]
    @test Optimizers.optimize_minmax_span(zeros(3, 0), [3, 0, 0],
                                          limited=[1 2
                                                   0 1
                                                   1 0]) ≈ [0.75, -0.75, -0.75]
    @test Optimizers.optimize_minmax_span(zeros(3, 0), [3, 0, 0],
                                          limited=[1 2
                                                   0 1
                                                   1 0] .* 0.5) ≈ [1.5, -0.5, -0.5]
end

@testset "minmax" begin
    @test Optimizers.optimize_minmax([1 -1], [1]) ≈ [0.5, -0.5]
    @test Optimizers.optimize_minmax([1 -1
                                      1 1], [1, 1]) ≈ [1, 0]

    @test Optimizers.optimize_minmax([1 -1 1], [1]) ≈ [1/3, -1/3, 1/3]
    @test Optimizers.optimize_minmax([1 -1 1], [1 2]) ≈ [1/3 2/3
                                                          -1/3 -2/3
                                                          1/3 2/3]
end
