#

using TrapVoltages: Units

using Test

@testset "TrapUnits" begin
    @test Units.SI == Units.TrapUnits(1, 1)
    @test Units.SI == Units.TrapUnits(1, 1)
    @test Units.euriqa == Units.Yb171

    @test Units.Yb171.V_unit_uV ≈ 525.3 atol=0.2
    @test Units.Yb171.l_unit_um ≈ 2.741 atol=0.002

    u171 = Units.MHz_unit(171e-3 / 6.022e23)
    @test Units.Yb171.V_unit ≈ u171.V_unit rtol=2e-3
    @test Units.Yb171.l_unit ≈ u171.l_unit rtol=2e-3

    u172 = Units.MHz_unit(172e-3 / 6.022e23)
    @test Units.Yb172.V_unit ≈ u172.V_unit rtol=2e-3
    @test Units.Yb172.l_unit ≈ u172.l_unit rtol=2e-3

    u174 = Units.MHz_unit(174e-3 / 6.022e23)
    @test Units.Yb174.V_unit ≈ u174.V_unit rtol=2e-3
    @test Units.Yb174.l_unit ≈ u174.l_unit rtol=2e-3

    u133 = Units.MHz_unit(133e-3 / 6.022e23)
    @test Units.Ba133.V_unit ≈ u133.V_unit rtol=2e-3
    @test Units.Ba133.l_unit ≈ u133.l_unit rtol=2e-3

    u137 = Units.MHz_unit(137e-3 / 6.022e23)
    @test Units.Ba137.V_unit ≈ u137.V_unit rtol=2e-3
    @test Units.Ba137.l_unit ≈ u137.l_unit rtol=2e-3

    u138 = Units.MHz_unit(138e-3 / 6.022e23)
    @test Units.Ba138.V_unit ≈ u138.V_unit rtol=2e-3
    @test Units.Ba138.l_unit ≈ u138.l_unit rtol=2e-3
end
