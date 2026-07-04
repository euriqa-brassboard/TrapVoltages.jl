#

using Test
using TrapVoltages

@testset "opt deps" begin
    pkgid = Base.PkgId(TrapVoltages)
    flag = Ref(false)
    TrapVoltages.load_optdep(pkgid, flag, "")
    @test flag[]

    # This has the wrong UUID
    pkgid2 = Base.PkgId(Base.UUID("c2eab58b-df5a-4bb0-a5d5-4505d39c1d02"), "TrapVoltages")
    flag2 = Ref(false)
    @test_throws ArgumentError TrapVoltages.load_optdep(pkgid2, flag2, "")
    @test !flag2[]
    @test_throws ArgumentError TrapVoltages.load_optdep(pkgid2, flag2, "")
    @test !flag2[]
    flag2[] = true
    TrapVoltages.load_optdep(pkgid2, flag, "")
end
