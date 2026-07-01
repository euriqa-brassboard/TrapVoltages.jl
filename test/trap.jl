#

using Test
using TrapVoltages
using TrapVoltages: find_electrodes

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

@testset "find_electrodes" begin
    px = TrapDesc("phoenix")
    pr = TrapDesc("peregrine")
    hoa = TrapDesc("hoa")

    @test px.ele_region_pos[1].inner == pr.ele_region_pos[1].inner
    @test px.ele_region_pos[1].outer == pr.ele_region_pos[1].outer

    found20 = false
    found21 = false
    for p in px.ele_region_pos[1].inner
        if p.name == "Q20"
            found20 = true
            @test p.left == -35
            @test p.right == 35
            @test p.up
        elseif p.name == "Q21"
            found21 = true
            @test p.left == -35
            @test p.right == 35
            @test !p.up
        end
    end
    @assert found20
    @assert found21

    found54 = false
    found55 = false
    for p in px.ele_region_pos[1].outer
        if p.name == "Q54"
            found54 = true
            @test p.left == -70
            @test p.right == 70
            @test p.up
        elseif p.name == "Q55"
            found55 = true
            @test p.left == -70
            @test p.right == 70
            @test !p.up
        end
    end
    @assert found54
    @assert found55

    @test_throws ArgumentError find_electrodes("hoa", hoa.ele_indices, 0)
    @test_throws ArgumentError find_electrodes("phoenix", px.ele_indices, 0, region=2)
    @test_throws ArgumentError find_electrodes(pr, pr.ele_indices, 0, region=0)

    e1 = find_electrodes("phoenix", px.ele_indices, 0)
    @test sort!([px.ele_names[e] for e in e1]) == ["Q20", "Q21", "Q54", "Q55"]
    e2 = find_electrodes("phoenix", px.ele_indices, 0, min_dist=5000, ignore_id=())
    @test length(e2) == 90
    @test !(1 in e2)

    e3 = find_electrodes("phoenix", px.ele_indices, 0, min_num=6, ignore_id=())
    @test sort!([px.ele_names[e] for e in e3]) == ["Q18", "Q19", "Q20", "Q21", "Q22", "Q23", "Q54", "Q55"]

    e4 = find_electrodes("phoenix", px.ele_indices, 1, min_num=6, ignore_id=())
    @test sort!([px.ele_names[e] for e in e4]) == ["Q20", "Q21", "Q22", "Q23", "Q54", "Q55"]

    e5 = find_electrodes("phoenix", px.ele_indices, 0, min_num=9, ignore_id=())
    @test sort!([px.ele_names[e] for e in e5]) == ["Q18", "Q19", "Q20", "Q21", "Q22", "Q23", "Q52", "Q53", "Q54", "Q55", "Q56", "Q57"]

    e6 = find_electrodes("phoenix", px.ele_indices, 0, min_dist=40, ignore_id=())
    @test sort!([px.ele_names[e] for e in e6]) == ["Q18", "Q19", "Q20", "Q21", "Q22", "Q23", "Q52", "Q53", "Q54", "Q55", "Q56", "Q57"]

    id_q18 = px.ele_indices["Q18"]
    id_q57 = px.ele_indices["Q57"]

    e7 = find_electrodes("phoenix", px.ele_indices, 0, min_dist=40, ignore_id=(id_q18, id_q57))
    @test sort!([px.ele_names[e] for e in e7]) == ["Q19", "Q20", "Q21", "Q22", "Q23", "Q52", "Q53", "Q54", "Q55", "Q56"]

    @test_throws ErrorException find_electrodes("phoenix", px.ele_indices, 0, min_num=100)
end
