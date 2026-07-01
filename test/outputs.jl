#

using Test

using TrapVoltages.Outputs

@testset "MapFile" begin
    mktempdir() do d
        mf = MapFile(["AA", "BB", "CC", "DDDD"])
        write_file(joinpath(d, "map.txt"), mf)
        open(joinpath(d, "map.txt"), "a") do fh
            println(fh, "\n   # aaaa bbbb CCCC = 123")
        end
        mf2 = load_file(joinpath(d, "map.txt"), MapFile)
        @test mf2.names == ["AA", "BB", "CC", "DDDD"]
    end
end
