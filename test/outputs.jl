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

        @test Dict(mf2) == Dict("AA"=>1, "BB"=>2, "CC"=>3, "DDDD"=>4)

        io = IOBuffer()
        write_file(io, mf)
        seek(io, 0)
        mf3 = load_file(io, MapFile)
        @test mf3.names == ["AA", "BB", "CC", "DDDD"]
    end
end
