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

@testset "CompensationFile" begin
    mktempdir() do d
        mf = MapFile(["AA", "BB", "CC", "DDDD"])
        mf2 = MapFile(["AA", "BB", "CC"])

        cf = CompensationFile(mf, ["X", "Y", "Z"],
                              [rand(4), rand(4), [rand(), rand(), 1.0, -1.0]])

        write_file(joinpath(d, "comp.txt"), cf)
        cf2 = load_file(joinpath(d, "comp.txt"), CompensationFile)

        @test cf2.map.names == mf.names
        @test cf2.term_names == cf.term_names
        @test cf2.term_values == cf.term_values

        io = IOBuffer()
        write_file(io, cf)
        seek(io, 0)
        cf3 = load_file(io, CompensationFile, mapfile=mf)

        @test cf3.map.names == mf.names
        @test cf3.term_names == cf.term_names
        @test cf3.term_values == cf.term_values

        @test_throws ArgumentError load_file(joinpath(d, "comp.txt"), CompensationFile,
                                             mapfile=mf2)

        write_file(joinpath(d, "comp_ele_mismatch.txt"),
                   CompensationFile(mf, ["X", "Y", "Z"], [rand(4), rand(4), rand(6)]))
        @test_throws ArgumentError load_file(joinpath(d, "comp_ele_mismatch.txt"),
                                             CompensationFile)

        write_file(joinpath(d, "comp_term_mismatch.txt"),
                   CompensationFile(mf, ["X", "Y", "Z"], [rand(4), rand(4)]))
        @test_throws ArgumentError load_file(joinpath(d, "comp_term_mismatch.txt"),
                                             CompensationFile)
    end
end

@testset "TransferFile" begin
    mktempdir() do d
        mf = MapFile(["AA", "BB", "CC", "DDDD"])
        mf2 = MapFile(["AA", "BB", "CC"])

        tf = TransferFile(mf, [rand(4), rand(4), rand(4), [rand(), rand(), 1.0, -1.0]])

        write_file(joinpath(d, "trans.txt"), tf)
        tf2 = load_file(joinpath(d, "trans.txt"), TransferFile)

        @test tf2.map.names == mf.names
        @test tf2.line_values == tf.line_values

        io = IOBuffer()
        write_file(io, tf)
        seek(io, 0)
        tf3 = load_file(io, TransferFile, mapfile=mf)

        @test tf3.map.names == mf.names
        @test tf3.line_values == tf.line_values

        @test_throws ArgumentError load_file(joinpath(d, "trans.txt"), TransferFile,
                                             mapfile=mf2)

        write_file(joinpath(d, "trans_ele_mismatch.txt"),
                   TransferFile(mf, [rand(4), rand(4), rand(6)]))
        @test_throws ArgumentError load_file(joinpath(d, "trans_ele_mismatch.txt"),
                                             TransferFile)
    end
end
