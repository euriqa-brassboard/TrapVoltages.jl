#

using TrapVoltages: Potentials, TrapDesc

using Test

const data_dir = joinpath(@__DIR__, "test_data")

function write_v0(fh, electrodes; nxyz, stride, origin)
    write(fh, Int32(0))
    write(fh, Int32(electrodes))
    nx, ny, nz = nxyz
    write(fh, Int32(nx))
    write(fh, Int32(ny))
    write(fh, Int32(nz))
    write(fh, Int32(0))
    write(fh, Float64(stride[1]))
    write(fh, Float64(stride[2]))
    write(fh, Float64(stride[3]))
    write(fh, Float64(origin[1]))
    write(fh, Float64(origin[2]))
    write(fh, Float64(origin[3]))
    write(fh, Int32(0))
    write(fh, Int32(0))
    for i in 1:electrodes
        write(fh, Int32(i - 1))
    end
    data = rand(nz, ny, nx, electrodes)
    write(fh, data)
    return data
end

function write_v1(fh, electrodes; nxyz, stride, origin)
    write(fh, Int32(0))
    write(fh, Int32(electrodes))
    nx, ny, nz = nxyz
    write(fh, Int32(nx))
    write(fh, Int32(ny))
    write(fh, Int32(nz))
    write(fh, Int32(0))
    write(fh, 1.0)
    write(fh, 0.0)
    write(fh, 0.0)
    write(fh, 0.0)
    write(fh, 1.0)
    write(fh, 0.0)
    write(fh, Float64(stride[1]))
    write(fh, Float64(stride[2]))
    write(fh, Float64(stride[3]))
    write(fh, Float64(origin[1]))
    write(fh, Float64(origin[2]))
    write(fh, Float64(origin[3]))
    write(fh, Int32(0))
    write(fh, Int32(0))
    for i in 1:electrodes
        write(fh, Int32(i - 1))
    end
    data = rand(nz, ny, nx, electrodes)
    write(fh, data)
    return data
end

function write_64(fh, electrodes; nxyz, stride, origin)
    write(fh, Int64(0))
    write(fh, Int64(electrodes))
    nx, ny, nz = nxyz
    write(fh, Int64(nx))
    write(fh, Int64(ny))
    write(fh, Int64(nz))
    write(fh, Int64(0))
    write(fh, 1.0)
    write(fh, 0.0)
    write(fh, 0.0)
    write(fh, 0.0)
    write(fh, 1.0)
    write(fh, 0.0)
    write(fh, Float64(stride[1]))
    write(fh, Float64(stride[2]))
    write(fh, Float64(stride[3]))
    write(fh, Float64(origin[1]))
    write(fh, Float64(origin[2]))
    write(fh, Float64(origin[3]))
    write(fh, Int64(0))
    write(fh, Int64(0))
    for i in 1:electrodes
        write(fh, Int64(i - 1))
    end
    data = rand(nz, ny, nx, electrodes)
    write(fh, data)
    return data
end

function test_raw_potential(electrodes; nxyz)
    mktempdir() do d
        for (ver, f_wr, f_im) in (("v0", write_v0, Potentials.import_pillbox_v0_raw),
                                  ("v1", write_v1, Potentials.import_pillbox_v1_raw),
                                  ("64", write_64, Potentials.import_pillbox_64_raw))
            stride = (rand(), rand(), rand()) .* 1e-5 .+ 1e-7
            origin = (rand(), rand(), rand()) .* 2e-4 .- 1e-4

            file = joinpath(d, "data_$(ver).bin")
            data = open(file, "w") do fh
                return f_wr(fh, electrodes; nxyz=nxyz, stride=stride, origin=origin)
            end
            p = open(f_im, file)
            @test p.electrodes == electrodes
            @test p.nx == nxyz[1]
            @test p.ny == nxyz[2]
            @test p.nz == nxyz[3]
            @test all(p.stride_um .≈ stride .* 1e6)
            @test all(p.origin_um .≈ origin .* 1e6)
            @test size(p.data) == (p.nz, p.ny, p.nx, electrodes)
            @test p.data == data

            for i in 1:10
                v = Potentials.x_index_to_axis(p, i)
                @test v ≈ (i - 1) * stride[1] * 1e6 + origin[1] * 1e6
                @test Potentials.x_axis_to_index(p, v) ≈ i

                v = Potentials.y_index_to_axis(p, i)
                @test v ≈ (i - 1) * stride[2] * 1e6 + origin[2] * 1e6
                @test Potentials.y_axis_to_index(p, v) ≈ i

                v = Potentials.z_index_to_axis(p, i)
                @test v ≈ (i - 1) * stride[3] * 1e6 + origin[3] * 1e6
                @test Potentials.z_axis_to_index(p, v) ≈ i
            end

            fsz = stat(file).size

            open(file, "a") do fh
                truncate(fh, fsz + 10)
            end
            @test_throws ArgumentError open(f_im, file)

            open(file, "a") do fh
                truncate(fh, fsz - 10)
            end
            @test_throws ArgumentError open(f_im, file)

            open(file, "a") do fh
                truncate(fh, fsz)
            end
            open(f_im, file)
        end
    end
end

@testset "Raw Potentials" begin
    test_raw_potential(92, nxyz=(5, 4, 3))
    test_raw_potential(96, nxyz=(2, 3, 7))
end

function test_trap_potential(trap; nxyz)
    trap_desc = TrapDesc(trap)
    q1_id = trap_desc.ele_indices["Q1"]
    q2_id = trap_desc.ele_indices["Q2"]
    l1_id = trap_desc.ele_indices["L1"]
    l2_id = trap_desc.ele_indices["L2"]

    mktempdir() do d
        for (ver, f_wr, f_im) in (("v0", write_v0, Potentials.import_pillbox_v0),
                                  ("v1", write_v1, Potentials.import_pillbox_v1),
                                  ("64", write_64, Potentials.import_pillbox_64))
            electrodes = length(TrapDesc(trap).ele_names)
            stride = (rand(), rand(), rand()) .* 1e-5 .+ 1e-7
            origin = (rand(), rand(), rand()) .* 2e-4 .- 1e-4

            file = joinpath(d, "data_$(ver).bin")
            data = open(file, "w") do fh
                return f_wr(fh, electrodes; nxyz=nxyz, stride=stride, origin=origin)
            end
            data = permutedims(data, (3, 2, 1, 4))
            for p in [f_im(file, trap=trap), open(fh->f_im(fh, trap=trap), file)]
                @test p.electrodes == electrodes
                @test p.nx == nxyz[1]
                @test p.ny == nxyz[2]
                @test p.nz == nxyz[3]
                @test all(p.stride_um .≈ stride .* 1e6)
                @test all(p.origin_um .≈ origin .* 1e6)
                @test size(p.data) == (p.nx, p.ny, p.nz, electrodes)
                @test p.data == data
                @test p.trap === trap_desc

                for i in 1:10
                    v = Potentials.x_index_to_axis(p, i)
                    @test v ≈ (i - 1) * stride[1] * 1e6 + origin[1] * 1e6
                    @test Potentials.x_axis_to_index(p, v) ≈ i

                    v = Potentials.y_index_to_axis(p, i)
                    @test v ≈ (i - 1) * stride[2] * 1e6 + origin[2] * 1e6
                    @test Potentials.y_axis_to_index(p, v) ≈ i

                    v = Potentials.z_index_to_axis(p, i)
                    @test v ≈ (i - 1) * stride[3] * 1e6 + origin[3] * 1e6
                    @test Potentials.z_axis_to_index(p, v) ≈ i
                end

                pot1 = Potentials.get_potential(p, Dict("Q1"=>0.2, l2_id=>-0.1))
                @test pot1 ≈ @view(data[:, :, :, q1_id]) .* 0.2 .+ @view(data[:, :, :, l2_id]) .* -0.1
            end

            file_p1 = joinpath(d, "data_p1_$(ver).bin")
            open(file_p1, "w") do fh
                return f_wr(fh, electrodes + 1; nxyz=nxyz, stride=stride, origin=origin)
            end
            @test_throws ArgumentError f_im(file_p1, trap=trap)

            file_m1 = joinpath(d, "data_m1_$(ver).bin")
            open(file_m1, "w") do fh
                return f_wr(fh, electrodes - 1; nxyz=nxyz, stride=stride, origin=origin)
            end
            @test_throws ArgumentError f_im(file_m1, trap=trap)

            p1 = f_im(file, trap=trap, aliases=Dict("Q1"=>"GND"))
            @test p1.electrodes == electrodes - 1
            @test p1.electrode_names[1] == ["GND", "Q1"]
            @test p1.nx == nxyz[1]
            @test p1.ny == nxyz[2]
            @test p1.nz == nxyz[3]
            @test size(p1.data) == (p1.nx, p1.ny, p1.nz, electrodes - 1)
            @test @view(p1.data[:, :, :, 1]) ≈ @view(data[:, :, :, 1]) .+ @view(data[:, :, :, q1_id])
            @test @view(p1.data[:, :, :, 2:q1_id - 1]) ≈ @view(data[:, :, :, 2:q1_id - 1])
            @test @view(p1.data[:, :, :, q1_id:end]) ≈ @view(data[:, :, :, q1_id + 1:end])

            p2 = f_im(file, trap=trap, aliases=Dict(q1_id=>1))
            @test p2.electrodes == electrodes - 1
            @test p2.electrode_names[1] == ["GND", "Q1"]
            @test p2.nx == nxyz[1]
            @test p2.ny == nxyz[2]
            @test p2.nz == nxyz[3]
            @test p2.data == p1.data

            p3 = f_im(file, trap=trap, electrode_names=[["Q1", "Q2"], ["L1", "L2"]])
            @test p3.electrodes == 2
            @test p3.nx == nxyz[1]
            @test p3.ny == nxyz[2]
            @test p3.nz == nxyz[3]
            @test @view(p3.data[:, :, :, 1]) == @view(data[:, :, :, q1_id]) .+ @view(data[:, :, :, q2_id])
            @test @view(p3.data[:, :, :, 2]) == @view(data[:, :, :, l1_id]) .+ @view(data[:, :, :, l2_id])

            @test_throws ArgumentError f_im(file, trap=trap, aliases=Dict("Q1"=>"GND", "Q2"=>"Q1"))
            @test_throws ArgumentError f_im(file, trap=trap, aliases=Dict("Q2"=>"Q1", "Q1"=>"GND"))
            @test_throws ArgumentError f_im(file, trap=trap, aliases=Dict(q1_id=>1, q2_id=>q1_id))
            @test_throws ArgumentError f_im(file, trap=trap, aliases=Dict("Q1"=>"GND"),
                                            electrode_names=[["Q1", "Q2"], ["L1", "L2"]])
        end
    end
end

@testset "Potentials" begin
    test_trap_potential("hoa", nxyz=(5, 4, 3))
    test_trap_potential("phoenix", nxyz=(2, 3, 7))
    test_trap_potential("peregrine", nxyz=(2, 1, 7))
end

@testset "Fitting" begin
    trap = TrapDesc("phoenix")
    q1_id = trap.ele_indices["Q1"]
    q2_id = trap.ele_indices["Q2"]
    q3_id = trap.ele_indices["Q3"]

    io = IOBuffer()
    nx, ny, nz = (20, 10, 7)
    write_64(io, length(trap.ele_names), nxyz=(nx, ny, nz),
             stride=(1e-3, 1e-3, 1e-3), origin=(-10e-3, -5e-3, 20e-3))
    seek(io, 0)
    p = Potentials.import_pillbox_64(io, trap=trap)
    for x in 1:nx
        for y in 1:ny
            for z in 1:nz
                p.data[x, y, z, q1_id] = (x - 0.5)^2 - (y + 0.4)^2 + (z + 1.2) - (x - 0.5)^2 * (z + 1.2) - 0.2 * (x - 0.5)^3 + 0.2
                p.data[x, y, z, q2_id] = -(x + 0.7)^2 + 0.5 * (y - 1.3)^2 - (z + 2.1) + (x + 0.7) * (z + 2.1)^2 * (y - 1.3) + (x + 0.7)^4
                p.data[x, y, z, q3_id] = 0.1 * (x - 0.5)^2 - (y + 0.4) + (z + 1.2)^2 - (x - 0.5)^2 * (z + 1.2)^2 - 0.2 * (x - 0.5)^4 - 0.2
            end
        end
    end

    fitting = Potentials.Fitting(p, orders=(4, 2, 2), sizes=(6, 4, 4))
    @test get(fitting, q1_id) === get(fitting, "Q1")
    @test get(fitting, "Q2") === get(fitting, "Q2")

    fit1 = get(fitting, q1_id, (0.5, -0.4, -1.2))
    for xo in 0:4
        for yo in 0:2
            for zo in 0:2
                c = fit1[xo, yo, zo]
                @test Potentials.get_single(fitting, "Q1", (0.5, -0.4, -1.2), (xo, yo, zo)) ≈ c atol=1e-10
                if (xo, yo, zo) == (2, 0, 0)
                    @test c ≈ 1
                elseif (xo, yo, zo) == (0, 2, 0)
                    @test c ≈ -1
                elseif (xo, yo, zo) == (0, 0, 1)
                    @test c ≈ 1
                elseif (xo, yo, zo) == (2, 0, 1)
                    @test c ≈ -1
                elseif (xo, yo, zo) == (3, 0, 0)
                    @test c ≈ -0.2
                elseif (xo, yo, zo) == (0, 0, 0)
                    @test c ≈ 0.2
                else
                    @test c ≈ 0 atol=1e-9
                end
            end
        end
    end

    fit2 = get(fitting, "Q2", (-0.7, 1.3, -2.1))
    for xo in 0:4
        for yo in 0:2
            for zo in 0:2
                c = fit2[xo, yo, zo]
                @test Potentials.get_single(fitting, q2_id, (-0.7, 1.3, -2.1), (xo, yo, zo)) ≈ c atol=1e-10
                if (xo, yo, zo) == (2, 0, 0)
                    @test c ≈ -1
                elseif (xo, yo, zo) == (0, 2, 0)
                    @test c ≈ 0.5
                elseif (xo, yo, zo) == (0, 0, 1)
                    @test c ≈ -1
                elseif (xo, yo, zo) == (1, 1, 2)
                    @test c ≈ 1
                elseif (xo, yo, zo) == (4, 0, 0)
                    @test c ≈ 1
                else
                    @test c ≈ 0 atol=1e-9
                end
            end
        end
    end

    fit3 = get(fitting, q3_id, (0.5, -0.4, -1.2))
    for xo in 0:4
        for yo in 0:2
            for zo in 0:2
                c = fit3[xo, yo, zo]
                @test Potentials.get_single(fitting, "Q3", (0.5, -0.4, -1.2), (xo, yo, zo)) ≈ c atol=1e-10
                if (xo, yo, zo) == (2, 0, 0)
                    @test c ≈ 0.1
                elseif (xo, yo, zo) == (0, 1, 0)
                    @test c ≈ -1
                elseif (xo, yo, zo) == (0, 0, 2)
                    @test c ≈ 1
                elseif (xo, yo, zo) == (2, 0, 2)
                    @test c ≈ -1
                elseif (xo, yo, zo) == (4, 0, 0)
                    @test c ≈ -0.2
                elseif (xo, yo, zo) == (0, 0, 0)
                    @test c ≈ -0.2
                else
                    @test c ≈ 0 atol=1e-9
                end
            end
        end
    end

    fit4 = Potentials.get_electrodes(fitting, [q3_id=>0.5, ("Q1", -2)], (0.5, -0.4, -1.2))
    for xo in 0:4
        for yo in 0:2
            for zo in 0:2
                c1 = fit1[xo, yo, zo]
                c3 = fit3[xo, yo, zo]
                c4 = fit4[xo, yo, zo]
                @test c4 ≈ -2 * fit1[xo, yo, zo] + 0.5 * fit3[xo, yo, zo] atol=1e-10
                @test Potentials.get_electrodes(fitting, Dict("Q1"=>-2, "Q3"=>0.5),
                                                (0.5, -0.4, -1.2), (xo, yo, zo)) ≈ c4 atol=1e-10
            end
        end
    end

    fit0 = Potentials.get_electrodes(fitting, [], (0.5, -0.4, -1.2))
    for xo in 0:4
        for yo in 0:2
            for zo in 0:2
                @test fit0[xo, yo, zo] == 0
            end
        end
    end
end

@testset "Short map" begin
    mktempdir() do d
        f = joinpath(d, "short.csv")
        open(f, "w") do io
            println(io, "A,B")
            println(io, "C,B")
            println(io, "XYZ,X")
        end
        @test Potentials.load_short_map(f) == Dict("A"=>"B", "C"=>"B", "XYZ"=>"X")
    end
end
