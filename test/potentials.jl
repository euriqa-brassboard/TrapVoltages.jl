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
            @test all(p.stride .≈ stride .* 1000)
            @test all(p.origin .≈ origin .* 1000)
            @test size(p.data) == (p.nz, p.ny, p.nx, electrodes)
            @test p.data == data

            for i in 1:10
                v = Potentials.x_index_to_axis(p, i)
                @test v ≈ (i - 1) * stride[1] * 1000 + origin[1] * 1000
                @test Potentials.x_axis_to_index(p, v) ≈ i

                v = Potentials.y_index_to_axis(p, i)
                @test v ≈ (i - 1) * stride[2] * 1000 + origin[2] * 1000
                @test Potentials.y_axis_to_index(p, v) ≈ i

                v = Potentials.z_index_to_axis(p, i)
                @test v ≈ (i - 1) * stride[3] * 1000 + origin[3] * 1000
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
            for p in [f_im(file, trap=trap), open(fh->f_im(fh, trap=trap), file)]
                @test p.electrodes == electrodes
                @test p.nx == nxyz[1]
                @test p.ny == nxyz[2]
                @test p.nz == nxyz[3]
                @test all(p.stride .≈ stride .* 1000)
                @test all(p.origin .≈ origin .* 1000)
                @test size(p.data) == (p.nz, p.ny, p.nx, electrodes)
                @test p.data == data
                @test p.trap === trap_desc

                for i in 1:10
                    v = Potentials.x_index_to_axis(p, i)
                    @test v ≈ (i - 1) * stride[1] * 1000 + origin[1] * 1000
                    @test Potentials.x_axis_to_index(p, v) ≈ i

                    v = Potentials.y_index_to_axis(p, i)
                    @test v ≈ (i - 1) * stride[2] * 1000 + origin[2] * 1000
                    @test Potentials.y_axis_to_index(p, v) ≈ i

                    v = Potentials.z_index_to_axis(p, i)
                    @test v ≈ (i - 1) * stride[3] * 1000 + origin[3] * 1000
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
            @test size(p1.data) == (p1.nz, p1.ny, p1.nx, electrodes - 1)
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
