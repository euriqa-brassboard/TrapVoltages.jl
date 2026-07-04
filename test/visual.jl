#

using TrapVoltages.Visual

using EzXML
using Test

const svg_ns = ["svg" => "http://www.w3.org/2000/svg"]

@testset "Visual" begin
    v1 = Visual.SVGInfo("phoenix")

    io = IOBuffer()
    show(io, MIME"image/svg+xml"(), v1)
    s = String(take!(io))
    @test !contains(s, "Title")

    t1 = Visual.set_title!(v1, "AAA BBB CCC DDD")
    t1["data-x-attr-111"] = "aabbccasxkjfs"
    @test t1["data-x-attr-111"] == "aabbccasxkjfs"
    io = IOBuffer()
    show(io, MIME"image/svg+xml"(), v1)
    s = root(parsexml(String(take!(io))))
    st1 = findfirst("//svg:text[@id='Title']", s, svg_ns)
    @test st1 !== nothing
    @test st1.content == "AAA BBB CCC DDD"
    @test st1["data-x-attr-111"] == "aabbccasxkjfs"

    @test Visual.set_title!(v1, "") === t1
    io = IOBuffer()
    show(io, MIME"image/svg+xml"(), v1)
    s = String(take!(io))
    @test !contains(s, "Title")
    @test !contains(s, "AAA BBB CCC DDD")
    @test !contains(s, "data-x-attr-111=\"aabbccasxkjfs\"")

    function ele_nodes(s, ele)
        nodes = findall("//svg:path[starts-with(@id, '$(ele)-') or @id='$(ele)']",
                        s, svg_ns)
        @test !isempty(nodes)
        return nodes
    end

    Visual.fill_electrodes!(v1, Dict("Q1"=>"red", "Q2"=>0.0, "Q3"=>(-2, 5, 300), "Q4"=>(-0.1, 1.1, 0.9), "Q5"=>-1.1, "Q6"=>1.1, "Q7"=>0.7, "Q8"=>-0.9, "S1"=>"blue"))
    io = IOBuffer()
    show(io, MIME"image/svg+xml"(), v1)
    s = root(parsexml(String(take!(io))))
    for n in ele_nodes(s, "Q1")
        @test contains(n["style"], "fill:red")
    end
    for n in ele_nodes(s, "Q2")
        @test contains(n["style"], "fill:#f7f7f7")
    end
    for n in ele_nodes(s, "Q3")
        @test contains(n["style"], "fill:#0005ff")
    end
    for n in ele_nodes(s, "Q4")
        @test contains(n["style"], "fill:#00ffe6")
    end
    for n in ele_nodes(s, "Q5")
        @test contains(n["style"], "fill:#053061")
    end
    for n in ele_nodes(s, "Q6")
        @test contains(n["style"], "fill:#67001f")
    end
    for n in ele_nodes(s, "Q7")
        @test contains(n["style"], "fill:#c43c3c")
    end
    for n in ele_nodes(s, "Q8")
        @test contains(n["style"], "fill:#134b87")
    end
    for n in ele_nodes(s, "S1")
        @test contains(n["style"], "fill:blue")
    end

    c1 = Visual.add_circle!(v1, 0, 0, 10)
    c1["id"] = "circle1-12345"
    @test c1["id"] == "circle1-12345"
    c2 = Visual.add_circle!(v1, -1000, -20, 13.5318)
    c2["id"] = "circle2-akljafs"
    @test c2["id"] == "circle2-akljafs"
    io = IOBuffer()
    show(io, MIME"image/svg+xml"(), v1)
    s = root(parsexml(String(take!(io))))
    sc1 = findfirst("//svg:circle[@id='circle1-12345']", s, svg_ns)
    @test sc1["cx"] == "297.524"
    @test sc1["cy"] == "38.803"
    @test sc1["r"] == "0.739"
    sc2 = findfirst("//svg:circle[@id='circle2-akljafs']", s, svg_ns)
    @test sc2["cx"] == "223.636"
    @test sc2["cy"] == "37.326"
    @test sc2["r"] == "1"

    l1 = Visual.add_line!(v1, [-10, 0, 10], [0.4, 1.1, -0.1])
    l1["id"] = "line-124u3"
    @test l1["id"] == "line-124u3"
    l2 = Visual.add_line!(v1, [-10, 0, 1, 2, 3, 4], [0.4, 0.12, 0.12, 0.12, -0.1, -0.1])
    l2["id"] = "line2-lkajsdf"
    @test l2["id"] == "line2-lkajsdf"
    io = IOBuffer()
    show(io, MIME"image/svg+xml"(), v1)
    s = root(parsexml(String(take!(io))))
    sl1 = findfirst("//svg:polyline[@id='line-124u3']", s, svg_ns)
    @test sl1["points"] == "296.785,84.88 297.524,62.2 298.262,100"
    sl2 = findfirst("//svg:polyline[@id='line2-lkajsdf']", s, svg_ns)
    @test sl2["points"] == "296.785,84.88 297.524,95.464 297.671,95.464 297.745,100 297.819,100"
end
