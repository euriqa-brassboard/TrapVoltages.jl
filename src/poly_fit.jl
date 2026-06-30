#

module PolyFit

# N dimensional polynomial fitting
"""
    Fitter{N}

Fitter for N-dimensional/N-variable polynomials. Contains precomputed coefficients
for fitting a N dimentional array of certain size.
"""
struct Fitter{N}
    orders::NTuple{N,Int}
    sizes::NTuple{N,Int}
    coefficient::Matrix{Float64}
    scales::Vector{Float64}
    # center is the origin of the polynomial in index (1-based)
    function Fitter(orders::Vararg{Integer,N};
                    sizes=orders .+ 1, center=(sizes .+ 1) ./ 2) where N
        @assert all(sizes .> orders)
        nterms = prod(orders .+ 1)
        npoints = prod(sizes)

        coefficient = Matrix{Float64}(undef, npoints, nterms)
        pos_lidxs = LinearIndices(sizes)
        pos_cidxs = CartesianIndices(sizes)
        ord_lidxs = LinearIndices(orders .+ 1)
        ord_cidxs = CartesianIndices(orders .+ 1)
        scales = Vector{Float64}(undef, nterms)
        scale_max = max.((sizes .- 1) ./ 2, 1.0)
        for iorder in ord_lidxs
            order = Tuple(ord_cidxs[iorder]) .- 1
            scales[iorder] = 1 / prod(scale_max.^order)
        end
        # Index for position
        for ipos in pos_lidxs
            # Position of the point, with the origin in the middle of the grid.
            pos = Tuple(pos_cidxs[ipos]) .- Float64.(center)
            # Index for the polynomial order
            for iorder in ord_lidxs
                order = Tuple(ord_cidxs[iorder]) .- 1
                coefficient[ipos, iorder] = prod(pos.^order) * scales[iorder]
            end
        end
        return new{N}(orders, sizes, coefficient, scales)
    end
end

struct Result{N}
    orders::NTuple{N,Int}
    coefficient::Vector{Float64}
end

function assert_same_orders(u::Result{N}, v::Result{N}) where N
    @assert u.orders == v.orders
end

Base.:+(u::Result) = u
function Base.:+(u::Result{N}, v::Result{N}) where N
    assert_same_orders(u, v)
    Result(u.orders, u.coefficient .+ v.coefficient)
end
Base.:-(u::Result) = Result(u.orders, .-u.coefficient)
function Base.:-(u::Result{N}, v::Result{N}) where N
    assert_same_orders(u, v)
    Result(u.orders, u.coefficient .- v.coefficient)
end

Base.:*(u::Result, v::Number) = Result(u.orders, u.coefficient .* v)
Base.:*(v::Number, u::Result) = u * v

Base.:/(u::Result, v::Number) = Result(u.orders, u.coefficient ./ v)
Base.:\(v::Number, u::Result) = u / v

function (res::Result{N})(pos::Vararg{Any,N}) where N
    sizes = res.orders .+ 1
    lindices = LinearIndices(sizes)
    cindices = CartesianIndices(sizes)
    v = 0.0
    for iorder in lindices
        order = Tuple(cindices[iorder]) .- 1
        v += res.coefficient[iorder] * prod(pos.^order)
    end
    return v
end

end
