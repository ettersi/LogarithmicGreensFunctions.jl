module LogarithmicGreensFunctions

export greens,domain

csqrt(z) = sqrt(complex(z))

# Need to turn domains into abstract arrays to make the dot syntax work
# (e.g. greens(domain(-1,1),[-2,0,2]) )
abstract Domain{T} <: AbstractVector{T}
Base.size(::Domain) = (1,)
Base.getindex(D::Domain,i::Int) = D
Base.linearindexing{T<:Domain}(::Type{T}) = Base.LinearFast()


# One interval

domain(a,b) = Domain1I(promote(a,b)...)
immutable Domain1I{T} <: Domain{Domain1I{T}}
    a::T
    b::T
end

function greens(z)
    s = signbit(real(z)) ? -1 : 1
    return log(abs(z + s*csqrt(z^2-1)))
end
function greens(D::Domain1I,z) 
    a,b = D.a,D.b
    greens(2/(b-a)*(z - (b+a)/2))
end


# Two intervals

using FastGaussQuadrature
using StaticArrays

function domain(a,b,c,d) 
    σ = SVector(a,b,c,d)
    Domain2I(σ,compute_s(σ))
end
immutable Domain2I <: Domain{Domain2I}
    σ::SVector{4,Float64}
    s::Float64
end

# Estimate number of quadrature points needed to reach machine precision
# Actually, the convergence rate is 2*greens(...), but dropping the factor 
# 2 balances that we are not estimating the prefactor
nquad(D::Domain1I,σ) = ceil(Int,-log(eps(Float64))/minimum(greens.(D,σ)))

function compute_s(σ)
    q = nquad(domain(σ[2],σ[3]), SVector(σ[1],σ[4]))
    x,w = gaussjacobi(q, -0.5,-0.5)
    x = (σ[3]+σ[2])/2 + (σ[3]-σ[2])/2 * x
    denom = sqrt(-(x-σ[1]).*(x-σ[4]))
    return sum(w.*x./denom)/sum(w./denom)
end

function greens(D::Domain2I,z)
    σ = D.σ
    s = D.s

    i = indmin(abs(z - σ))
    q = nquad(domain(σ[i],z), SVector(σ[mod1(i-1,4)],σ[mod1(i+1,4)]))
    x,w = gaussjacobi(q,0,-0.5)::Tuple{Vector{Float64},Vector{Float64}}
    x = (z + σ[i])/2 + (z - σ[i])/2 * x
    w *= csqrt((z - σ[i])/2)

    return real(sum([w*(x-s)/prod([csqrt(x-σ[j]) for j = [1:i-1;i+1:4]]) for (x,w) in zip(x,w)]))
end

end # module
