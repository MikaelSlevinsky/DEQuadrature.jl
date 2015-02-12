function DEconv{T<:Number}(f::Function,g::Function,a::T,b::T,α::T,β::T,z::T)
    val = zero(T)
    if z < zero(T)
        ϕv = z/2 + (z-2a)/2*x
        ϕpv = (z-2a)/2*w
        temp = f(ϕv).*g(z-ϕv).*ϕpv.*singularities(a,z-a,α,α).*((b-ϕv).*(b+ϕv-z)).^β
        temp = temp[!isinf(temp).*!isnan(temp)]
        val = sum(temp)
    elseif z > zero(T)
        ϕv = z/2 + (2b-z)/2*x
        ϕpv = (2b-z)/2*w
        temp = f(ϕv).*g(z-ϕv).*ϕpv.*singularities(z-b,b,β,β).*((ϕv-a).*(z-a-ϕv)).^α
        temp = temp[!isinf(temp).*!isnan(temp)]
        val = sum(temp)
    end
    return val
end
DEconv{T<:Number}(f::Function,g::Function,a::T,b::T,α::T,β::T,z::Vector{T}) = T[DEconv(f,g,a,b,α,β,zk) for zk in z]
DEconv{T<:Number}(f::Function,a::T,b::T,α::T,β::T,z::Union(T,Vector{T})) = DEconv(f,f,a,b,α,β,z)

function DEconv2{T<:Number}(f::Function,g::Function,a::T,b::T,α::T,β::T,z::T)
    val = zero(T)
    if z < -one(T)
        ϕv = (z-a)/2 + (z-3a)/2*x
        ϕpv = (z-3a)/2*w
        temp = f(ϕv).*g(z-ϕv).*ϕpv.*singularities(a,z-2a,α,zero(α)).*(b-ϕv).^β
        temp = temp[!isinf(temp).*!isnan(temp)]
        val = real(sum(temp))
    elseif z < one(T)
        ϕv = (a+z)/2 + (z-a)/2*x
        ϕpv = (z-a)/2*w
        temp = f(ϕv).*g(z-ϕv).*ϕpv.*singularities(a,z,α,zero(α)).*(b-ϕv).^β
        temp = temp[!isinf(temp).*!isnan(temp)]
        val = real(sum(temp))
        ϕv = (z+b)/2 + (b-z)/2*x
        ϕpv = (b-z)/2*w
        temp = f(ϕv).*g(z-ϕv).*ϕpv.*singularities(z,b,zero(β),β).*(-a+ϕv).^α
        temp = temp[!isinf(temp).*!isnan(temp)]
        val += real(sum(temp))
    else
        ϕv = (z-b)/2 + (3b-z)/2*x
        ϕpv = (3b-z)/2*w
        temp = f(ϕv).*g(z-ϕv).*ϕpv.*singularities(z-2b,b,zero(β),β).*(-a+ϕv).^α
        temp = temp[!isinf(temp).*!isnan(temp)]
        val = real(sum(temp))
    end
    return val
end
DEconv2{T<:Number}(f::Function,g::Function,a::T,b::T,α::T,β::T,z::Vector{T}) = T[DEconv2(f,g,a,b,α,β,zk) for zk in z]

function DEconv2{T<:Number}(f::Function,g::Matrix{T},a::T,b::T,α::T,β::T,z::T)
    val = zero(T)
    if -one(T) ≤ z ≤ one(T)
        ϕv = (a+z)/2 + (z-a)/2*x
        ϕpv = (z-a)/2*w
        temp = f(ϕv).*g[:,1].*ϕpv.*singularities(a,z,α,zero(α)).*(b-ϕv).^β
        temp = temp[!isinf(temp).*!isnan(temp)]
        val = real(sum(temp))
        ϕv = (z+b)/2 + (b-z)/2*x
        ϕpv = (b-z)/2*w
        temp = f(ϕv).*g[:,2].*ϕpv.*singularities(z,b,zero(β),β).*(-a+ϕv).^α
        temp = temp[!isinf(temp).*!isnan(temp)]
        val += real(sum(temp))
    end
    return val
end
