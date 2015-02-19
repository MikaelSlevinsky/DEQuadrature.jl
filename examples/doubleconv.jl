# This file calculates the Chebyshev coefficients of the continuous function
# on (-1,1) whose 3-fold autoconvolution is constant on [-1,1], i.e. f*f*f = 1.
# f is even and has inverse square root endpoint singularities.

using ApproxFun, DEQuadrature

include("DEconv.jl")

T = Float64

n,ga = 50,convert(T,1.5)
h = log(2oftype(ga,π)*n)/ga/n
exp2sinhv,DATA = exp(sinh(h*[-n:n])*π),zeros(T,2n+1,2)
x,w = DENodesAndWeights(Complex{typeof(ga)}[],n;ga=ga)

a,b,α,β = -one(T),one(T),-one(T)/2,-one(T)/2
singularities{T<:Number}(a::T,b::T,α::T,β::T) = exp2sinhv.^α.*((b-a)./(exp2sinhv+1)).^(α+β)

f(c,x) = clenshaw(ApproxFun.interlace(c,zero(c)),x)
f3(c,x) = DEconv2(y->f(c,y),z->DEconv(y->f(c,y),a,b,α,β,z),a,b,α,β,x)
function ∂f3∂c!(c,z,DATA)
    DATA[:,1] = DEconv(y->f(c,y),a,b,α,β,(z-a)/2 + (a-z)/2*x)
    DATA[:,2] = DEconv(y->f(c,y),a,b,α,β,(z-b)/2 + (z-b)/2*x)
    ret = T[3DEconv2(y->cos((2j-2)*acos(y)),DATA,a,b,α,β,z) for j=1:length(c)]
end

N = 17
ceilN = div(N,2)+1
cf = pad(T[.5],ceilN)
#cf = pad(T[0.48238627247064414,-0.006602127192777267,-0.00013764938683027218,-3.568125070206212e-6,-9.636520776027288e-8,-2.6361598562903124e-9,-7.263559855979343e-11,-2.0118206133720452e-12,-5.420498916004464e-14],ceilN)
pts = chebyshevpoints(T,N)
halfpts = pts[ceilN:N]

function fgNL!(c,fvec,fjac,DATA)
    for i=1:ceilN
        fjac[i,:] = ∂f3∂c!(c,halfpts[i],DATA)
    end
    fvec[:] = fjac*c/3-1
end

@time println(f(cf,T[0]))
@time println(f3(cf,T[0]))

fvec,fjac,normupd = zeros(T,ceilN),zeros(T,ceilN,ceilN),one(T)
while normupd > sqrt(eps(T))
    @time fgNL!(cf,fvec,fjac,DATA)
    upd = fjac\fvec
    cf -= upd
    normupd = norm(upd)
    println("This is the updated coefficient vector: ",cf)
    println("This is the norm of the update: ",normupd)
end


xf = linspace(-one(T),one(T),101)
xf3 = linspace(-3one(T),3one(T),301)
fx = f(cf,xf)
f3x = f3(cf,xf3)
f3pts = f3(cf,pts)
K = f3(cf/cf[1]/π,zero(T))
