module DEQuadrature
#
# A module by Richard Mikaël Slevinsky,
# Department of Mathematical & Statistical Sciences,
# University of Alberta, July 2014.
#
# The primary function of this module computes the nodes and weights
# of the trapezoidal rule after a variable transformation induces
# double exponential endpoint decay. In addition, the variable transformations
# maximize the convergence rate despite singularities near the solution interval.
#
# The secondary function of this module computes the parameters of the
# conformal map h(t) in Eq. (3.14) of the reference [1]. This module requires
# the use of the Julia package Ipopt for solving the nonlinear program.
#
#
# References:
#
#   1.	R. M. Slevinsky and S. Olver. On the use of conformal maps
#		for the acceleration of convergence of the trapezoidal rule
#		and Sinc numerical methods, SIAM J. Sci. Comput., 37:A676--A700, 2015.
#
using SincFun, Ipopt

import SincFun: digits, singularities

export DENodesAndWeights, DEMapValues

include("sincpade.jl")
include("ipoptfunctions.jl")

function DENodesAndWeights{T<:Number}(h::ConformalMap{T},n::Int;b2factor::T=one(T),ga::T=one(T),domain::Domain{T}=Finite(T),Hint::Int=10,obj_scaling_factor::Float64=-1.0)
    #
    # On entry:
    #
    # A ConformalMap object that stores the data for
    # h(t) = u₀sinh(t) + u₁ + u₂t + u₃ t^2 + ⋯ + uₙt^(n-1).
    # u0 and u are the parameters of the map h(t) in Eq. (3.14), and
    # x are the x-coordinates of the pre-images x +/- i pi/2γ of the singularities, and
    #
    # n determines the number of quadrature nodes and weights (N = 2n+1).
    #
    # Optionally:
    #
    # b2factor controls the proportionality of b2opt with u0,
    # ga specifies the factor γ,
    # domain specifies the domain, with the default set to Finite,
    # Hint controls the homotopy solution process for the nonlinear program, and
    # obj_scaling_factor controls the scaling of the objective function.
    #
    # On return:
    #
    # x, w are the nodes and weights of the double exponential quadrature.
    #

    b2opt = h.u0*b2factor
    global gaopt = ga
    dDEopt = convert(T,pi)/2gaopt

    hs = log(2*convert(T,pi)*dDEopt*gaopt*n/b2opt)/gaopt/n
    hsk=linspace(-hs*n,hs*n,2n+1);hhsk=hfast(h,hsk)

    x,w = domain.ψ(hhsk),hs*domain.ψp(hhsk).*singularities(domain,hhsk).*hpfast(h,hsk)
    cutoff = !isinf(x).*!isnan(x).*!isinf(w).*!isnan(w)
    return x[cutoff],w[cutoff]
end

function DENodesAndWeights{T<:Number}(z::Vector{Complex{T}},n::Int;b2factor::T=one(T),ga::T=one(T),domain::Domain{T}=Finite(T),Hint::Int=10,obj_scaling_factor::Float64=-1.0)
    h = DEMapValues(z;ga=ga,domain=domain,Hint=Hint,obj_scaling_factor=obj_scaling_factor)
    x,w = DENodesAndWeights(h,n;b2factor=b2factor,ga=ga,domain=domain,Hint=Hint,obj_scaling_factor=obj_scaling_factor)
    return x,w
end

function DEMapValues{T<:Number}(z::Vector{Complex{T}};ga::T=one(T),domain::Domain{T}=Finite(T),Hint::Int=10,obj_scaling_factor::Float64=-1.0)
    #
    # On entry:
    #
    # z is a vector of singular points of the integrand.
    # If the integrand has no complex singularities, then
    # let z = Complex{T}[]::Vector{Complex{T}}, a vector of type Complex{T} with length 0.
    #
    # Optionally:
    #
    # domain specifies the domain, with the default set to Finite,
    # Hint controls the homotopy solution process for the nonlinear program, and
    # obj_scaling_factor controls the scaling of the objective function.
    #
    # On return:
    #
    # A ConformalMap object that stores the data for
    # h(t) = u₀sinh(t) + u₁ + u₂t + u₃ t^2 + ⋯ + uₙt^(n-1).
    # u0 and u are the parameters of the map h(t) in Eq. (3.14), and
    # x are the x-coordinates of the pre-images x +/- i pi/2γ of the singularities.
    #

    ψinvz = domain.ψinv(z)
    global n = length(z)
    global dat = convert(Vector{Float64},real(ψinvz))
    global ept = convert(Vector{Float64},abs(imag(ψinvz)))
    global gaopt = convert(Float64,ga)
    global spg = sinpi(1/2gaopt)
    global cpg = cospi(1/2gaopt)

    if n == 0
        u0,u,x = convert(T,pi)/2,zeros(T,1),T[]
        return ConformalMap(u0,u,x)
    end

    eptbar,mindex = findmin(ept)
    eptbar /= spg
    datbar = dat[mindex]

    if n == 1
        u0,u,x = convert(T,eptbar),[convert(T,datbar)],zeros(T,1)
        return ConformalMap(u0,u,x)
    end

    datexact,eptexact = dat,ept
    dattest,epttest = fill(datbar,n),ept

    x_U = [fill(30.0,n),fill(10.0,n)]
    x_L = -x_U

    g_U = zeros(Float64,2n)
    g_U[end] = 20.0#/gaopt#20.0
    g_L = -g_U

    prob = createProblem(2n, x_L, x_U, 2n, g_L, g_U, 4n^2, n*(2n+1),eval_f, eval_g, eval_grad_f, eval_jac_g,)#, eval_h)
    addOption(prob, "print_level", 0);addOption(prob, "obj_scaling_factor", obj_scaling_factor);addOption(prob, "hessian_approximation", "limited-memory")
    prob.x = [real(asinh(complex(sign(dat-datbar),ept)/eptbar)),datbar,zeros(Float64,n-1)]
    for j=0:Hint
        global dat = (1-j/Hint).*dattest.+j/Hint.*datexact
        global ept = (1-j/Hint).*epttest.+j/Hint.*eptexact
        status = solveProblem(prob)
        if Ipopt.ApplicationReturnStatus[status] != :Solve_Succeeded && Ipopt.ApplicationReturnStatus[status] != :Solved_To_Acceptable_Level
            error("There was a problem with the convergence.\n Try a larger homotopy or a smaller objective scaling factor.\n The Ipopt return status is ",Ipopt.ApplicationReturnStatus[status],".")
        end
    end

    u0,u,x = convert(T,prob.obj_val),convert(Vector{T},prob.x[n+1:2n]),convert(Vector{T},prob.x[1:n])

    return ConformalMap(u0,u,x)
end

#
# These are fast methods to evaluate a ConformalMap and its derivative
# at a vector of equally spaced points using recurrence relations for hyperbolic functions.
#

function hfast{S,T<:Number}(h::ConformalMap{S},t::Vector{T})
    nu,ntm1d2 = length(h.u),int((length(t)-1)/2)
    schrec(ntm1d2,zero(T),h.u0*sinh(t[ntm1d2+2]),2cosh(t[ntm1d2+2])) .+ t.^([0:nu-1]')*h.u
end
function hpfast{S,T<:Number}(h::ConformalMap{S},t::Vector{T})
    nu,ntm1d2 = length(h.u),int((length(t)-1)/2)
    schrec(ntm1d2,h.u0,h.u0*cosh(t[ntm1d2+2]),2cosh(t[ntm1d2+2])) .+ t.^([0:nu-2]')*([1:nu-1].*h.u[2:nu])
end

function schrec{T<:Number}(nt::Int,v1::T,v2::T,v3::T)
    rec = Array(T,2nt+1)
    rec[nt+1] = v1
    rec[nt+2] = v2
    rec[nt] = v3*v1-v2
    @inbounds for i=2:nt
        rec[nt+i+1] = v3*rec[nt+i] - rec[nt+i-1]
        rec[nt-i+1] = v3*rec[nt-i+2] - rec[nt-i+3]
    end
    return rec
end

#
# Below are 400-digit accurate results for the examples in the readme and tests.
# To use, simply apply BigFloat to the strings to return the desired accuracy
# in the current BigFloat precision.
#
example4p1 = "-2.04645081160694748690442050179886173463698400851312978159495108281833923937599915411239665314732909414150138401559892277212580897054853305044117124435121918549124604645789190248721167573993297627891884242759929922470299331035480062586543433278436544374899377193299256405496172603459957078880800178958528771642962893695508317394167405216487683527096100475182610564645687173201068270992551244010304597e+00"
example4p2 =  "1.50133619876062770101030470326173553208854739646240081225845195322624377330867094144192464692844931351956249439446687259575908821074467954700173593866783934590906808530844595015500541680530415264580194453133658724574536164866695611699188173914168980057456962540582381978462005370937587010012338985572909021270993977527355156445163416743376725727256310809095229447802636632237845095430431700882602931e+01"
example4p4 =  "1.25561272649571457524072745777324565745811577731244208918556803079860974103212280154092098290346761945467623458393541899340910959432986144854308183981183773821087234316281865748406315573642170223984437961159361930528002471311495840420142211103620828044445220101847580610875452719310048279384346206258445806770989182361981757265589400465406224583252123717264572443063396429861634681780593162842211179e+01"

end #module
