module DEQuadrature
#
# A module by Richard Mikaël Slevinsky,
# Department of Mathematical & Statistical Sciences,
# University of Alberta, July 2014.
#
# The primary function of this module computes the nodes and weights
# of the trapezoidal rule, dot(f(x),w), after a variable transformation induces 
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
#		to accelerate the convergence of the trapezoidal rule
#		and Sinc numerical methods, arXiv:1406.3320, 2014.
#
using Ipopt

export DENodesAndWeights,DEMapValues
export Domain,Finite,Infinite,SemiInfinite1,SemiInfinite2

include("SincPade.jl")
export SincPade,PadeVal,PolyRoots

#
# The type Domain is used to select from the outer maps (3.9)--(3.12).
#
type Domain
	psi::Function
	psiinv::Function
	psip::Function
end

#For functions on a Finite domain with algebraic and logarithmic endpoint singularities, use the numbers α, β, γ, and δ to compute the singularities in a stable way.
Finite{T<:Number}(α::T,β::T,γ::T,δ::T) = Domain(t->tanh(t),t->atanh(t),t->sech(t).^2.*(2./(exp(2t).+1)).^α.*(2./(exp(-2t).+1)).^β.*log(2./(exp(2t).+1)).^γ.*log(2./(exp(-2t).+1)).^δ)
Infinite = Domain(t->sinh(t),t->asinh(t),t->cosh(t))
SemiInfinite1 = Domain(t->log(exp(t)+1),t->log(exp(t)-1),t->1./(1+exp(-t)))
SemiInfinite2 = Domain(t->exp(t),t->log(t),t->exp(t))

function DENodesAndWeights{T<:Number}(z::Array{Complex{T},1},n::Integer;ga::T=one(T),digits::Integer=77,domain::Domain=Finite(zero(T),zero(T),zero(T),zero(T)),Hint::Integer=10,obj_scaling_factor::Float64=-1.0)
	#
	# On entry:
	#
	# z is an array of singular points of the integrand, and
	# n determines the number of quadrature nodes and weights (2n+1).
	#
	# Optionally:
	#
	# digits controls the precision of BigFloat,
	# domain specifies the domain, with the default set to Finite,
	# Hint controls the homotopy solution process for the nonlinear program, and
	# obj_scaling_factor controls the scaling of the objective function.
	#
	# On return:
	#
	# x, w are the nodes and weights of the double exponential quadrature.
	#
	if T <: BigFloat
		bits = convert(Int64,ceil(digits*log2(10)))
		set_bigfloat_precision(bits)
	end

	u0,u,xpre = DEMapValues(z;ga=ga,digits=digits,domain=domain,Hint=Hint,obj_scaling_factor=obj_scaling_factor)
	
	b2opt = u0
	dDEopt = convert(T,pi)/2ga
	gaopt = ga
		
	hs = log(2*convert(T,pi)*dDEopt*gaopt*n/b2opt)/gaopt/n
	hsk=linspace(-hs*n,hs*n,2n+1);hhsk=hfast(hsk,u0,u)#hhsk=h(hsk,u0,u)

	x,w = domain.psi(hhsk),hs*domain.psip(hhsk).*hpfast(hsk,u0,u)#hp(hsk,u0,u)
	cutoff = !isinf(x).*!isnan(x).*!isinf(w).*!isnan(w)
	return x[cutoff],w[cutoff]
end

function DENodesAndWeights{T<:Number}(u0::T,u::Array{T,1},n::Integer;ga::T=one(T),digits::Integer=77,domain::Domain=Finite(zero(T),zero(T),zero(T),zero(T)),Hint::Integer=10,obj_scaling_factor::Float64=-1.0)
	#
	# On entry:
	#
	# u0,u are the predetermined map values, and
	# n determines the number of quadrature nodes and weights (2n+1).
	#
	# Optionally:
	#
	# digits controls the precision of BigFloat,
	# domain specifies the domain, with the default set to Finite,
	# Hint controls the homotopy solution process for the nonlinear program, and
	# obj_scaling_factor controls the scaling of the objective function.
	#
	# On return:
	#
	# x, w are the nodes and weights of the double exponential quadrature.
	#
	if T <: BigFloat
		bits = convert(Int64,ceil(digits*log2(10)))
		set_bigfloat_precision(bits)
	end
	
	b2opt = u0
	dDEopt = convert(T,pi)/2ga
	gaopt = ga
		
	hs = log(2*convert(T,pi)*dDEopt*gaopt*n/b2opt)/gaopt/n
	hsk=linspace(-hs*n,hs*n,2n+1);hhsk=hfast(hsk,u0,u)#hhsk=h(hsk,u0,u)

	x,w = domain.psi(hhsk),hs*domain.psip(hhsk).*hpfast(hsk,u0,u)#hp(hsk,u0,u)
	cutoff = !isinf(x).*!isnan(x).*!isinf(w).*!isnan(w)
	return x[cutoff],w[cutoff]
end

function DEMapValues{T<:Number}(z::Array{Complex{T},1};ga::T=one(T),digits::Integer=77,domain::Domain=Finite(zero(T),zero(T),zero(T),zero(T)),Hint::Integer=10,obj_scaling_factor::Float64=-1.0)
	#
	# On entry:
	#
	# z is an array of singular points of the integrand.
	#
	# Optionally:
	#
	# digits controls the precision of BigFloat,
	# domain specifies the domain, with the default set to Finite,
	# Hint controls the homotopy solution process for the nonlinear program, and
	# obj_scaling_factor controls the scaling of the objective function.
	#
	# On return:
	#
	# u0 and u are the parameters of the map h(t) in Eq. (3.14), and
	# x are the x-coordinates of the pre-images x +/- i pi/2γ of the singularities.
	#
	if T <: BigFloat
		bits = convert(Int64,ceil(digits*log2(10)))
		set_bigfloat_precision(bits)
	end

	psiinvz = domain.psiinv(z)
	global n = length(z)
	global dat = convert(Array{Float64,1},real(psiinvz))
	global ept = convert(Array{Float64,1},abs(imag(psiinvz)))
	global gaopt = ga

	x_U = [fill(30.0,n),fill(10.0,n)]
	x_L = -x_U
	
	g_U = zeros(Float64,2n)
	g_U[end] = 20.0
	g_L = -g_U

	prob = createProblem(2n, x_L, x_U, 2n, g_L, g_U, 4n^2, n*(2n+1),eval_f, eval_g, eval_grad_f, eval_jac_g,)#, eval_h)
	addOption(prob, "print_level", 0)
	addOption(prob, "obj_scaling_factor", obj_scaling_factor)
	addOption(prob, "hessian_approximation", "limited-memory")
	
	eptbar,mindex = findmin(ept)
	datbar = dat[mindex]

	datexact = dat
	eptexact = ept

	dattest = fill(datbar,n)
	epttest = ept
	
	prob.x = [sign(dat.-datbar).*acosh(ept./eptbar),datbar,zeros(Float64,n-1)]
	for j=0:Hint
		global dat = (1-j/Hint).*dattest.+j/Hint.*datexact
		global ept = (1-j/Hint).*epttest.+j/Hint.*eptexact
		status = solveProblem(prob)
		if Ipopt.ApplicationReturnStatus[status] != :Solve_Succeeded && Ipopt.ApplicationReturnStatus[status] != :Solved_To_Acceptable_Level
			#error("There was a problem with the convergence.\n Try a larger homotopy or a smaller objective scaling factor.\n The Ipopt return status is ",Ipopt.ApplicationReturnStatus[status],".")
		end
	end
	
	u0 = convert(T,prob.obj_val)
	u = convert(Array{T,1},prob.x[n+1:2n])
	x = convert(Array{T,1},prob.x[1:n])
	
	return u0,u,x
end

#
# These are auxiliary functions used to create the map h(t) in Eq. (3.14) and its derivative.
#
function h{T<:Real}(t::Union(T,Complex{T}),u0::T,u::Array{T,1})
	nu = length(u)
	u0*sinh(t) + dot(u,t.^[0:nu-1])
end
function h{T<:Real}(t::Union(Array{T,1},Array{Complex{T},1}),u0::T,u::Array{T,1})
	nu = length(u)
	u0*sinh(t) .+ t.^([0:nu-1]')*u
end
function h{T<:Real}(t,u0::T,u::Array{T,1})
	nu = length(u)
	ret = u0*sinh(t)
	for j=1:nu
		ret .+= u[j]*t.^(j-one(T))
	end
	return ret
end

function hp{T<:Real}(t::Union(T,Complex{T}),u0::T,u::Array{T,1})
	nu = length(u)
	u0*cosh(t) + dot(u[2:nu],[1:nu-1].*t.^[0:nu-2])
end
function hp{T<:Real}(t::Union(Array{T,1},Array{Complex{T},1}),u0::T,u::Array{T,1})
	nu = length(u)
	u0*cosh(t) .+ t.^([0:nu-2]')*([1:nu-1].*u[2:nu])
end
function hp{T<:Real}(t,u0::T,u::Array{T,1})
	nu = length(u)
	ret = u0*cosh(t)
	for j=2:nu
		ret .+= u[j]*(j-one(T))*t.^(j-2one(T))
	end
	return ret
end

function hfast{T<:Number}(t::Array{T,1},u0::T,u::Array{T,1})
	nu,ntm1d2 = length(u),int((length(t)-1)/2)
	schrec(ntm1d2,zero(T),u0*sinh(t[ntm1d2+2]),2cosh(t[ntm1d2+2])) .+ t.^([0:nu-1]')*u
end

function hpfast{T<:Number}(t::Array{T,1},u0::T,u::Array{T,1})
	nu,ntm1d2 = length(u),int((length(t)-1)/2)
	schrec(ntm1d2,u0,u0*cosh(t[ntm1d2+2]),2cosh(t[ntm1d2+2])) .+ t.^([0:nu-2]')*([1:nu-1].*u[2:nu])
end

function schrec{T<:Number}(nt::Integer,v1::T,v2::T,v3::T)
	rec = Array(T,2nt+1)
	rec[nt+1] = v1
	rec[nt+2] = v2
	rec[nt] = v3*v1-v2
	for i=2:nt
		rec[nt+i+1] = v3*rec[nt+i] - rec[nt+i-1]
		rec[nt-i+1] = v3*rec[nt-i+2] - rec[nt-i+3]
	end
	return rec
end

#
# The remainder of the functions are used by Ipopt.
#
function eval_f(x)
  temp1 = 0.0
  temp2 = 0.0
  for k=1:n
	temp3=0.0
	for j=1:n
		temp3+=x[n+j]*complex(x[k],pi/2)^(j-1)
	end
	temp1+=ept[k]-imag(temp3)
	temp2+=cosh(x[k])
  end
  return temp1/temp2
end

function eval_g(x, g)
  # Bad: g    = zeros(2)  # Allocates new array
  # OK:  g[:] = zeros(2)  # Modifies 'in place'
  g[:] = zeros(2n)
  f = eval_f(x)
  for k=1:n
	temp1=0.0
	for j=1:n
		temp1+=x[n+j]*complex(x[k],pi/2)^(j-1)
	end
	g[k] = real(temp1)-dat[k]
	g[n+k] = f - (ept[k] - imag(temp1))/cosh(x[k])
  end
  g[2n] = x[1] + x[n]
end

function eval_grad_f(x, grad_f)
  # Bad: grad_f    = zeros(4)  # Allocates new array
  # OK:  grad_f[:] = zeros(4)  # Modifies 'in place'
  
  grad_f[:] = zeros(2n)
  for r=1:n
	temp1=0.0
	temp2=0.0
	temp4=0.0
	temp5=0.0
	for k=1:n
		temp3=0.0
		for j=1:n
			temp3+=x[n+j]*complex(x[k],pi/2)^(j-1)
		end
		temp1+=ept[k]-imag(temp3)
		temp2+=cosh(x[k])
		temp4+=x[n+k]*(k-1)*complex(x[r],pi/2)^(k-2)
		temp5+=complex(x[k],pi/2)^(r-1)
	end
	grad_f[r] = -(temp2*imag(temp4) + sinh(x[r])*temp1)/temp2^2
	grad_f[n+r] = -imag(temp5)/temp2
  end
end

function eval_jac_g(x, mode, rows, cols, values)
  if mode == :Structure
    idx = 1
    for row = 1:2n
      for col = 1:2n
        rows[idx] = row
        cols[idx] = col
        idx += 1
      end
    end
  else
	grad_f = zeros(2n)
	values[:]=zeros(4n^2)
	eval_grad_f(x,grad_f)
	for k=1:n
	
		temp1=0.0
		temp2=0.0
		for r=1:n
			values[2n*(k-1)+r]=0.0
			values[2n^2+2n*(k-1)+r]=grad_f[r]
			values[2n^2+2n*(k-1)+n+r]=grad_f[n+r]
			
			temp1+=x[n+r]*(r-1)*complex(x[k],pi/2)^(r-2)
			values[2n*(k-1)+n+r] = real(complex(x[k],pi/2)^(r-1))
			
			temp2+=x[n+r]*complex(x[k],pi/2)^(r-1)
			
			temp3=0.0
			for j=1:n
				temp3+=x[n+j]*complex(x[k],pi/2)^(j-1)
			end
			values[2n^2+2n*(k-1)+n+r] += imag(complex(x[k],pi/2)^(r-1))/cosh(x[k])
		end
		values[2n*(k-1)+k] = real(temp1)
		values[2n^2+2n*(k-1)+k] += ((ept[k]-imag(temp2))*sinh(x[k])+cosh(x[k])*imag(temp1))/cosh(x[k])^2
	
	end
	for k =1:n
		values[4n^2-k+1] = 0.0
		values[4n^2-n-k+1] = 0.0
	end
	values[4n^2-2n+1] = 1.0
	values[4n^2-n] = 1.0
  end
end

#
# Below are 400-digit accurate results for the examples in the readme. To use, simply apply BigFloat to the strings to return the desired accuracy in the current BigFloat precision.
#
example4p1 = "-2.04645081160694748690442050179886173463698400851312978159495108281833923937599915411239665314732909414150138401559892277212580897054853305044117124435121918549124604645789190248721167573993297627891884242759929922470299331035480062586543433278436544374899377193299256405496172603459957078880800178958528771642962893695508317394167405216487683527096100475182610564645687173201068270992551244010304597e+00"
example4p2 =  "1.50133619876062770101030470326173553208854739646240081225845195322624377330867094144192464692844931351956249439446687259575908821074467954700173593866783934590906808530844595015500541680530415264580194453133658724574536164866695611699188173914168980057456962540582381978462005370937587010012338985572909021270993977527355156445163416743376725727256310809095229447802636632237845095430431700882602931e+01"
example4p4 =  "1.25561272649571457524072745777324565745811577731244208918556803079860974103212280154092098290346761945467623458393541899340910959432986144854308183981183773821087234316281865748406315573642170223984437961159361930528002471311495840420142211103620828044445220101847580610875452719310048279384346206258445806770989182361981757265589400465406224583252123717264572443063396429861634681780593162842211179e+01"

end #module



#This set of auxiliary functions does not work.
#This is the set I'm currently working on.
#=

#
# The remainder of the functions are used by Ipopt.
#
function eval_f(x)
  temp1 = 0.0
  temp2 = 0.0
  for k=1:n
	temp3=0.0
	for j=1:n
		temp3+=x[n+j]*complex(x[k],pi/2gaopt)^(j-1)
	end
	temp1+=ept[k]-imag(temp3)
	temp2+=cosh(x[k])*sin(pi/2gaopt)
  end
  println(temp1/temp2)
  return temp1/temp2
end

function eval_g(x, g)
  # Bad: g    = zeros(2)  # Allocates new array
  # OK:  g[:] = zeros(2)  # Modifies 'in place'
  g[:] = zeros(2n)
  f = eval_f(x)
  for k=1:n
	temp1=0.0
	for j=1:n
		temp1+=x[n+j]*complex(x[k],pi/2gaopt)^(j-1)
	end
	g[k] = f*sinh(x[k])*cos(pi/2gaopt) + real(temp1)-dat[k]
	g[n+k] =  f*cosh(x[k])*sin(pi/2gaopt) + imag(temp1)-ept[k]
	println(g[k],"        ",g[n+k])
  end
  g[2n] = x[1] + x[n]
end

function eval_grad_f(x, grad_f)
  # Bad: grad_f    = zeros(4)  # Allocates new array
  # OK:  grad_f[:] = zeros(4)  # Modifies 'in place'
  
  grad_f[:] = zeros(2n)
  for r=1:n
	temp1=0.0
	temp2=0.0
	temp4=0.0
	temp5=0.0
	for k=1:n
		temp3=0.0
		for j=1:n
			temp3+=x[n+j]*complex(x[k],pi/2gaopt)^(j-1)
		end
		temp1+=ept[k]-imag(temp3)
		temp2+=cosh(x[k])*sin(pi/2gaopt)
		temp4+=x[n+k]*(k-1)*complex(x[r],pi/2gaopt)^(k-2)
		temp5+=complex(x[k],pi/2gaopt)^(r-1)
	end
	grad_f[r] = -(temp2*imag(temp4) + sinh(x[r])*sin(pi/2gaopt)*temp1)/temp2^2
	grad_f[n+r] = -imag(temp5)/temp2
  end
end

function eval_jac_g(x, mode, rows, cols, values)
  if mode == :Structure
    idx = 1
    for row = 1:2n
      for col = 1:2n
        rows[idx] = row
        cols[idx] = col
        idx += 1
      end
    end
  else
	grad_f = zeros(2n)
	values[:]=zeros(4n^2)
	f = eval_f(x)
	eval_grad_f(x,grad_f)
	for k=1:n
	
		temp1=0.0
		temp2=0.0
		for r=1:n
			values[2n*(k-1)+r]=0.0
			values[2n^2+2n*(k-1)+r]=grad_f[r]*sinh(x[r])*cos(pi/2gaopt) + f*cosh(x[r])*cos(pi/2gaopt)
			values[2n^2+2n*(k-1)+n+r]=grad_f[n+r]*cosh(x[r])*sin(pi/2gaopt)
			
			temp1+=x[n+r]*(r-1)*complex(x[k],pi/2gaopt)^(r-2)
			values[2n*(k-1)+n+r] = real(complex(x[k],pi/2gaopt)^(r-1))
			
			temp2+=x[n+r]*complex(x[k],pi/2gaopt)^(r-1)
			
			temp3=0.0
			for j=1:n
				temp3+=x[n+j]*complex(x[k],pi/2gaopt)^(j-1)
			end
			values[2n^2+2n*(k-1)+n+r] += imag(complex(x[k],pi/2gaopt)^(r-1))
		end
		#values[2n*(k-1)+k] = real(temp1)
		#values[2n^2+2n*(k-1)+k] += ((ept[k]-imag(temp2))*sinh(x[k])+cosh(x[k])*imag(temp1))/cosh(x[k])^2
	
	end
	for k =1:n
		values[4n^2-k+1] = 0.0
		values[4n^2-n-k+1] = 0.0
	end
	values[4n^2-2n+1] = 1.0
	values[4n^2-n] = 1.0
  end
end
=#