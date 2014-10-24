# DEQuadrature

A module by Richard Mikaël Slevinsky,
Department of Mathematical & Statistical Sciences,
University of Alberta, July 2014.

The primary function of this module computes the nodes and weights
of the trapezoidal rule, dot(f(x),w), after a variable transformation induces 
double exponential endpoint decay. In addition, the variable transformations
maximize the convergence rate despite singularities near the solution interval.

The secondary function of this module computes the parameters of the
conformal map h(t) in Eq. (3.14) of the reference [1]. This module requires
the use of the Julia package Ipopt for solving the nonlinear program.


	using DEQuadrature


### Example 4.1 from [1]

Suppose we are interested in calculating the integral of:


	f(x) = exp(1./((x.-z[1].re).^2.+z[1].im.^2))./((x.-z[2].re).^2.+z[2].im.^2)


on the interval [-1,1] with the singularities:


	z = [complex(-0.5,1.0),complex(0.5,0.5)]


as well as a square root singularity at the left endpoint and a logarithmic singularity at the right endpoint. We use the package function DEMapValues to calculate the optimized map and the function DENodesAndWeights to calculate nodes and weights. Looping over a geometrically increasing order, we can approximate the integral very accurately:


	u0,u,xpre = DEMapValues(z;domain=Finite(-0.5,0.0,0.0,1.0))
	for i = 1:6
		x,w = DENodesAndWeights(u0,u,2^i;domain=Finite(-0.5,0.0,0.0,1.0))
		val = dot(f(x),w)
		err = abs(val-convert(Float64,BigFloat(DEQuadrature.example4p1)))
		println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
	end


### Example 4.2 from [1]

The package has equal support for BigFloats, making high precision calculations a breeze! Suppose we are interested in calculating the integral of:


	f(x) = exp(10./((x.-z[1].re).^2.+z[1].im.^2)).*cos(10./((x.-z[2].re).^2.+z[2].im.^2))./((x.-z[3].re).^2.+z[3].im.^2)./sqrt((x.-z[4].re).^2.+z[4].im.^2)


on the real line with the singularities:


	z = [complex(BigFloat("-2.0"),BigFloat("1.0")),complex(BigFloat("-1.0"),BigFloat("0.5")),complex(BigFloat("1.0"),BigFloat("0.25")),complex(BigFloat("2.0"),BigFloat("1.0"))]


We use the package function DEMapValues to calculate the optimized map and the function DENodesAndWeights to calculate nodes and weights. Looping over a geometrically increasing order, we can approximate the integral very accurately:


	u0,u,xpre = DEMapValues(z;digits=100,domain=Infinite1)
	for i = 1:10
		x,w = DENodesAndWeights(u0,u,2^i;digits=100,domain=Infinite1)
		val = dot(f(x),w)
		err = abs(val-BigFloat(DEQuadrature.example4p2))
		println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
	end


### Example 4.4 from [1]

Suppose we are interested in calculating the integral of:


	f(x) = x./sqrt((x.-z[1].re).^2.+z[1].im.^2)./((x.-z[2].re).^2.+z[2].im.^2)./((x.-z[3].re).^2.+z[3].im.^2)


on [0,∞) with the singularities:


	z = [complex(BigFloat("1.0"),BigFloat("1.0")),complex(BigFloat("2.0"),BigFloat("0.5")),complex(BigFloat("3.0"),BigFloat("1.0")/BigFloat("3.0"))]


We use the package functions SincPade and PolyRoots to compute the approximate locations of the singularities adaptively. Then, we use the package function DEMapValues to calculate the optimized map and the function DENodesAndWeights to calculate nodes and weights. Looping over a geometrically increasing order, we can approximate the integral very accurately:


	x = zeros(BigFloat,5);
	for i = 1:4
		x,w = DENodesAndWeights(Complex{BigFloat}[],2^i;digits=100,domain=SemiInfinite2)
		val = dot(f(x),w)
		err = abs(val-BigFloat(DEQuadrature.example4p4))
		println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
	end
	for i = 5:8
		(p,q) = SincPade(f(x),x,(length(x)-1)/2,i-2,i+2);
		rootvec = PolyRoots(q);
		x,w = DENodesAndWeights(convert(Vector{Complex{BigFloat}},rootvec[end-4:2:end]),2^i;digits=100,domain=SemiInfinite2,Hint=25)
		val = dot(f(x),w)
		err = abs(val-BigFloat(DEQuadrature.example4p4))
		println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
	end




# References:

 
   1.	R. M. Slevinsky and S. Olver. On the use of conformal maps
		for the acceleration of convergence of the trapezoidal rule
		and Sinc numerical methods, arXiv:1406.3320, 2014.
