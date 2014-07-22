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

Example 4.1 from [1]


	using DEQuadrature


Suppose with the singularities:


	z = [complex(BigFloat("-0.5"),BigFloat("1.0")),complex(BigFloat("0.5"),BigFloat("0.5"))]


we are interested in calculating the integral of:


	f(x) = exp(1./((x.-z[1].re).^2.+z[1].im.^2))./((x.-z[2].re).^2.+z[2].im.^2)

with a square root singularity at the left endpoint, and a logarithmic singularity at the right endpoint. We use the package function DEMapValues to calculate the optimized map and the function DENodesAndWeights to calculate nodes and weights. Looping over a geometrically increasing order, we can approximate the integral very accurately:


	u0,u,xpre = DEMapValues(z;digits=100,domain=Finite(BigFloat("0.0"),BigFloat("-0.5"),BigFloat("1.0"),BigFloat("0.0")))
	for i = 1:8
		x,w = DENodesAndWeights(u0,u,2^i;digits=100,domain=Finite(BigFloat("0.0"),BigFloat("-0.5"),BigFloat("1.0"),BigFloat("0.0")));
		println(dot(f(x),w))
	end

Example 4.2 from [1]


	using DEQuadrature


Suppose with the singularities:


	z = [complex(BigFloat("-2.0"),BigFloat("1.0")),complex(BigFloat("-1.0"),BigFloat("0.5")),complex(BigFloat("1.0"),BigFloat("0.25")),complex(BigFloat("2.0"),BigFloat("1.0"))]


we are interested in calculating the integral of:


	f(x) = exp(10./((x.-z[1].re).^2.+z[1].im.^2)).*cos(10./((x.-z[2].re).^2.+z[2].im.^2))./((x.-z[3].re).^2.+z[3].im.^2)./sqrt((x.-z[4].re).^2.+z[4].im.^2)


We use the package function DEMapValues to calculate the optimized map and the function DENodesAndWeights to calculate nodes and weights. Looping over a geometrically increasing order, we can approximate the integral very accurately:


	u0,u,xpre = DEMapValues(z;digits=100,domain=Infinite)
	for i = 1:10
		x,w = DENodesAndWeights(u0,u,2^i;digits=100,domain=Infinite);
		println(dot(f(x),w))
	end


References:
 
   1.	R. M. Slevinsky and S. Olver. On the use of conformal maps
		to accelerate the convergence of the trapezoidal rule
		and Sinc numerical methods, arXiv:1406.3320, 2014.
