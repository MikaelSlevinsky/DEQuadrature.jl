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

Example


	using DEQuadrature


Suppose we are interested in calculating the integral of:


	f(x) = 1.0./((x.+0.5).^2.+1.0.^2)./((x.-0.5).^2.+0.5.^2)


on the real line to high accuracy (which equals 12pi/13). We start by recording the singularities:


	z = [complex(BigFloat("-0.5"),BigFloat("1.0")),complex(BigFloat("0.5"),BigFloat("0.5"))]


and we use the package function DEMapValues to calculate the optimized map and the function DENodesAndWeights to calculate nodes and weights. Looping over a geometrically increasing order, we can approximate the integral very accurately:


	u0,u,xpre = DEMapValues(z;digits=400,domain=Infinite)
	for i = 1:11
		x,w = DENodesAndWeights(u0,u,2^i;digits=400,domain=Infinite);
		println(dot(f(x),w)/pi*13/12)
	end


References:
 
   1.	R. M. Slevinsky and S. Olver. On the use of conformal maps
		to accelerate the convergence of the trapezoidal rule
		and Sinc numerical methods, arXiv:1406.3320, 2014.
