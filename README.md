# DEQuadrature.jl

The aims of this `Julia` package are: to be the fastest general purpose quadrature package,
supporting the canonical interval, and semi-infinite and infinite domains in arithmetic up to
and including `BigFloat`; and, to provide support for maximizing convergence rates when complex
singularities are present near the contour of integration. Since the package can handle integrable
algebraic and logarithmic endpoint singularities, and since it allows the user to consider other
domains by declaring a new instance of the type `Domain`, the package is general purpose.

The primary function of this module computes the nodes and weights
of the trapezoidal rule, `dot(f(x),w)`, after a variable transformation induces
double exponential endpoint decay. In addition, the variable transformations
maximize the convergence rate despite complex singularities near the solution interval.

The secondary function of this module computes the parameters of the
conformal map `h(t)` in Eq. (3.14) of <a href="http://arxiv.org/abs/1406.3320">[1]</a>.
This module requires the use of the Julia package Ipopt for solving the nonlinear program.

```julia
using DEQuadrature
```

### Example 4.1 from <a href="http://arxiv.org/abs/1406.3320">[1]</a>

Suppose we are interested in calculating the integral of:

```julia
f(x) = exp(1./abs2(x-z[1]))./abs2(x-z[2])
```

on the interval `[-1,1]` with the singularities:

```julia
z = [complex(-0.5,1.0),complex(0.5,0.5)]
```

as well as a square root singularity at the left endpoint and a logarithmic singularity at the right endpoint. We use the package function `DEMapValues` to calculate the optimized map and the function `DENodesAndWeights` to calculate nodes and weights. Looping over a geometrically increasing order, we can approximate the integral very accurately:

```julia
u0,u,xpre = DEMapValues(z;domain=Finite(-0.5,0.0,0.0,1.0))
for i = 1:6
	x,w = DENodesAndWeights(u0,u,2^i;b2factor=0.5,domain=Finite(-0.5,0.0,0.0,1.0))
	val = dot(f(x),w)
	err = abs(val-convert(Float64,BigFloat(DEQuadrature.example4p1)))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end
```

### Example 4.2 from <a href="http://arxiv.org/abs/1406.3320">[1]</a>

The package has equal support for `BigFloat`s, making high precision calculations a breeze! Suppose we are interested in calculating the integral of:

```julia
f(x) = exp(10./abs2(x-z[1])).*cos(10./abs2(x-z[2]))./abs2(x-z[3])./abs(x-z[4])
```

on the real line with the singularities:

```julia
z = [complex(big(-2.0),1.0),complex(-1.0,.5),complex(1.0,0.25),complex(2.0,1.0)]
```

We start by setting the `digits` we desire:

```julia
DEQuadrature.digits(100)
```

Then, we use the package function `DEMapValues` to calculate the optimized map and the function `DENodesAndWeights` to calculate nodes and weights. Looping over a geometrically increasing order, we can approximate the integral very accurately:

```julia
u0,u,xpre = DEMapValues(z;domain=Infinite1)
for i = 1:10
	x,w = DENodesAndWeights(u0,u,2^i;domain=Infinite1)
	val = dot(f(x),w)
	err = abs(val-BigFloat(DEQuadrature.example4p2))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end
```

### Example 4.4 from <a href="http://arxiv.org/abs/1406.3320">[1]</a>

Suppose we are interested in calculating the integral of:

```julia
f(x) = x./abs(x-z[1])./abs2(x-z[2])./abs2(x-z[3])
```

on `[0,∞)` with the singularities:

```julia
z = [complex(big(1.0),1.0),complex(2.,.5),complex(3,1//3)]
```

We use the package functions `SincPade` and `PolyRoots` to compute the approximate locations of the singularities adaptively. Then, we use the package function `DENodesAndWeights` to calculate the optimized nodes and weights. Looping over a geometrically increasing order, we can approximate the integral very accurately:

```julia
x = zeros(BigFloat,5);
for i = 1:4
	x,w = DENodesAndWeights(Complex{BigFloat}[],2^i;domain=SemiInfinite2)
	val = dot(f(x),w)
	err = abs(val-BigFloat(DEQuadrature.example4p4))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end
for i = 5:8
	(p,q) = SincPade(f(x),x,(length(x)-1)/2,i-2,i+2)
	rootvec = PolyRoots(q)
	x,w = DENodesAndWeights(convert(Vector{Complex{BigFloat}},rootvec[end-4:2:end]),2^i;domain=SemiInfinite2,Hint=25)
	val = dot(f(x),w)
	err = abs(val-BigFloat(DEQuadrature.example4p4))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end
```


# References:


   1.	R. M. Slevinsky and S. Olver. <a href="http://arxiv.org/abs/1406.3320">On the use of conformal maps
		for the acceleration of convergence of the trapezoidal rule
		and Sinc numerical methods</a>, arXiv:1406.3320, to appear in SIAM SISC, 2014.
