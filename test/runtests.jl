using SincFun, DEQuadrature, Base.Test

println("Testing Example 4.1")

f(x) = exp(1./abs2(x-z[1]))./abs2(x-z[2])
z = [complex(-0.5,1.0),complex(0.5,0.5)]

h = DEMapValues(z;domain=Finite(-1.0,1.0,-0.5,0.0,0.0,1.0))
@time h = DEMapValues(z;domain=Finite(-1.0,1.0,-0.5,0.0,0.0,1.0))
@time for i = 1:6
	x,w = DENodesAndWeights(h,2^i;b2factor=0.5,domain=Finite(-1.0,1.0,-0.5,0.0,0.0,1.0))
	val = dot(f(x),w)
	err = abs(val-parse(Float64,DEQuadrature.example4p1))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end

@test_approx_eq h.u0 0.13912017017717093
@test_approx_eq h.u [0.19081242047193833,0.2193815287134606]
@test_approx_eq h.x [-1.9586438486791273,0.9642883741220253]

println("Testing Example 4.2")

DEQuadrature.digits(100)

f(x) = exp(10./abs2(x-z[1])).*cos(10./abs2(x-z[2]))./abs2(x-z[3])./abs(x-z[4])
z = [complex(big(-2.0),1.0),complex(-1.0,.5),complex(1.0,0.25),complex(2.0,1.0)]

@time h = DEMapValues(z;domain=Infinite1(BigFloat))
@time for i = 1:10
	x,w = DENodesAndWeights(h,2^i;domain=Infinite1(BigFloat))
	val = dot(f(x),w)
	err = abs(val-parse(BigFloat,DEQuadrature.example4p2))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end

@test_approx_eq Float64(h.u0) 5.771509561879045e-6
@test_approx_eq map(Float64,h.u) [0.25431456735504016,0.14935754761346481,-0.004543340891142737,9.988000449253805e-5]
@test_approx_eq map(Float64,h.x) [-9.061574715065113,-6.528358357285027,4.864488580186121,11.534505102993966]

println("Testing Example 4.4")

f(x) = x./abs(x-z[1])./abs2(x-z[2])./abs2(x-z[3])
z = [complex(big(1.0),1.0),complex(2.,.5),complex(3,1//3)]

x = zeros(BigFloat,5);
err = 0.0
for i = 1:4
	x,w = DENodesAndWeights(Complex{BigFloat}[],2^i;domain=SemiInfinite2(BigFloat))
	val = dot(f(x),w)
	err = abs(val-parse(BigFloat,DEQuadrature.example4p4))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end
@time for i = 5:9
	(p,q) = sincpade(f(x),x,div(length(x)-1,2),i-2,i+2);
	rootvec = polyroots(q);
	x,w = DENodesAndWeights(convert(Vector{Complex{BigFloat}},rootvec[end-4:2:end]),2^i;domain=SemiInfinite2(BigFloat),Hint=25)
	val = dot(f(x),w)
	err = abs(val-parse(BigFloat,DEQuadrature.example4p4))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end

@test err <= -3log(eps(BigFloat))*eps(BigFloat)

println("Testing γ > 1")
z = [complex(-0.5,1.0),complex(0.0,0.5),complex(0.5,0.75)]
ga=2.0
h = DEMapValues(z;ga=ga)
@test norm(z - tanh(h[complex(h.x,pi/2/ga)])) <= sqrt(eps())

println("Testing from Townsend, Trogdon and Olver, arXiv:1410.5286, 2014")

DEQuadrature.digits(56)

val = zeros(BigFloat,10)
z = [complex(big(0.0),big(0.2))]
f(x) = exp(-x.^8/2+cos(10x))./(1+25x.^2)

@time h = DEMapValues(z;ga=big(8.0),domain=Infinite2(BigFloat))
@time for i = 1:10
	x,w = DENodesAndWeights(h,2^i;b2factor=h.u0^7/2^9,ga=big(8.0),domain=Infinite2(BigFloat))
	val[i] = dot(f(x),w)
	println("Order: ",i," Value: ",val[i])
end
abs(val[10]-val[9]) <= -3log(eps(BigFloat))*eps(BigFloat)
