using DEQuadrature
using Base.Test

println("Testing Example 4.1")

f(x) = exp(1./abs2(x-z[1]))./abs2(x-z[2])
z = [complex(-0.5,1.0),complex(0.5,0.5)]

u0,u,xpre = DEMapValues(z;domain=Finite(-0.5,0.0,0.0,1.0))
for i = 1:6
	x,w = DENodesAndWeights(u0,u,2^i;b2factor=0.5,domain=Finite(-0.5,0.0,0.0,1.0))
	val = dot(f(x),w)
	err = abs(val-convert(Float64,BigFloat(DEQuadrature.example4p1)))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end

@test (u0,u,xpre) == (0.13912017017717093,[0.19081242047193833,0.2193815287134606],[-1.9586438486791273,0.9642883741220253])

println("Testing Example 4.2")

DEQuadrature.digits(100)

f(x) = exp(10./abs2(x-z[1])).*cos(10./abs2(x-z[2]))./abs2(x-z[3])./abs(x-z[4])
z = [complex(big(-2.0),1.0),complex(-1.0,.5),complex(1.0,0.25),complex(2.0,1.0)]

u0,u,xpre = DEMapValues(z;domain=Infinite1)
for i = 1:10
	x,w = DENodesAndWeights(u0,u,2^i;domain=Infinite1)
	val = dot(f(x),w)
	err = abs(val-BigFloat(DEQuadrature.example4p2))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end

@test (u0,u,xpre) == (BigFloat("5.771509561879041460999255341857150369833107106387615203857421875e-06"),BigFloat[2.54314567355082010724487417974160052835941314697265625e-01,1.4935754761346553554091087789856828749179840087890625e-01,-4.54334089114279658649930837555075413547456264495849609375e-03,9.98800044925044036743522202215217475895769894123077392578125e-05],BigFloat[-9.0615747150653209729398440686054527759552001953125e+00,-6.52835835728521995946493916562758386135101318359375e+00,4.86448858018575069905864438624121248722076416015625e+00,1.1534505102994042857744716457091271877288818359375e+01])

println("Testing Example 4.4")

f(x) = x./abs(x-z[1])./abs2(x-z[2])./abs2(x-z[3])
z = [complex(big(1.0),1.0),complex(2.,.5),complex(3,1//3)]

x = zeros(BigFloat,5);
for i = 1:4
	x,w = DENodesAndWeights(Complex{BigFloat}[],2^i;domain=SemiInfinite2)
	val = dot(f(x),w)
	err = abs(val-BigFloat(DEQuadrature.example4p4))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end
for i = 5:8
	(p,q) = SincPade(f(x),x,(length(x)-1)/2,i-2,i+2);
	rootvec = PolyRoots(q);
	x,w = DENodesAndWeights(convert(Vector{Complex{BigFloat}},rootvec[end-4:2:end]),2^i;domain=SemiInfinite2,Hint=25)
	val = dot(f(x),w)
	err = abs(val-BigFloat(DEQuadrature.example4p4))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end

@test x[2^8+1] == BigFloat("2.55307299802806318558093347887260855495042920224216067365244613360324258988776205014979791343391738471e+00")


println("Testing γ > 1")
z = [complex(-0.5,1.0),complex(0.0,0.5),complex(0.5,0.75)]
ga=2.0
u0,u,xpre = DEMapValues(z;ga=ga)
@test norm(z - tanh(DEQuadrature.h(complex(xpre,pi/2/ga),u0,u))) <= sqrt(eps())

println("Testing from Townsend, Trogdon and Olver, arXiv:1410.5286, 2014")

DEQuadrature.digits(56)

val = zeros(BigFloat,10)
z = [complex(big(0.0),big(0.2))]
f(x) = exp(-x.^8/2+cos(10x))./(1+25x.^2)

u0,u,xpre = DEMapValues(z;ga=big(8.0),domain=Infinite2)
for i = 1:10
	x,w = DENodesAndWeights(u0,u,2^i;b2factor=u0^7/2^9,ga=big(8.0),domain=Infinite2)
	val[i] = dot(f(x),w)
	println("Order: ",i," Value: ",val[i])
end
abs(val[10]-val[9]) <= -3log(eps(BigFloat))*eps(BigFloat)