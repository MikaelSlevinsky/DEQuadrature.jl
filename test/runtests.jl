using DEQuadrature
using Base.Test

println("Testing Example 4.1")

f(x) = exp(1./((x.-z[1].re).^2.+z[1].im.^2))./((x.-z[2].re).^2.+z[2].im.^2)
z = [complex(-0.5,1.0),complex(0.5,0.5)]

u0,u,xpre = DEMapValues(z;domain=Finite(-0.5,0.0,0.0,1.0))
for i = 1:6
	x,w = DENodesAndWeights(u0,u,2^i;domain=Finite(-0.5,0.0,0.0,1.0));
	val = dot(f(x),w)
	err = abs(val-convert(Float64,BigFloat(DEQuadrature.example4p1)))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end

@test (u0,u,xpre) == (0.139120170177171,[0.19081242047193833,0.21938152871346062],[-1.9586438486791269,0.9642883741220253])

println("Testing Example 4.2")

f(x) = exp(10./((x.-z[1].re).^2.+z[1].im.^2)).*cos(10./((x.-z[2].re).^2.+z[2].im.^2))./((x.-z[3].re).^2.+z[3].im.^2)./sqrt((x.-z[4].re).^2.+z[4].im.^2)
z = [complex(BigFloat("-2.0"),BigFloat("1.0")),complex(BigFloat("-1.0"),BigFloat("0.5")),complex(BigFloat("1.0"),BigFloat("0.25")),complex(BigFloat("2.0"),BigFloat("1.0"))]

u0,u,xpre = DEMapValues(z;digits=100,domain=Infinite)
for i = 1:10
	x,w = DENodesAndWeights(u0,u,2^i;digits=100,domain=Infinite);
	val = dot(f(x),w)
	err = abs(val-BigFloat(DEQuadrature.example4p2))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end

@test (u0,u,xpre) == (BigFloat("5.771509561879044849131044359058506643123109824955463409423828125e-06"),BigFloat[2.54314567355185927599592332626343704760074615478515625e-01,1.493575476134672286310234312622924335300922393798828125e-01,-4.54334089114293883382433847373249591328203678131103515625e-03,9.9880004492420689714109183210410947140189819037914276123046875e-05],BigFloat[-9.0615747150658449982074671424925327301025390625e+00,-6.5283583572857022403468363336287438869476318359375e+00,4.86448858018483409892951385700143873691558837890625e+00,1.15345051029942311515696928836405277252197265625e+01])

println("Testing Example 4.4")

f(x) = x./sqrt((x.-z[1].re).^2.+z[1].im.^2)./((x.-z[2].re).^2.+z[2].im.^2)./((x.-z[3].re).^2.+z[3].im.^2)
z = [complex(BigFloat("1.0"),BigFloat("1.0")),complex(BigFloat("2.0"),BigFloat("0.5")),complex(BigFloat("3.0"),BigFloat("1.0")/BigFloat("3.0"))]

x = zeros(BigFloat,5);
for i = 1:4
	x,w = DENodesAndWeights(Complex{BigFloat}[],2^i;digits=100,domain=SemiInfinite2);
	val = dot(f(x),w)
	err = abs(val-BigFloat(DEQuadrature.example4p4))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end
for i = 5:8
	(p,q) = SincPade(f(x),x,(length(x)-1)/2,i-2,i+2);
	rootvec = PolyRoots(q);
	x,w = DENodesAndWeights(convert(Vector{Complex{BigFloat}},rootvec[end-4:2:end]),2^i;digits=100,domain=SemiInfinite2,Hint=25);
	val = dot(f(x),w)
	err = abs(val-BigFloat(DEQuadrature.example4p4))
	println(@sprintf("Order: %2i Value: %19.16e Relative error: %6.2e",i,val,err))
end

@test x[2^8+1] == BigFloat("2.55307299802869725885222070506009540485755286623571924850837139733934172282612773853050817193850434842e+00")


println("Testing γ > 1")
z = [complex(-0.5,1.0),complex(0.0,0.5),complex(0.5,0.75)]
ga=2.0
u0,u,xpre = DEMapValues(z;ga=ga)
@test norm(z - tanh(DEQuadrature.h(complex(xpre,pi/2/ga),u0,u))) <= sqrt(eps())