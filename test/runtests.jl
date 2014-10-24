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

@test (u0,u,xpre) == (0.13912017017964132,[0.19081241387776557,0.21938152342304876],[-1.9586438658470111,0.9642884274124068])

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

@test (u0,u,xpre) == (BigFloat("5.7715093553931392047519344234984828290180303156375885009765625e-06"),BigFloat[2.5431505642314800041958733345381915569305419921875e-01,1.4935755483790635889107534239883534610271453857421875e-01,-4.54334154322944862303135238335016765631735324859619140625e-03,9.98796175635166386987118247731132214539684355258941650390625e-05],BigFloat[-9.0615771916260126062070412444882094860076904296875e+00,-6.528360642803523461452641640789806842803955078125e+00,4.8644842749738881337862039799802005290985107421875e+00,1.15345060142379356449282568064518272876739501953125e+01])

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

@test x[2^8+1] == BigFloat("2.55307299568262473892690721206267000683898880351530600212408745088542760245128823939982928414699501625e+00")