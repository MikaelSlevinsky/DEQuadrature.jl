export sincpade,padeval,polyroots

function sincpade(fv,phiv,n,r,s)
# This function computes the Pade coefficients of f by sampling at the Sinc nodes on the interval (-\infty,+\infty).
	k = [0:r+s+1].-floor((r+s)/2)
	xk = convert(Array{Float64},phiv[k.+n.+1])
	fvec = convert(Array{Float64},fv[k.+n.+1])
	PM = [xk.^([0:r]') -fvec.*xk.^([1:s]')]
	pq = PM\fvec
	return (pq[1:r+1],pq[r+2:end])
end

function padeval(p,q,x)
#This function computes the Pade approximant of numerator coefficients p and denomiator coefficients q.
	s = length(q)
	r = length(p)-1
	return sum(p'.*x.^([0:r]'),2)./(1.0.+sum(q'.*x.^([1:s]'),2))
end

function polyroots(q)
#This function computes the roots of the polynomial 1 + q[1] x + ... + q[s] x^s by the inverse eigenvalues of the companion matrix.
	s = length(q)
	M = diagm(ones(s-1),-1)
	M[1,:] -= q'
	D = eigvals(M)
	return 1.0./D
end
