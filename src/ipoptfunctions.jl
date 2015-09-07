#
# These functions are required in the use of Ipopt.
#
function eval_f(x)
    temp1 = 0.0
    temp2 = 0.0
    xpg = complex(x[1:n],pi/2gaopt)
    @inbounds for k=1:n
        temp3=0.0
        @inbounds for j=1:n
            temp3+=x[n+j]*xpg[k]^(j-1)
        end
        temp1+=ept[k]-imag(temp3)
        temp2+=cosh(x[k])*spg
    end
    return temp1/temp2
end

function eval_g(x, g)
    # Bad: g    = zeros(2)  # Allocates new array
    # OK:  g[:] = zeros(2)  # Modifies 'in place'
    g[:] = zeros(2n)
    f = eval_f(x)
    xpg = complex(x[1:n],pi/2gaopt)
    @inbounds for k=1:n
        temp1=0.0
        @inbounds for j=1:n
            temp1+=x[n+j]*xpg[k]^(j-1)
        end
        g[k] = f*sinh(x[k])*cpg + real(temp1)-dat[k]
        g[n+k] =  f*cosh(x[k])*spg + imag(temp1)-ept[k]
    end
    @inbounds g[2n] = x[1] + x[n]
end

function eval_grad_f(x, grad_f)
    # Bad: grad_f    = zeros(4)  # Allocates new array
    # OK:  grad_f[:] = zeros(4)  # Modifies 'in place'
    grad_f[:] = zeros(2n)
    xpg = complex(x[1:n],pi/2gaopt)
    sinhx,coshx = sinh(x[1:n]),cosh(x[1:n])
    @inbounds for r=1:n
        temp1=0.0
        temp2=0.0
        temp4=0.0
        temp5=0.0
        @inbounds for k=1:n
            temp3=0.0
            @inbounds for j=1:n
                temp3+=x[n+j]*xpg[k]^(j-1)
            end
            temp1+=ept[k]-imag(temp3)
            temp2+=coshx[k]*spg
            temp4+=x[n+k]*(k-1)*xpg[r]^(k-2)
            temp5+=xpg[k]^(r-1)
        end
        grad_f[r] = -(temp2*imag(temp4) + sinhx[r]*spg*temp1)/temp2^2
        grad_f[n+r] = -imag(temp5)/temp2
    end
end

function eval_jac_g(x, mode, rows, cols, values)
    if mode == :Structure
        idx = 1
        for row = 1:2n, col = 1:2n
            rows[idx] = row
            cols[idx] = col
            idx += 1
        end
    else
        grad_f = zeros(2n)
        values[:]=zeros(4n^2)
        f = eval_f(x)
        eval_grad_f(x,grad_f)
        xpg = complex(x[1:n],pi/2gaopt)
        sinhx,coshx = sinh(x[1:n]),cosh(x[1:n])
        @inbounds for k=1:n
            temp1=0.0
            @inbounds for r=1:n
                values[2n*(k-1)+r]=grad_f[r]*sinhx[k]*cpg
                values[2n*(k-1)+n+r]=grad_f[n+r]*sinhx[k]*cpg

                values[2n^2+2n*(k-1)+r]=grad_f[r]*coshx[k]*spg
                values[2n^2+2n*(k-1)+n+r]=grad_f[n+r]*coshx[k]*spg

                values[2n*(k-1)+n+r] += real(xpg[k]^(r-1))
                values[2n^2+2n*(k-1)+n+r] += imag(xpg[k]^(r-1))

                temp1+=x[n+r]*(r-1)*xpg[k]^(r-2)
            end

            values[2n*(k-1)+k] += f*coshx[k]*cpg + real(temp1)
            values[2n^2+2n*(k-1)+k] += f*sinhx[k]*spg + imag(temp1)
        end
        @inbounds for k =1:n
            values[4n^2-k+1] = 0.0
            values[4n^2-n-k+1] = 0.0
        end
        @inbounds values[4n^2-2n+1] = 1.0
        @inbounds values[4n^2-n] = 1.0
    end
end

function eval_h(x, mode, rows, cols, obj_factor, lambda, values)
    if mode == :Structure
        idx = 1
        for row = 1:2n, col = 1:row
            rows[idx] = row
            cols[idx] = col
            idx += 1
        end
    else
        values[:] = zeros(n*(2n+1))
        constraints = zeros(n*(2n+1))
        grad_f = zeros(2n)
        f = eval_f(x)
        eval_grad_f(x,grad_f)
        xpg = complex(x[1:n],pi/2gaopt)
        sinhxs,coshxs = sinh(x[1:n])*spg,cosh(x[1:n])*spg
        sinhxc,coshxc = sinh(x[1:n])*cpg,cosh(x[1:n])*cpg
        
        # The objective values broken into three seperate cases:
        # case 1) (r<=n, p<r) 
        # case 2) (r<=n , p=r)
        # case 3) (r>n , p<=n)
        
        # case 1) ∂^2 f / ∂x_r ∂x_p (r<=n, p<r)
        for r=1:n
            temp7=0.0 
            for p=1:(r-1)
                temp1=0.0
                temp2=0.0
                temp4=0.0
                for k=1:n
                    temp3=0.0
                    for j=1:n
                        temp3+=x[n+j]*xpg[k]^(j-1)
                    end # for j
                    temp1+=ept[k]-imag(temp3)
                    temp2+=cosh(x[k])*spg
                    temp4+=x[n+k]*(k-1)*xpg[r]^(k-2)
                    temp7+=x[n+k]*(k-1)*xpg[p]^(k-2)                    
                end # for k
            values[int(r*(r-1)/2+p)] = (imag(temp4)*sinhxs[p])/temp2^2 + sinhxs[r]*(temp2*imag(temp7)+2*temp1*sinhxs[p])/temp2^3 
            end
        end

        # case 2) ∂^2 f / ∂x_r^2 (r<=n , p=r)
        for r=1:n
            temp1=0.0
            temp2=0.0
            temp4=0.0
            temp6=0.0
            for k=1:n
                temp3=0.0
                for j=1:n
                    temp3+=x[n+j]*xpg[k]^(j-1)
                end # for j
                temp1+=ept[k]-imag(temp3)
                temp2+=cosh(x[k])*spg
                temp4+=x[n+k]*(k-1)*xpg[r]^(k-2)
                temp6+=x[n+k]*(k-1)*(k-2)*xpg[r]^(k-3)
            end # for k
        values[int(r*(r-1)/2+r)] = (imag(temp4)*sinhxs[r])/temp2^2 + sinhxs[r]*(temp2*imag(temp4)+2*temp1*sinhxs[r])/temp2^3 - (imag(temp6)/temp2 + coshxs[r]*(temp1/temp2^2))          
        end
        
        # case 3) ∂^2 f / ∂u_r ∂ x_p = ∂^2 f / ∂x_{n+r} ∂ x_p (r>n , p<=n)
        for r=n+1:2n
            for p=1:n
                temp2 = 0.0
                temp5 = 0.0
                for k=1:n
                    temp2+=cosh(x[k])*spg
                    temp5+=xpg[k]^(r-1)
                end # for k
            values[int(r*(r-1)/2+p)] = (-temp2*imag((r-1)*xpg[p]^(r-2)) + sinhxs[p]*imag(temp5))/temp2^2
            end
        end
        # case 4) ∂^2 f / ∂u_r ∂ u_p = ∂^2 f / ∂x_{n+r} ∂ x_{n+p} == 0 (r>n,p>n)
        
# THE CONSTRAINTS
        for k =1:n-1
            for r=1:2n
                for p=1:r
                    if r<=n
                        #a loop for temp6
                        temp6 = 0.0
                        for j=1:n
                            temp6+=x[n+j]*(j-1)*(j-2)*xpg[r]^(j-3)
                        end # for j
                        # ∂^2 g_{k} / ∂x_r ∂x_p  && ∂^2 g_{n+k} / ∂x_r ∂x_p
                        constraints[int(r*(r-1)/2+p)] += lambda[k]*(values[int(r*(r-1)/2+p)]*sinh(x[k])*cpg + (k==p?  grad_f[r]*cosh(x[p])*cpg : 0.0) + (k==r?  grad_f[p]*cosh(x[r])*cpg : 0.0) + (k==r==p?  f*sinh(x[p])*cpg+real(temp6) : 0.0) ) 
                        constraints[int(r*(r-1)/2+p)] += lambda[n+k]*( values[int(r*(r-1)/2+p)]*cosh(x[n+k])*spg + ((k+n)==p?  grad_f[r]*sinh(x[p])*spg : 0.0) + ((k+n)==r?  grad_f[p]*sinh(x[r])*spg : 0.0) + ((k+n)==r==p?  f*cosh(x[p])*spg+imag(temp6) : 0.0) ) 
                    elseif r>n && p<n 
                        # ∂^2 g_{k} / ∂u_r ∂ x_p = ∂^2 g_{k} / ∂x_{n+r} ∂ x_p   &&  ∂^2 g_{n+k} / ∂x_{n+r} ∂x_p
                        constraints[int(r*(r-1)/2+p)] += lambda[k]*(values[int(r*(r-1)/2+p)]*sinh(x[k])*cpg + (k==p?  grad_f[r]*cosh(x[p])*cpg + real((r-1)*xpg[p]^(r-2)) : 0.0) ) 
                        constraints[int(r*(r-1)/2+p)] += lambda[n+k]*( values[int(r*(r-1)/2+p)]*cosh(x[(k+n)])*spg + ((k+n)==p?  grad_f[r]*sinh(x[p])*spg + imag((r-1)*xpg[p]^(r-2)) : 0.0) ) 
                    else
                        # ∂^2 g_{k} / ∂u_r ∂ u_p = ∂^2 g_{k} / ∂x_{n+r} ∂ x_{n+p}  &&    ∂^2 g_{n+k} / ∂x_{n+r} ∂ x_{n+p} 
                        constraints[int(r*(r-1)/2+p)] += lambda[k]*0.0
                        constraints[int(r*(r-1)/2+p)] += lambda[n+k]*0.0
                    end   # if loop
                end # for p
            end # for r
        end # for k
# case k = n
        for r=1:2n
                for p=1:r
                    if r<=n
                        #a loop for temp6
                        temp6 = 0.0
                        for j=1:n
                            temp6+=x[n+j]*(j-1)*(j-2)*xpg[r]^(j-3)
                        end # for j
                        # ∂^2 g_{k} / ∂x_r ∂x_p  && ∂^2 g_{n+k} / ∂x_r ∂x_p
                        constraints[int(r*(r-1)/2+p)] += lambda[n]*(values[int(r*(r-1)/2+p)]*sinh(x[n])*cpg + (n==p?  grad_f[r]*cosh(x[p])*cpg : 0.0) + (n==r?  grad_f[p]*cosh(x[r])*cpg : 0.0) + (n==r==p?  f*sinh(x[p])*cpg+real(temp6) : 0.0) ) 
                        constraints[int(r*(r-1)/2+p)] += lambda[2n]*(0.0)
                    elseif r>n && p<n 
                        # ∂^2 g_{k} / ∂u_r ∂ x_p = ∂^2 g_{k} / ∂x_{n+r} ∂ x_p   &&  ∂^2 g_{n+k} / ∂x_{n+r} ∂x_p
                        constraints[int(r*(r-1)/2+p)] += lambda[n]*(values[int(r*(r-1)/2+p)]*sinh(x[n])*cpg + (n==p?  grad_f[r]*cosh(x[p])*cpg + real((r-1)*xpg[p]^(r-2)) : 0.0) ) 
                        constraints[int(r*(r-1)/2+p)] += lambda[2n]*(0.0)
                    else
                        # ∂^2 g_{k} / ∂u_r ∂ u_p = ∂^2 g_{k} / ∂x_{n+r} ∂ x_{n+p}  &&    ∂^2 g_{n+k} / ∂x_{n+r} ∂ x_{n+p} 
                        constraints[int(r*(r-1)/2+p)] += lambda[n]*0.0
                        constraints[int(r*(r-1)/2+p)] += lambda[2n]*(0.0)
                    end   # if loop
                end # for p
            end # for r
            
        values[:] = obj_factor*values[:] + constraints[:]
    end # if loop
end # function

#= Fun example
using SincFun, DEQuadrature
z = [complex((-2.0),1.0),complex(-1.0,.5),complex(1.0,0.25),complex(2.0,1.0)]
ga =1.0
ψinvz = ψinv(Finite(-1.0,1.0,-0.5,0.0,0.0,1.0),z)
    global n = length(z)
    global dat = convert(Vector{Float64},real(ψinvz))
    global ept = convert(Vector{Float64},abs(imag(ψinvz)))
    global gaopt = convert(Float64,ga)
    global spg = sinpi(1/2gaopt)
    global cpg = cospi(1/2gaopt)

srand(1234); x = rand(2*n)
eval_h(x,2.0,zeros(2n),zeros(2n),2.0,rand(2n),zeros(n*(2n+1)))
=#
