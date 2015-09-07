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
        for r=1:n
            temp7=0.0 # depends on p, right? (I'm not sure.)

            for p=1:r
                temp1=0.0
                temp2=0.0
                temp4=0.0
                temp5=0.0
                temp6=0.0
                for k=1:n
                    temp3=0.0
                    for j=1:n
                        temp3+=x[n+j]*xpg[k]^(j-1)
                    end # for j
                    temp1+=ept[k]-imag(temp3)
                    temp2+=cosh(x[k])*spg
                    temp4+=x[n+k]*(k-1)*xpg[r]^(k-2)
                    temp5+=xpg[k]^(r-1)
                    temp6+=x[n+k]*(k-1)*(k-2)*xpg[r]^(k-3)
                    temp7+=x[n+k]*(k-1)*xpg[p]^(k-2)
                end # for k
                #values[int(r*(r-1)/2+p)] =  (r<=n? ((imag(temp4)*sinh(x[p])*spg)/temp2^2 + sinh(x[r])*spg*(temp2*imag(temp7)+2*temp1*sinh(x[p])*spg)/temp2^3 - (r == p? imag(temp6)/temp2 + cosh(x[r])*spg*(temp1/temp2^2) : 0.0 )) :  (p>n? 0.0: (-temp2*imag((r-1)*xpg[p]^(r-2)) + sinh(x[p])*spg*imag(temp5))/temp2^2) )

                # ∂^2 f / ∂x_r ∂ x_p
                #values[int(r*(r-1)/2+p)] =

                # ∂^2 f / ∂u_r ∂ x_p = ∂^2 f / ∂x_{n+r} ∂ x_p
                #values[int((r+n)*(r+n-1)/2+p)] =
            end # for p

            # ∂^2 f / ∂x_r^2
            #values[int(r*(r-1)/2+r)] =
        end # for r

        for k =1:n
            for r=1:n
                for p=1:r

                    # requires a loop for temp6

                    #constraints[int(r*(r-1)/2+p)] += lambda[k]*(r<=n? values[int(r*(r-1)/2+p)]*sinh(x[kk])*cpg + (kk==p?  grad_f[r]*cosh(x[p])*cpg : 0.0) + (kk==r?  grad_f[p]*cosh(x[r])*cpg : 0.0) + (kk==r==p?  f*sinh(x[p])*cpg+real(temp6) : 0.0)  : (p>n? 0.0 : values[int((n+r)*(n+r-1)/2+p)]*sinh(x[kk])*cpg + (kk==p?  grad_f[r]*cosh(x[p])*cpg + real((r-1)*xpg[p]^(r-2)) : 0.0) ) )

                    # ∂^2 g_k / ∂x_r ∂x_p
                    #constraints[int(r*(r-1)/2+p)] += lambda[k]*

                    # ∂^2 g_{n+k} / ∂x_r ∂x_p
                    if k < n
                        #constraints[int(r*(r-1)/2+p)] += lambda[n+k]*
                    else
                        #constraints[int(r*(r-1)/2+p)] += lambda[2n]*0.0
                    end

                    # ∂^2 g_k / ∂u_r ∂x_p = ∂^2 g_k / ∂x_{n+r} ∂x_p
                    #constraints[int((n+r)*(n+r-1)/2+p)] += lambda[k]*

                    # ∂^2 g_{n+k} / ∂x_{n+r} ∂x_p
                    if k < n
                        #constraints[int((n+r)*(n+r-1)/2+p)] += lambda[n+k]*
                    else
                        #constraints[int((n+r)*(n+r-1)/2+p)] += lambda[2n]*0.0
                    end
                end #for p
                # ∂^2 g_k / ∂x_r ∂x_k
                #constraints[int(r*(r-1)/2+k)] += lambda[k]*

                # ∂^2 g_{n+k} / ∂x_r ∂x_k
                if k < n
                    #constraints[int(r*(r-1)/2+k)] += lambda[n+k]*
                else
                    #constraints[int(r*(r-1)/2+k)] += lambda[2n]*0.0
                end

                # ∂^2 g_k / ∂x_r ∂u_k = ∂^2 g_k / ∂x_r ∂x_{n+k}
                #constraints[int(r*(r-1)/2+n+k)] += lambda[k]*

                # ∂^2 g_{n+k} / ∂x_r ∂u_k = ∂^2 g_{n+k} / ∂x_r ∂u_{n+k}
                if k < n
                    #constraints[int(r*(r-1)/2+n+k)] += lambda[n+k]*
                else
                    #constraints[int(r*(r-1)/2+n+k)] += lambda[2n]*0.0
                end
            end # for r
            # ∂^2 g_k / ∂x_k^2
            #constraints[int(k*(k-1)/2+k)] += lambda[k]*
            # ∂^2 g_{n+k} / ∂x_k^2
            #constraints[int(k*(k-1)/2+k)] += lambda[n+k]*
        end # for k

        values[:] = obj_factor*values[:] + constraints[:]
    end # if loop
end # function
