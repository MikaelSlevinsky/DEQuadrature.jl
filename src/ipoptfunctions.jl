#
# These functions are required in the use of Ipopt.
#
function eval_f(x)
    temp1 = 0.0
    temp2 = 0.0
    @inbounds for k=1:n
        temp3=0.0
        @inbounds for j=1:n
            temp3+=x[n+j]*complex(x[k],pi/2gaopt)^(j-1)
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
    @inbounds for k=1:n
        temp1=0.0
        @inbounds for j=1:n
            temp1+=x[n+j]*complex(x[k],pi/2gaopt)^(j-1)
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
    @inbounds for r=1:n
        temp1=0.0
        temp2=0.0
        temp4=0.0
        temp5=0.0
        @inbounds for k=1:n
            temp3=0.0
            @inbounds for j=1:n
                temp3+=x[n+j]*complex(x[k],pi/2gaopt)^(j-1)
            end
            temp1+=ept[k]-imag(temp3)
            temp2+=cosh(x[k])*spg
            temp4+=x[n+k]*(k-1)*complex(x[r],pi/2gaopt)^(k-2)
            temp5+=complex(x[k],pi/2gaopt)^(r-1)
        end
        grad_f[r] = -(temp2*imag(temp4) + sinh(x[r])*spg*temp1)/temp2^2
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
        @inbounds for k=1:n
            temp1=0.0
            @inbounds for r=1:n
                values[2n*(k-1)+r]=grad_f[r]*sinh(x[k])*cpg
                values[2n*(k-1)+n+r]=grad_f[n+r]*sinh(x[k])*cpg

                values[2n^2+2n*(k-1)+r]=grad_f[r]*cosh(x[k])*spg
                values[2n^2+2n*(k-1)+n+r]=grad_f[n+r]*cosh(x[k])*spg

                values[2n*(k-1)+n+r] += real(complex(x[k],pi/2gaopt)^(r-1))
                values[2n^2+2n*(k-1)+n+r] += imag(complex(x[k],pi/2gaopt)^(r-1))

                temp1+=x[n+r]*(r-1)*complex(x[k],pi/2gaopt)^(r-2)
            end

            values[2n*(k-1)+k] += f*cosh(x[k])*cpg
            values[2n^2+2n*(k-1)+k] += f*sinh(x[k])*spg

            values[2n*(k-1)+k] += real(temp1)
            values[2n^2+2n*(k-1)+k] += imag(temp1)
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
    # Symmetric matrix, fill the lower left triangle only
    if mode == :Structure
        # Again, only lower left triangle
        # Objective
        idx = 1
        for row = 1:2n
            for col = 1:row
                rows[idx] = row
                cols[idx] = col
                idx += 1
            end
        end
    else
        values[:] = zeros(n*(2n+1))
        constraints = zeros(n*(2n+1))
        grad_f = zeros(2n)
        f = eval_f(x)
        eval_grad_f(x,grad_f)  
        @inbounds for r=1:2n
            @inbounds for p=1:r
            temp1=0.0
            temp2=0.0
            temp4=0.0
            temp5=0.0
            temp6=0.0
            temp7=0.0
            @inbounds   for k=1:n
                            temp3=0.0
                 @inbounds  for j=1:n
                                temp3+=x[n+j]*complex(x[k],pi/2gaopt)^(j-1)
                            end # for j
                            temp1+=ept[k]-imag(temp3)
                            temp2+=cosh(x[k])*spg
                            temp4+=x[n+k]*(k-1)*complex(x[r],pi/2gaopt)^(k-2)
                            temp5+=complex(x[k],pi/2gaopt)^(r-1)
                            temp6+=x[n+k]*(k-1)*(k-2)*complex(x[r],pi/2gaopt)^(k-3)
                            temp7+=x[n+k]*(k-1)*complex(x[p],pi/2gaopt)^(k-2)                
                        end # for k
            values[int(r*(r-1)/2+p)] =  (r<=n? ((imag(temp4)*sinh(x[p])*spg)/temp2^2 + sinh(x[r])*spg*(temp2*imag(temp7)+2*temp1*sinh(x[p])*spg)/temp2^3 - (r == p? imag(temp6)/temp2 + cosh(x[r])*spg*(temp1/temp2^2) : 0.0 )) :  (p>n? 0.0: (-temp2*imag((r-1)*complex(x[p],pi/2gaopt)^(r-2)) + sinh(x[p])*spg*imag(temp5))/temp2^2) )    
                for kk =1:2n        
                    if kk<=n
                        constraints[int(r*(r-1)/2+p)] += lambda[kk]*(r<=n? values[int(r*(r-1)/2+p)]*sinh(x[kk])*cpg + (kk==p?  grad_f[r]*cosh(x[p])*cpg : 0.0) + (kk==r?  grad_f[p]*cosh(x[r])*cpg : 0.0) + (kk==r==p?  f*sinh(x[p])*cpg+real(temp6) : 0.0)  : (p>n? 0.0 : values[int((n+r)*(n+r-1)/2+p)]*sinh(x[kk])*cpg + (k==p?  grad_f[r]*cosh(x[p])*cpg + real((r-1)*complex(x[p],pi/2gaopt)^(r-2)) : 0.0) ) )
                    elseif kk>n && kk<2n
                        constraints[int(r*(r-1)/2+p)] += lambda[kk]*(r<=n? values[int(r*(r-1)/2+p)]*cosh(x[kk])*spg + (kk==p?  grad_f[r]*sinh(x[p])*spg : 0.0) + (kk==r?  grad_f[p]*sinh(x[r])*spg : 0.0) + (kk==r==p?  f*cosh(x[p])*spg+imag(temp6) : 0.0)  : (p>n? 0.0 : values[int((n+r)*(n+r-1)/2+p)]*cosh(x[kk])*spg + (k==p?  grad_f[r]*sinh(x[p])*spg + imag((r-1)*complex(x[p],pi/2gaopt)^(r-2)) : 0.0) ) )
                    else
                        constraints[int(r*(r-1)/2+p)] += 0.0    
                    end
                end #for kk
                values[int(r*(r-1)/2+p)]  = obj_factor*values[int(r*(r-1)/2+p)] + constraints[int(r*(r-1)/2+p)]
            end  # for p
        end # for r
    end # if loop
end # function
