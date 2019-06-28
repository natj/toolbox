# Mathematical & Numerical functions
############################################
using DSP

export locate,
       interp,
       deriv,
       integ,
       cuminteg,
       smooth!,
       smooth,
       smooth_spline,
       smooth_plaw,
       gauss_laguerre_nw,
       gauss_legendre_nw,
       expi



############################################
# locate
############################################
# Binary search for vectors. Returns matching index.
# If value is outside of the array, return -1
############################################
function locate(arr::AbstractVector, #input vector
                x::Real) #input value

    n = length(arr)
    ascnd = (arr[n] >= arr[1])
    jl = 0
    ju = n+1
    while true
        if ju-jl <= 1
            break
        end
        #jm = int((ju+jl)/2)
        jm = round(Int,((ju+jl)/2))
        if ascnd == (x >= arr[jm])
            jl = jm
        else
            ju = jm
        end
    end

    if x == arr[1]
        return 1
    elseif x == arr[n]
        return n-1
    elseif ascnd && (x > arr[n] || x < arr[1])
        return -1
    elseif !ascnd && (x < arr[n] || x > arr[1])
        return -1
    else
        return jl
    end
end

############################################
# interp
############################################
# Interpolation function for tables
#
# Cubic interpolation & linear extrapolation (lerp)
############################################
function interp(xold::AbstractVector, #grid
                fold::AbstractVector, #values
                xnew::Real; #new value
                method=:cubic) #:cubic, :lin

    n=length(xold)

    #Determining if interpolation or extrapolation
    xmin=xold[1]
    xmax=xold[end]

    #Interpolation
    if xnew >= xmin && xnew <= xmax

        #locate the current index
        ind = locate(xold, xnew)
        #ind = ind < 3 ? 3 : ind+1

        #TODO: method keyword to toggle linear interpolation
        #Linear interpolation
        if method == :lin
	    ind = ind < 2 ? 2 : ind+1
            ind1=ind-1
            ind2=ind
	    x2 = xold[ind2]
	    x1 = xold[ind1]
            fnew=fold[ind2]+(fold[ind2]-fold[ind1])*(xnew-x2)/(x2-x1)
        elseif method == :cubic
            #Cubic interpolation
            ind = ind < 3 ? 3 : ind+1
	    x1=xold[ind-2]
            x2=xold[ind-1]
            x3=xold[ind]
            y1=fold[ind-2]
            y2=fold[ind-1]
            y3=fold[ind]

            fnew = ((x3-xnew)*((x2-x3)*(x2-xnew)*y1+(x3-x1)*(x1-xnew)*y2)+(x1-x2)*(x1-xnew)*(x2-xnew)*y3)/((x1-x2)*(x1-x3)*(x2-x3))
        else
            error("Check method")
        end
        
    #Extrapolation (inner)
    elseif xnew < xmin
        x1=xold[1]
        x2=xold[2]
        fnew = fold[1]+(fold[2]-fold[1])*(xnew-x1)/(x2-x1)

    #Extrapolation (outer)
    elseif xnew > xmax
        x1=xold[n-1]
        x2=xold[n]
        fnew = fold[n]+(fold[n]-fold[n-1])*(xnew-x2)/(x2-x1)

    end

    return fnew
end

#vectorized interp
function interp(xold::AbstractVector,
                fold::AbstractVector,
                xnew::AbstractVector; 
                method=:cubic)

    N = length(xnew)
    fnewarr = zeros(N)
    for i = 1:N
        fnewarr[i] = interp(xold, fold, xnew[i], method=method)
    end

    return fnewarr
end

############################################
# Deriv
############################################
# First derivative of vector
#
# Uses weighted average of two parabolas to compute the derivative
############################################
function deriv(x::AbstractVector, #grid
               f::AbstractVector) #values

    n=length(x)
    @assert n > 1
    @assert length(x) == length(f)

    dfdx=zeros(n)

    dfdx[1]=(f[2]-f[1])/(x[2]-x[1])
    n1=n-1
    dfdx[n]=(f[n]-f[n1])/(x[n]-x[n1])

    if(n==2)
        return dfdx
    end

    s=abs(x[2]-x[1])/(x[2]-x[1])
    for j in 2:n1
        scale=maximum([abs(f[j-1]),abs(f[j]),abs(f[j+1])])/abs(x[j])
        if(scale==0.0)
            scale=1.0
        end

        d1=(f[j+1]-f[j])/(x[j+1]-x[j])/scale
        d=(f[j]-f[j-1])/(x[j]-x[j-1])/scale
        tann1=d1/(s*sqrt(1.0+d1^2.0)+1.0)
        tann=d/(s*sqrt(1.0+d^2.0)+1.0)
        dfdx[j]=(tann1+tann)/(1.0-tann1*tann)*scale
    end

   return dfdx
end

############################################
# Integ
############################################
# Integration of arrays
#
# Uses weighted parabolas
############################################

#Weighted parabolas
function integ(x::AbstractVector,
               f::AbstractVector)

    function parcoe(f, x)
        a = similar(f)
        b = similar(f)
        c = similar(f)

        n = length(f)

        c[1] = 0.0
        b[1] = (f[2]-f[1])/(x[2]-x[1])
        a[1] = f[1]-x[1]*b[1]

        c[n] = 0.0
        b[n] = (f[n]-f[n-1])/(x[n]-x[n-1])
        a[n] = f[n]-x[n]*b[n]
        if n == 2; return a, b, c; end

        for j = 2:(n-1)
            d = (f[j]-f[j-1])/(x[j]-x[j-1])
            c[j] = f[j+1]/((x[j+1]-x[j])*(x[j+1]-x[j-1]))-f[j]/((x[j]-x[j-1])*(x[j+1]-x[j]))+f[j-1]/((x[j]-x[j-1])*(x[j+1]-x[j-1]))
            b[j] = d- (x[j]+x[j-1])*c[j]
            a[j] = f[j-1]-x[j-1]*d+x[j]*x[j-1]*c[j]
        end

        for j = 2:(n-2)
            if c[j] == 0.0; continue; end
            wt = abs(c[j+1])/(abs(c[j+1])+abs(c[j]))
            a[j] = a[j+1]+wt*(a[j]-a[j+1])
            b[j] = b[j+1]+wt*(b[j]-b[j+1])
            c[j] = c[j+1]+wt*(c[j]-c[j+1])
        end

        return a, b, c
    end

    n = length(x)
    if n == 1
        return x .* f
    end
    
    a, b, c = parcoe(f, x)
    fint = zeros(n)

    fint[1] = 0.0
    #finti = (a[1]+(b[1]/2.0+c[1]/3.0*x[2])*x[2])*x[2]-x[1]
    finti = (a[1]+b[1]/2.0*(x[2]+x[1])+c[1]/3.0*((x[2]+x[1])*x[2]+x[1]*x[1]))*(x[2]-x[1])


    #if unstable, use linear interpolation
    #TODO: check when this is needed
    finti = finti < 0.0 ? 0.5*(x[2]-x[1])*(f[1]+f[2]) : finti
    fint[2] = finti

    if n == 2; return fint; end

    for i = 2:(n-1)
        finti = (a[i]+b[i]/2.0*(x[i+1]+x[i])+c[i]/3.0*((x[i+1]+x[i])*x[i+1]+x[i]*x[i]))*(x[i+1]-x[i])

        #if unstable, use linear interpolation
        #TODO: check when this is needed
        finti = finti < 0.0 ? 0.5*(x[i+1]-x[i])*(f[i]+f[i+1]) : finti

        fint[i+1] = finti
    end
    return fint
end

#Cumulative integral
function cuminteg(x::AbstractVector,
                  f::AbstractVector;
                  extrapolate=:zlin) #:zlin, lin, :zplaw, :none

    fint = integ(x, f)
    
    #extrapolate from 0 to x1
    if extrapolate == :zlin #linear [0,x1]
        fint[1] = 0.5*x[1]*f[1]
#    elseif extrapolate_zero == :quad #quadratic
#        c = 0.0
#        b = (f[2]-f[1])/(x[2]-x[1])
#        a = f[1]-x[1]*b[1]
#        fint[1] = (a+(b/2.0+(c/3.0*x[1])*x[1])*x[1])
    elseif extrapolate == :lin #linear
        k = (f[2]-f[1])/(x[2]-x[1])
        b = f[1] - k*x[1]
        if b >= 0.0
            fint[1] = 0.5*k*x[1]^2+b*x[1]
        else
            x0 = -b/k
            fint[1] = 0.5*k*(x[1]^2.0-x0^2.0)+b*(x[1]-x0)
        end
    elseif extrapolate == :zplaw #powerlaw [0,x1]
        a=log(f[2]/f[1])/log(x[2]/x[1])
        c=f[1]/(x[1]^a)
        fint[1]=c*(x[1]^(a+1))/(a+1)
    elseif extrapolate == :none #no extrapolation
        fint[1] = 0.0
    else
        error("Unrecognized extrapolation method $(string(extrapolate))")
    end

    cfint = cumsum(fint)
    return cfint
end

#Simpsons rule (cumulative)
#may or may not work
function integ_Simpson(x,f)

    n = length(f)
    fint = Array(Float64,n)

    #Linear extrapolation between 0-x1
    fint[1] = (x[1]/(2.0*(x[1]-x[2])))*(-2.0*x[2]*f[1]+x[1]*(f[1]+f[2]))

    #Quadratic extrapolation between 0-x1
#    fint[1]=(x[1]*(6*x[2]*(x[2]-x[3])*x[3]*f[1]+x[1]^3.0*(f[2]-f[3])+3*x[1]*(x[3]^2*(f[1]+f[2])-x[2]^2.0*(f[1]+f[3]))+2*x[1]^2.0*(-x[3]*(f[1]+2*f[2])+x[2]*(f[1]+2*f[3]))))/(6*(x[1]-x[2])*(x[1]-x[3])*(x[2]-x[3]))
#    fint[1]=f[1]

    #Simpsons rule for everything else
    for i in 1:n-2

        d31=x[i+2]-x[i]
        d32=x[i+2]-x[i+1]
        d21=x[i+1]-x[i]

        x3=x[i+2]
        x2=x[i+1]
        x1=x[i]
        y3=f[i+2]
        y2=f[i+1]
        y1=f[i]

        integ=(d21/(6.0*d31*d32))*((3.0*x3^2.0)*(y1+y2)-2.0*x2*x3*(2.0*y1+y2)+(x2^2.0)*(y1-y3)+(x1^2.0)*(y2-y3)+2*x1*(x2*(y1+y2+y3)-x3*(y1+2.0*y2)))

#        println("i=$i $integ")

        fint[i+1]=fint[i]+abs(integ)

        if i==(n-2)
            integ=-(d32/(6.0*d21*d31))*((x3^2.0)*(y1-y2)+(x2^2.0)*(y1-y3)-(3.0*x1^2.0)*(y2+y3)+2.0*x1*x3*(2.0*y2+y3)-2.0*x2*(x3*(y1+y2+y3)-x1*(y2+2.0*y3)))
            fint[i+2]=fint[i+1]+integ
        end
    end

return fint
end

############################################
# smooth
############################################
#
# Smooth vector with Gaussian kernel
#
############################################
function smooth!(arr::AbstractVector, #input vector
                 N::Int=1; #number of smoothings
                 offs::Int=3) #offset of smoothing

    l = length(arr)
    offs = 3

    #Gaussian 7x1 kernel
    #XXX: accept other kernels too
    kernel = Float64[0.00443185,
                            0.05399100,
                            0.24197100,
                            0.39894200,
                            0.24197100,
                            0.05399100,
                            0.00443185]

    #XXX is this correct?
    logkernel = Float64[
                            0.00443185,
                            0.07342030,
                            0.27844800,
                            0.39894200,
                            0.16782600,
                            0.02927920,
                            0.00443185
                            ]


    #arr2 = zeros(l+2offs)
    for s = 1:N

        #Expand original array to avoid boundaries
        #XXX: extrapolate instead
        #arr2[1:offs] = ones(offs)*arr[1]
        #arr2[(offs+1):(end-offs)] = arr[:]
        #arr2[(end-offs+1):end] = ones(offs)*arr[l]
        #sarr = conv(arr2, kernel)
        #arr[:] = sarr[(1+offs+3):(end-offs-3)]

        #smooth only middle parts
        sarr = conv(arr, kernel)
        arr[(offs+1):(l-offs)] = sarr[(1+offs+3):(end-offs-3)]

    end

    return arr
end
smooth(arr, N=1; kvs...) = smooth!(deepcopy(arr), N; kvs...)

############################################
# smooth_spline
############################################
#
# Smooth vector with B-splines according to DeBoor's
# algorithm.
#
# Adapted from Fred Frigo's (Dec 8, 2001) Matlab source file
# that originates from http://www.psc.edu/~burkardt/src/splpak/splpak.f90
############################################
function smooth_spline(y::AbstractVector, dx::AbstractVector, smooth_factor::Real)

npoint = length(y)
p = smooth_factor

a = zeros(npoint, 4)
v = zeros(npoint, 7)

x = linspace(0.0, (npoint-1.0)/npoint, npoint)

# setupq

  v[1,4] = x[2]-x[1]

  for i = 2:npoint-1
    v[i,4] = x[i+1]-x[i]
    v[i,1] = dx[i-1]/v[i-1,4]
    v[i,2] = ((-1.0 .*dx[i])/v[i,4]) - (dx[i]/v[i-1,4])
    v[i,3] = dx[i+1]/v[i,4]
  end


  v[npoint,1] = 0.0
  for i = 2:npoint-1
    v[i,5] = (v[i,1]*v[i,1]) + (v[i,2]*v[i,2]) + (v[i,3]*v[i,3])
  end

  for i = 3:npoint-1
    v[i-1,6] = (v[i-1,2]*v[i,1]) + (v[i-1,3]*v[i,2])
  end

  v[npoint-1,6] = 0.0

  for i = 4: npoint-1
    v[i-2,7] = v[i-2,3]*v[i,1]
  end

  v[npoint-2,7] = 0.0
  v[npoint-1,7] = 0.0

#  Construct  q-transp. * y  in  qty.

  prev = (y[2]-y[1])/v[1,4]
  for i= 2:npoint-1
    diff = (y[i+1]-y[i])/v[i,4]
#    %qty(i) = diff-prev
    a[i,4] = diff - prev
    prev = diff
  end

# end setupq

#chol1d

#  Construct 6*(1-p)*q-transp .*(d**2)*q + p*r

  six1mp = 6.0 .*(1.0-p)
  twop = 2.0 .*p

  for i = 2: npoint-1
    v[i,1] = (six1mp .*v[i,5]) + (twop .*(v[i-1,4]) + v[i,4])
    v[i,2] = (six1mp .*v[i,6]) +( p .*v[i,4])
    v[i,3] = six1mp .*v[i,7]
  end

  if  npoint < 4
    u[1] = 0.0
    u[2] = a[2,4]/v[2,1]
    u[3] = 0.0

#  Factorization

  else
    for i = 2: npoint-2
      ratio = v[i,2]/v[i,1]
      v[i+1,1] = v[i+1,1]-(ratio .*v[i,2])
      v[i+1,2] = v[i+1,2]-(ratio .*v[i,3])
      v[i,2] = ratio
      ratio = v[i,3]./v[i,1]
      v[i+2,1] = v[i+2,1]-(ratio .*v[i,3])
      v[i,3] = ratio
    end

#  Forward substitution

    a[1,3] = 0.0
    v[1,3] = 0.0
    a[2,3] = a[2,4]
    for i = 2: npoint-2
      a[i+1,3] = a[i+1,4] - (v[i,2]*a[i,3]) - (v[i-1,3]*a[i-1,3])
    end

#%  Back substitution.

    a[npoint,3] = 0.0
    a[npoint-1,3] = a[npoint-1,3] / v[npoint-1,1]

    for i = npoint-2:-1:2
      a[i,3] = (a[i,3]/v[i,1]) - (a[i+1,3]*v[i,2]) - (a[i+2,3]*v[i,3])
    end

  end
#
#  Construct Q*U.
#
  prev = 0.0
  for i = 2: npoint
    a[i,1] = (a[i,3]-a[i-1,3])/v[i-1,4]
    a[i-1,1] = a[i,1]-prev
    prev = a[i,1]
  end

  a[npoint,1] = -1.0 .*a[npoint,1]

#end chol1d

  spline_sig = zeros(npoint)
  for i = 1: npoint
    spline_sig[i] = y[i]-(6.0 .*(1.0-p) .*dx[i] .*dx[i] .*a[i,1])
  end

    return spline_sig
end

#no weights
smooth_spline(y::AbstractVector, smooth_factor::Real) = smooth_spline(y, ones(length(y)), smooth_factor)


############################################
# smooth_plaw
############################################
#
# Smooth vector that has powerlaw behaviour by
# smoothing the ratio to some reference vector
# with either gaussian kernel or splines.
#
############################################
function smooth_plaw(x::AbstractVector,
                     ref::AbstractVector, #reference vector
                     N::Real=1; #number of smoothings OR smoothing factor
                     offs::Int=3, #offset of smoothing
                     method=:kernel) #:kernel, :spline

    @assert length(x) == length(ref)

    dr = x./ref

    if method == :kernel
        sdr = smooth(dr, int(N), offs=offs)
    elseif method == :spline
        sdr = smooth_spline(dr, N)
    else
        error("unrecognized method $(string(method))")
    end

    sx = ref .* sdr

    return sx
end
############################################
# gauss_laquerre_nw
############################################
#
# Gauss-Laguerra quadrature
#
# Calculates weights and nodes of Gaussian quadrature
# in interval [0,infty) with weight exp(-x)
#
# output: r (vector of nodes)
#         w (vector of weights)
############################################
function gauss_laguerre_nw(num::Real) #number of points
    #J=diagm([1:2:(2*num-1)])+diagm([1:num-1],1)+diagm([1:num-1],-1)
    J=diagm(0 => [1:2:(2*num-1)])+diagm(1 => [1:num-1])+diagm(-1 => [1:num-1])
    #J=spdiagm(sparsevec)+diagm(Pair{1,[1:num-1]})+diagm(Pair{-1,[1:num-1]})
    R,W=eig(J)
    r=sort(R)
    w=(W[1,:].^2.0)

    return r,w
end

############################################
# gauss_legendre_nw
############################################
#
# Gauss-Legendre quadrature
#
# calculates weights and nodes of Gaussian quadrature
# in interval [-1,1] with weight 1
#
# output: r (vector of nodes)
#         w (vector of weights)
############################################
function gauss_legendre_nw(num::Real)

    nlist=[1:num-1]
    J=diagm(nlist./sqrt(4 .*nlist.^2.0-1),1)+diagm(nlist./sqrt(4 .*nlist.^2.0-1),-1)
    R,W=eig(J)
    r=sort(R)
    w=2.0*(W[1,:].^2.0)

    return r,w
end

############################################
# expi
############################################
#
# Exponential integral for positive agruments
# after Cody and Thacher, Math. of Comp. 22,641 (1968)
#
############################################
function expi(n::Real,
              x::Real)

    @assert x >= 0.0

    #Table of predefined values
    X1=-1.0e20

    A0=-44178.5471728217
    A1=57721.7247139444
    A2=9938.31388962037
    A3=1842.11088668000
    A4=101.093806161906
    A5=5.03416184097568
    B0=76537.3323337614
    B1=32597.1881290275
    B2=6106.10794245759
    B3=635.419418378382
    B4=37.2298352833327

    C0=4.65627107975096e-07
    C1=0.999979577051595
    C2=9.04161556946329
    C3=24.3784088791317
    C4=23.0192559391333
    C5=6.90522522784444
    C6=0.430967839469389
    D1=10.0411643829054
    D2=32.4264210695138
    D3=41.2807841891424
    D4=20.4494785013794
    D5=3.31909213593302
    D6=0.103400130404874

    E0=-0.999999999998447
    E1=-26.6271060431811
    E2=-241.055827097015
    E3=-895.927957772937
    E4=-1298.85688746484
    E5=-545.374158883133
    E6=-5.66575206533869
    F1=28.6271060422192
    F2=292.310039388533
    F3=1332.78537748257
    F4=2777.61949509163
    F5=2404.01713225909
    F6=631.657483280800

    if x==X1
#        expint=ex1
        expint=0.0
        if n==1
            return expint
        end
        n1=n-1
        for i in 1:n1
            expint=(ex-x*expint)/(1.0*i)
        end
        return expint
    end

    ex=exp(-x)
    x1=x
    if x>4.0
        ex1=(ex+ex*(E0+(E1+(E2+(E3+(E4+(E5+E6/x)/x)/x)/x)/x)/x)/(x+F1+(F2+(F3+(F4+(F5+F6/x)/x)/x)/x)/x))/x
    elseif x>1.0
        ex1=ex*(C6+(C5+(C4+(C3+(C2+(C1+C0*x)*x)*x)*x)*x)*x)/(D6+(D5+(D4+(D3+(D2+(D1+x)*x)*x)*x)*x)*x)
    elseif x>0.0
        ex1=(A0+(A1+(A2+(A3+(A4+A5*x)*x)*x)*x)*x)/(B0+(B1+(B2+(B3+(B4+x)*x)*x)*x)*x)-log(x)
    else
        ex1=0.0
    end

    expint=ex1
    if n==1
        return expint
    end

    for i in 1:n-1
        expint=(ex-x*expint)/(1.0*i)
    end

    return expint
end
