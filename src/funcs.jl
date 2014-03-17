# Mathematical & Numerical functions
############################################

export locate,
       interp,
       deriv,
       integ,
       smooth!,
       smooth,
       gauss_laguerre_nw,
       gauss_legendre_nw,
       expi



############################################
# locate
############################################
# Binary search for vectors. Returns matching index
#
# input: arr, x
#
# output: mid (index)
############################################
function locate(arr, x)
    n = length(arr)
    ascnd = (arr[n] >= arr[1])
    jl = 0
    ju = n+1
    while true
        if ju-jl <= 1
            break
        end
        jm = int((ju+jl)/2)
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
# Takes only one x0 grid point and returns f(x0)
#
# input: xold (old grid), fold (old values), xnew (new grid)
#
# output: fnew (new values interpolated to xnew grid)
############################################
function interp(xold,fold,xnew)
    n=length(xold)

    #Determining if interpolation or extrapolation
    xmin=xold[1]
    xmax=xold[end]

    #Interpolation
    if xnew >= xmin && xnew <= xmax

        x1=xold[n-1]
        x2=xold[n]
        ind=n

        while !(xnew >= x1 && xnew <= x2) && ind>3
            ind=ind-1
            x1=xold[ind-1]
            x2=xold[ind]
        end

        #Linear interpolation
#        ind=ind-1
#        ind2=ind
#        fnew=fold[ind2]+(fold[ind2]-fold[ind])*(xnew-x2)/(x2-x1)

        #Cubic interpolation
        x1=xold[ind-2]
        x2=xold[ind-1]
        x3=xold[ind]
        y1=fold[ind-2]
        y2=fold[ind-1]
        y3=fold[ind]

        fnew = ((x3-xnew)*((x2-x3)*(x2-xnew)*y1+(x3-x1)*(x1-xnew)*y2)+(x1-x2)*(x1-xnew)*(x2-xnew)*y3)/((x1-x2)*(x1-x3)*(x2-x3))

    #Extrapolation (inner)
    elseif xnew < xmin
        x1=xold[1]
        x2=xold[2]
        fnew = fold[1]-(fold[2]-fold[1])*(xnew-x1)/(x2-x1)

    #Extrapolation (outer)
    elseif xnew > xmax
        x1=xold[n-1]
        x2=xold[n]
        fnew = fold[n]+(fold[n]-fold[n-1])*(xnew-x2)/(x2-x1)

    end

    return fnew
end


############################################
# Deriv
############################################
# First derivative of vector
#
# Uses weighted average of two parabolas to compute the derivative
#
# input: x,f,n
#
# output: fint
############################################

function deriv(x,f)

    n=length(x)
    if(n==1)
        error("Can not derivate a scalar")
    end

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
# Integ2
############################################
# Basic cumulative integration of arrays
# starting from 0.0
#
# Uses Simpsons rule
#
# input: x,f
#
# output: fint
############################################
#

#Weighted parabolas
function integ(x, f)

    function parcoe(f,x)
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
    a, b, c = parcoe(f, x)
    fint = similar(f)

    #extrapolate from 0 to x1
#    fint[1] = (a[1]+(b[1]/2.0+(c[1]/3.0*x[1])*x[1])*x[1]) #quadratic
    fint[1] = 0.5*x[1]*f[1] #linear

    #i=2
    finti = (a[1]+(b[1]/2.0+c[1]/3.0*x[2])*x[2])*x[2]

    #if unstable, use linear interpolation
    if finti < 0.0
        finti = 0.5*(x[2]-x[1])*(f[1]+f[2])
    end
    fint[2] = fint[1]+finti

    if n == 2; return fint; end

    for i = 2:(n-1)
        finti = (a[i]+b[i]/2.0*(x[i+1]+x[i])+c[i]/3.0*((x[i+1]+x[i])*x[i+1]+x[i]*x[i]))*(x[i+1]-x[i])

        #if unstable, use linear interpolation
        if finti < 0.0
            finti = 0.5*(x[i+1]-x[i])*(f[i]+f[i+1])
        end

        fint[i+1] = fint[i] + finti
    end
    return fint
end

#Simpsons rule
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
# input: arr,
#        N
#
# output: arr
############################################
#
function smooth!(arr, N::Int=1)

    l = length(arr)
    offs = 6

    #Gaussian 7x1 kernel
    #XXX: accept other kernels too
    const kernel = [0.00443185,
                    0.05399100,
                    0.24197100,
                    0.39894200,
                    0.24197100,
                    0.05399100,
                    0.00443185]

    arr2 = zeros(l+2offs)

    for s = 1:N

        #Expand original array to avoid boundaries
        #XXX: extrapolate instead
        arr2[1:offs] = ones(offs)*arr[1]
        arr2[(offs+1):(end-offs)] = arr[:]
        arr2[(end-offs+1):end] = ones(offs)*arr[l]

        sarr = conv(arr2, kernel)

        arr[:] = sarr[(1+offs+3):(end-offs-3)]
    end

    return arr
end
smooth(arr, N=1) = smooth!(deepcopy(arr), N)

############################################
# gauss_laquerre_nw
############################################
#
# Gauss-Laguerra quadrature
#
# Calculates weights and nodes of Gaussian quadrature
# in interval [0,infty) with weight exp(-x)
#
# input: num (number of points)
#
# output: r (vector of nodes)
#         w (vector of weights)
############################################
#
function gauss_laguerre_nw(num)

    J=diagm([1:2:(2*num-1)])+diagm([1:num-1],1)+diagm([1:num-1],-1)
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
# input: num (number of points)
#
# output: r (vector of nodes)
#         w (vector of weights)
############################################
#
function gauss_legendre_nw(num)

    nlist=[1:num-1]
    J=diagm(nlist./sqrt(4.*nlist.^2.0-1),1)+diagm(nlist./sqrt(4.*nlist.^2.0-1),-1)
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
# input: n
#
# output: expint
############################################
#
function expi(n,x)

    @assert x >= 0.0


    const X1=-1.0e20

    const A0=-44178.5471728217
    const A1=57721.7247139444
    const A2=9938.31388962037
    const A3=1842.11088668000
    const A4=101.093806161906
    const A5=5.03416184097568
    const B0=76537.3323337614
    const B1=32597.1881290275
    const B2=6106.10794245759
    const B3=635.419418378382
    const B4=37.2298352833327

    const C0=4.65627107975096e-07
    const C1=0.999979577051595
    const C2=9.04161556946329
    const C3=24.3784088791317
    const C4=23.0192559391333
    const C5=6.90522522784444
    const C6=0.430967839469389
    const D1=10.0411643829054
    const D2=32.4264210695138
    const D3=41.2807841891424
    const D4=20.4494785013794
    const D5=3.31909213593302
    const D6=0.103400130404874

    const E0=-0.999999999998447
    const E1=-26.6271060431811
    const E2=-241.055827097015
    const E3=-895.927957772937
    const E4=-1298.85688746484
    const E5=-545.374158883133
    const E6=-5.66575206533869
    const F1=28.6271060422192
    const F2=292.310039388533
    const F3=1332.78537748257
    const F4=2777.61949509163
    const F5=2404.01713225909
    const F6=631.657483280800

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
