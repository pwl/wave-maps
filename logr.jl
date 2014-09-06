using DASSL
include("discretization.jl")

function initial(u0::Function,r)
    npts = length(r)
    v1=map(u0,r)
    v2=zero(r)
    y0=reshape([v1 v2], 2*npts)
    return y0
end

function newF(x0,x1,d)
    function rhs(t,y,dy)
        res = zero(y)
        npts=int(length(y)/2)
        dv=reshape(dy,npts,2)
        v=reshape( y,npts,2)
        x=linspace(x0,x1,npts)
        dx=x[2]-x[1]

        ord = 2
         vxi= vx(v[:,1],dx,order=ord)
        vxxi=vxx(v[:,1],dx,order=ord)

        resv = reshape(res,npts,2)
        resv[:,1]=v[:,2]-dv[:,1]
        resv[:,2]=vxxi+(d-2)*vxi-(d-1)/2*sin(2*v[:,1])-dv[:,2].*exp(2*x)
        return res
    end
    return rhs
end

function vx(v,dx;order=2)
    vxi = zero(v)
    npts = length(v)
    if order == 2
        for i = 1:npts
            if i == 1
                # vxi[i]= -(3*v[i]-4*v[i+1]+v[i+2]) / dx / 2
                vxi[i]=(exp(-dx)*(exp(dx)*(exp(dx)-2)*v[1]+v[2]))/(exp(dx)-1)
            elseif i == npts
                vxi[i]= (3*v[i]-4*v[i-1]+v[i-2]) / dx / 2
            else
                # vxi[i]=(v[i+1]-v[i-1]) / dx / 2
                vxi[i]=1/2*(-1 + coth(dx))*(-v[i]+exp(2*dx)*(-v[-1 + i] + v[i]) + v[1 + i])
            end
        end
    elseif order == 4
        for i = 1:npts
            if i == 1
                vxi[i]= -csch(2*dx)*(1+cosh(dx)+cosh(2*dx)-(1+4*cosh(dx))*sinh(dx))*v[1]+(2*exp(-2*dx)*v[2])/(exp(dx)-1)-exp(-5*dx)*csch(2*dx)*v[3]
            elseif i == 2
                vxi[i]=(exp(-dx)*(-4+exp(dx)*(-5-3*exp(dx)+exp(3*dx)))*v[1])/(2*(sinh(dx)+sinh(2*dx)))+(exp(-dx)*(1+2*cosh(dx))*v[2])/(exp(dx)-1)-1/2*exp(-5*dx)*(1+2*cosh(dx))*(1+coth(dx))*v[3]+(exp(-5*dx)*v[4])/(-1+exp(3*dx))
            elseif i == npts-1
                vxi[i]=( -v[i-3]+ 6*v[i-2]-18*v[i-1]+10*v[i]   +3*v[i+1])/ dx / 12
            elseif i == npts
                vxi[i]=(3*v[i-4]-16*v[i-3]+36*v[i-2]-48*v[i-1]+25*v[i]  )/ dx / 12
            else
                vxi[i]=(  v[i-2]- 8*v[i-1]           +8*v[i+1]   -v[i+2])/ dx / 12
            end
        end
    end
    return vxi
end

function vxx(v,dx;order = 2)
    vxxi = zero(v)
    npts = length(v)
    if order == 2
        for i = 1:npts
            if i == 1
                # vxxi[i]=(v[i]-2v[i+1]+v[i+2])/dx^2
                vxxi[i]=(1-3/(exp(dx)-1))*v[1]+(3*exp(-dx)*v[2])/(exp(dx)-1)
            elseif i == npts
                vxxi[i]=(v[i-2]-2v[i-1]+v[i])/dx^2
            else
                # vxxi[i]=(v[i-1]-2v[i]+v[1+i])/dx^2
                vxxi[i]=(-exp(2*dx)*(-3 + exp(dx))*v[-1 + i] +  2*exp(dx)*(1 + exp(dx))*(-2 + cosh(dx))*v[i] + (-1 + 3*exp(dx))*v[1 + i])/((-1 + exp(dx))^2*(1 + exp(dx)))
            end
        end
    elseif order == 4
        for i = 1:npts
            if i == 1
                vxxi[i]=((3+exp(dx)*(3+exp(dx)*(-6+exp(dx)*(-5+exp(dx))*(-2+exp(dx)))))v[1])/((-1+exp(dx))^2*(1+exp(dx))*(1+exp(2dx)))+(2*(-5+3*exp(dx))*v[2])/((-1+exp(dx))^2*(1+exp(dx)))+1/4*exp(-4*dx)*csch(dx/2)*(csch(dx/2)+2*sech(dx/2))*sech(dx)*v[3]
            elseif i == 2
                vxxi[i]=v[1]-exp(dx)*v[1]+2*v[2]+exp(-4dx)*v[3]+(v[1]-2v[2]+v[3])/(-1+exp(dx))^2-exp(-dx)*(2v[2]+v[3]-4v[4])+(5*(-1+exp(dx))*(v[1]-v[4]))/(3*(1+exp(dx)+exp(2*dx)))+(2*(v[2]-v[4]))/(1+exp(dx))+(v[1]-3v[2]+3v[3]-v[4])/(3*(-1+exp(dx)))-exp(-5*dx)*v[4]+exp(-3*dx)*(v[3]+4v[4])-exp(-2*dx)*(v[2]+2*v[3]+5*v[4])
            elseif i == npts-1
                vxxi[i]=-(v[i-3]-4v[i-2]-6v[i-1]+20v[i]-11v[1+i])/(12*dx^2)
            elseif i == npts
                vxxi[i]=(11v[i-4]-56v[i-3]+114v[i-2]-104v[i-1]+35v[i])/(12*dx^2)
            else
                vxxi[i]=-(v[i-2]-16v[i-1]+30v[i]-16v[1+i]+v[2+i])/(12*dx^2)
            end
        end
    end
    return vxxi
end

function extract(t,y,dy)
    npts = int(length(y[1])/2)
     x = linspace(x0,x1,npts)
     r = exp(x)
     v = map(yi->reshape(yi,npts,2)[:,1],y)
    du = map(vi->exp(-x).*vx(vi),v)
    return r,v,du             # v=u
end

function wavelogsolve(uinit :: Function;
                      npts = 100,
                      x0 = -5,
                      x0min = -10,
                      stopcondition = u->(u[2,1] > 1e-1),
                      d = 4,
                      T = Float64,
                      args...)

    const x1 = one(T)
    const dx = (x1-x0)/npts
    r = exp(linspace(x0,x1,npts))
    # prepare initial data
    y0 = initial(uinit,r)
    u0 = reshape(y0,npts,2)
    ur0 = zero(u0)
    for i = 1:2
        ur0[:,i]=dur(r,u0[:,i])
    end

    # prepare the output arrays
      tout    = Array(T,1)
      rout    = Array(Array{T,1},1)
      uout    = Array(Array{T,2},1)
     urout    = Array(Array{T,2},1)
      tout[1] = zero(T)
      rout[1] = r
      uout[1] = u0
     urout[1] = ur0

    dy0 = zero(y0)

    while(true)

        if x0 < x0min
            break
        end

        F  = newF(x0,x1,d)

        for (t,y,dy) in dasslIterator(F, y0, tout[end]; dy0 = dy0, args...)
            print("$t\r")
            u  = reshape(y,npts,2)
            ur = zero(u)
            for j=1:2
                ur[:,j]=dur(r,u[:,j])
            end

            push!(tout,   t  )
            push!(rout,   r  )
            push!(uout,   u  )
            push!(urout,  ur )

            if stopcondition(u)
                break
            end
        end

        x0 = x0-dx
        npts = npts+1
        r  = [exp(x0), rout[end]]
        # r = exp(linspace(x0,x1,npts))
        u  = uout[end]
        du = uout[end]
        # interpolate
        u  = [[interpolate( u[:,1],dx), u[:,1]] [interpolate( u[:,2],dx), u[:,2]]]
        du = [[interpolate(du[:,1],dx),du[:,1]] [interpolate(du[:,2],dx),du[:,2]]]
        y0  = reshape( u,2*npts)
        dy0 = reshape(du,2*npts)
        println("npts=$npts, t=$(tout[end]), log(r[1])=$(log(r[1]))")
    end

    return tout, rout, uout, urout

end


function interpolate(u,dx)
    # second order interpolation of the newly created point
    return exp(-3*dx)*(exp(dx)*(1 + exp(dx))*u[1] - u[2])
end
