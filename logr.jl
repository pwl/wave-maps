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

        vxi=vx(v[:,1],dx)
        vxxi=vxx(v[:,1],dx)

        resv = reshape(res,npts,2)
        resv[:,1]=v[:,2]-dv[:,1]
        resv[:,2]=vxxi+(d-2)*vxi-(d-1)/2*sin(2*v[:,1])-dv[:,2].*exp(2*x)
        # The tricky part: boundary condition coming from expansion around
        # r=0 and replacing the derivatives ur and urrr with their
        # discretizations on a mesh r=[0,exp(x0),exp(x0+dx)] with
        # u=[0,v[1],v[2]]
        resv[1,2]=2*(d-1)*v[1,1]^3-3*(2+d)*exp(-dx)*(coth(dx)-1)*(exp(dx)*v[1,1]-v[2,1])-3*exp(2*x0)*dv[1,2]
        return res
    end
    return rhs
end

function vx(v,dx)
    vxi = zero(v)
    npts = length(v)
    for i = 1:npts
        if i == 1
            vxi[i]=0
        elseif i == npts
            vxi[i]= (3*v[i]-4*v[i-1]+v[i-2])/(2dx)
        else
            vxi[i]=(v[i+1]-v[i-1])/(2dx)
        end
    end
    return vxi
end

function vxx(v,dx)
    vxxi = zero(v)
    npts = length(v)
    vxxi[2:end-1] = (v[1:end-2]-2*v[2:end-1]+v[3:end])/dx^2
    vxxi[1]=0
    vxxi[end]=vxxi[end-1]
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
                      args...)

    const x1 = 1
    const dx = (x1-x0)/npts
    T=Float64
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

    while(r[1] > 1e-10)

        if x0 < x0min
            break
        end

        x0 = x0-dx
        npts = npts+1
        r  = [exp(x0), rout[end]]
        # r = exp(linspace(x0,x1,npts))
        u  = uout[end]
        # interpolate
        u = [[exp(-dx)*u[1,1],u[:,1]] [exp(-dx)*u[1,2],u[:,2]]]
        y0 = reshape(u,2*npts)
        F  = newF(x0,x1,d)
        println("npts=$npts, t=$(tout[end]), r[1]=$(r[1])")

        for (t,y,dy) in dasslIterator(F, y0, tout[end]; args...)
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
    end

    return tout, rout, uout, urout


end
