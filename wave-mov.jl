using Winston
using DASSL
using ArrayViews

npde = 2

type Rhs
    fun  :: Function
    npde :: Int
end

function fun(r,u)
    dudt = zero(u)
    d = 3
    Lu = L(d,r,view(u,:,1))
    for i=1:size(u,1)
        dudt[i,1] = u[i,2]
        dudt[i,2] = Lu[i]-(d-1)/2*sin(2*u[i,1])/r[i]^2
    end

    return dudt
end

const myRhs=Rhs(fun,2)

function monitorarclength(r,u,ur)
    M=sqrt(1.+ur.^2)
end

function gnull(u,ur)
    1/sqrt(ur[1]^2+1)
end

function newF(rhs::Rhs;epsilon=1e-5,monitor=monitorarclength,g=gnull,args...)

    function F(tau,y,dy)


        npde = rhs.npde

        ny   = length(y)
        npts = int((ny-1)/(npde+1))
        dxi  = 1/(npts-1)

        res  = zero(y)

        # warning, t and dt are _copies_ of y[1] and dy[1]
         t  =   y[1]
        dt  =  dy[1]

         r  = view(  y,2:npts+1)
        dr  = view( dy,2:npts+1)
        rr  = view(res,2:npts+1)

         u  = reshape_view( view(  y,npts+2:ny), (npts, npde))
        du  = reshape_view( view( dy,npts+2:ny), (npts, npde))
        ru  = reshape_view( view(res,npts+2:ny), (npts, npde))

        dudr = zero(u)
        dudt = zero(u)
        for j = 1:npde
            dudr[:,j] = dur(r,view(u,:,j))
            dudt[:,j] = view(du,:,j)-dr.*dudr[:,j]
        end
        rhsval = rhs.fun(r,u)

        M    = monitor(r,view(u,:,1),view(dudr,:,1))
        gval = g(view(u,:,1),view(dudr,:,1))

        # sundman transform equations
        res[1] = dt-gval

        for i = 2:npts-1
            # moving mesh equations
            rr[i]   = 1/dxi^2*(dr[i-1]+dr[i+1]-2*dr[i])-gval / epsilon / dxi^2*( (M[i+1]+M[i])*(r[i+1]-r[i])-(M[i]+M[i-1])*(r[i]-r[i-1]) )
            # physical equations
            for j = 1:npde
                ru[i,j] = -dudt[i,j] + gval*rhsval[i,j]
            end
        end

        # mesh boundary conditions
        rr[    1] = dr[    1]
        rr[  end] = dr[  end]

        # function boundary conditions
        for j = 1:npde
            ru[1,  j] = du[1,  j]
            ru[end,j] = du[end,j]
        end

        return res
    end

    return F

end

function meshinit(r0,
                  uinit;
                  monitor = monitorarclength,
                  epsilon = 1e-5,
                  args...)

    info("Mesh initialization: Start")
    T = eltype(r0)
    npts = length(r0)
    dxi = T(1/(npts-1))

    function F(t,r,dr)
        u    = uinit(r)
        M    = monitor(r,u,dur(r,u))
        res  = copy(dr)
        for i = 2:npts-1
            res[i]=-epsilon*(dr[i-1]+dr[i+1]-2*dr[i])-( (M[i+1]+M[i])*(r[i+1]-r[i])-(M[i]+M[i-1])*(r[i]-r[i-1]) )
        end
        return res
    end

    err = Inf
    # for (t,r,dr) in dasslIterator(F, r0, zero(epsilon); reltol = dxi*epsilon, abstol = 1e-2*dxi*epsilon, args...)
    for (t,r,dr) in dasslIterator(F, r0, zero(T); reltol = dxi, abstol = 1e-2*dxi)
        err = norm(dr)*dxi
        if err < eps(T)^(2/3)
            info("Mesh initialization: Done! Total error = $err")
            return r
        end
    end

    warn("Mesh initialization: Failed! Total error = $err")
    return r0
end

function wavesolve(rhs   :: Rhs,
                   uinit :: Function;
                   npts  = 50,
                   rspan = [0,pi],
                   T     = Float64,
                   args...)
    npde = rhs.npde

    # TODO remove the generation of initial data to another function

    # generate initial data
    t0 = zero(T)
    dr = T((rspan[2]-rspan[1])/(npts-1))
    r0 = [rspan[1]:dr:rspan[2]]
    r0 = meshinit(r0, uinit; args...)
    u0 = zeros(npts,npde)
    u0[:,1]=uinit(r0)
    ur0=zero(u0)
    for j=1:npde
        ur0[:,j]=dur(r0,view(u0,:,j))
    end

    # generate y0
    ny=npts*(npde+1)+1
    y0=zeros(T,ny)
    y0[1]=t0
    y0[2:npts+1]=r0
    u0view=reshape_view(view(y0,npts+2:ny),(npts,npde))
    for j=1:npde, i=1:npts
        u0view[i,j]=u0[i,j]
    end

    F=newF(rhs; args...)

    # prepare the output arrays
    tauout = Array(T,1)
      tout = Array(T,1)
      rout = Array(Array{T,1},1)
      uout = Array(Array{T,2},1)
     urout = Array(Array{T,2},1)
    tauout[1] = zero(T)
      tout[1] = t0
      rout[1] = r0
      uout[1] = u0
     urout[1] = ur0

    for (tau,y,dy) in dasslIterator(F, y0, tauout[1]; args...)
        print("$tau\r")
        t  = y[1]
        r  = y[2:npts+1]
        u  = reshape(y[npts+2:end],npts,npde)
        ur = zero(u)
        for j=1:npde
            ur[:,j]=dur(r,view(u,:,j))
        end

        push!(tauout, tau)
        push!(tout,   t  )
        push!(rout,   r  )
        push!(uout,   u  )
        push!(urout,  ur )

        if u[2,1] > 1/2
            break
        end
    end

    return tauout, tout, rout, uout, urout

end


# dur and L assume:
# - r[1]=0
# - the function u is symmetric around r=0, i.e. around r[1]

function dur(r,u)
    ur=zero(u)
    npts=length(r)
    for i = 1:npts
        if i == 1
            ur[i]=(r[2+i]^2*(u[i]-u[1+i])+2r[i]*(r[2+i]*(-u[i]+u[1+i])+r[1+i]*(u[i]-u[2+i]))+r[1+i]^2*(-u[i]+u[2+i])+r[i]^2*(-u[1+i]+u[2+i]))/((r[i]-r[1+i])*(r[i]-r[2+i])*(r[1+i]-r[2+i]))
        elseif i == npts
            ur[i]=(r[i]^2*(-u[-2+i]+u[-1+i])+2r[-1+i]*r[i]*(u[-2+i]-u[i])+r[-2+i]^2*(u[-1+i]-u[i])+r[-1+i]^2*(-u[-2+i]+u[i])+2r[-2+i]*r[i]*(-u[-1+i]+u[i]))/((r[-2+i]-r[-1+i])*(r[-2+i]-r[i])*(r[-1+i]-r[i]))
        else
            ur[i]=(r[1+i]^2*(u[-1+i]-u[i])-2r[i]*(r[1+i]*(u[-1+i]-u[i])+r[-1+i]*(u[i]-u[1+i]))+r[i]^2*(u[-1+i]-u[1+i])+r[-1+i]^2*(u[i]-u[1+i]))/((r[-1+i]-r[i])*(r[-1+i]-r[1+i])*(r[i]-r[1+i]))
        end
    end

    return ur
end

function L(d,r,u;order=2)
    npts = length(r)
    Lu=zero(u)
    if order == 1
        for i = 1:npts

        end
    elseif order == 2
        for i = 1:npts
            if i == 1
                Lu[i]=0
            elseif i == 2
                Lu[i]=-(r[i]^(-2-d)*(r[i]^d*(2r[i]-r[1+i])*r[1+i]*(-2r[i]*r[1+i]*u[i]+r[1+i]^2*u[i]+r[i]^2*u[1+i])+(r[i]^4*r[1+i]^d*(r[2+i]^2*(u[i]-u[1+i])-2r[1+i]*(r[2+i]*(u[i]-u[1+i])+r[i]*(u[1+i]-u[2+i]))+r[1+i]^2*(u[i]-u[2+i])+r[i]^2*(u[1+i]-u[2+i])))/((r[i]-r[2+i])*(r[1+i]-r[2+i]))))/((r[i]-r[1+i])^2*r[1+i]^2)
            elseif i == npts-1
                Lu[i]=-((r[i]^(1-d)*((r[i-1]^(-1+d)*(r[i-1]-r[1+i])*(-r[i]+r[1+i])^3*(r[i]^2*(u[i-2]-u[i-1])-2r[i-1]*(r[i]*(u[i-2]-u[i-1])+r[i-2]*(u[i-1]-u[i]))+r[i-1]^2*(u[i-2]-u[i])+r[i-2]^2*(u[i-1]-u[i])))/((r[i-2]-r[i-1])*(r[i-2]-r[i])*(r[i-1]-r[i]))+r[i]^(-1+d)*(-1+(r[i]-r[1+i])/(r[i-1]-r[i]))*(r[i-1]-r[1+i])*(r[1+i]^2*(u[i-1]-u[i])-2r[i]*(r[1+i]*(u[i-1]-u[i])+r[i-1]*(u[i]-u[1+i]))+r[i]^2*(u[i-1]-u[1+i])+r[i-1]^2*(u[i]-u[1+i]))+(r[i-1]-r[i])*r[1+i]^(-1+d)*(r[1+i]^2*(-u[i-1]+u[i])+2r[i]*r[1+i]*(u[i-1]-u[1+i])+r[i-1]^2*(u[i]-u[1+i])+r[i]^2*(-u[i-1]+u[1+i])+2r[i-1]*r[1+i]*(-u[i]+u[1+i]))))/((r[i-1]-r[i])*(r[i-1]-r[1+i])^2*(r[i]-r[1+i])^2))
            elseif i == npts
                Lu[i]=(r[i]^(1-d)*(-((r[i-2]^(d-1)*(r[i-1]-r[i])^3*(r[i-1]^2*(u[i-3]-u[i-2])-2r[i-2]*(r[i-1]*(u[i-3]-u[i-2])+r[i-3]*(u[i-2]-u[i-1]))+r[i-2]^2*(u[i-3]-u[i-1])+r[i-3]^2*(u[i-2]-u[i-1])))/((r[i-3]-r[i-2])*(r[i-3]-r[i-1])))+r[i-1]^(-1+d)*(r[i-2]-r[i])*(r[i]^2*(u[i-2]-u[i-1])-2r[i-1]*(r[i]*(u[i-2]-u[i-1])+r[i-2]*(u[i-1]-u[i]))+r[i-1]^2*(u[i-2]-u[i])+r[i-2]^2*(u[i-1]-u[i]))-1/(r[i-2]-r[i])*(r[i-2]-r[i-1])*(r[i-2]+r[i-1]-2r[i])*r[i]^(d-1)*(r[i]^2*(-u[i-2]+u[i-1])+2r[i-1]*r[i]*(u[i-2]-u[i])+r[i-2]^2*(u[i-1]-u[i])+r[i-1]^2*(-u[i-2]+u[i])+2r[i-2]r[i]*(-u[i-1]+u[i]))))/((r[i-2]-r[i-1])^2*(r[i-2]-r[i])*(r[i-1]-r[i])^2)
            else
                Lu[i]=-(1/((r[i]-r[1+i])^2))*r[i]^(1-d)*(-((r[i-1]^(-1+d)*(r[i]-r[1+i])^3*(r[i]^2*(u[i-2]-u[i-1])-2r[i-1]*(r[i]*(u[i-2]-u[i-1])+r[i-2]*(u[i-1]-u[i]))+r[i-1]^2*(u[i-2]-u[i])+r[i-2]^2*(u[i-1]-u[i])))/((r[i-2]-r[i-1])*(r[i-2]-r[i])*(r[i-1]-r[i])^2*(r[i-1]-r[1+i])))-(r[i]^(-1+d)*(r[i-1]-2r[i]+r[1+i])*(r[1+i]^2*(u[i-1]-u[i])-2r[i]*(r[1+i]*(u[i-1]-u[i])+r[i-1]*(u[i]-u[1+i]))+r[i]^2*(u[i-1]-u[1+i])+r[i-1]^2*(u[i]-u[1+i])))/((r[i-1]-r[i])^2*(r[i-1]-r[1+i]))+((-r[i-1]+r[i])r[1+i]^(-1+d)*(r[2+i]^2*(u[i]-u[1+i])-2r[1+i]*(r[2+i]*(u[i]-u[1+i])+r[i]*(u[1+i]-u[2+i]))+r[1+i]^2*(u[i]-u[2+i])+r[i]^2*(u[1+i]-u[2+i])))/((-r[i-1]+r[1+i])*(r[i]-r[2+i])*(r[1+i]-r[2+i])))
            end
        end
    end
    return Lu
end
