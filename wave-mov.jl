using Winston
using DASSL
using ArrayViews

include("discretization.jl")

type Rhs
    fun  :: Function
    npde :: Int
end

function fun(r,u)
    npts = length(r)
    dudt = zero(u)
    d = 10
    Lu = L(d,r,u[:,1])
    dudt[1,:]=0
    for i=2:npts-1
        dudt[i,1] = u[i,2]
        # dudt[i,2] = Lu[i]-(d-1)/2*sin(2*u[i,1])/r[i]^2
        dudt[i,2] = Lu[i]-(d-1)*u[i,1]/r[i]^2
    end
    return dudt
end

const myRhs=Rhs(fun,2)

function monitorarclength(r,u,ur)
    # @todo add mesh smoothing or some normalization?
    M=sqrt(1.+abs(ur).^(2.1))
end

function gnull(u,ur)
    1/sqrt(ur[1]^2+1)
end

function newF(rhs::Rhs;
              epsilon=1e-5,
              monitor=monitorarclength,
              g=gnull,
              gamma=2,
              args...)

    function F(tau,y,dy)
        npde = rhs.npde

        ny   = length(y)
        npts = int((ny-1)/(npde+1))
        dxi  = 1/(npts-1)

         t  =   y[1]
        dt  =  dy[1]

         r  =   y[2:npts+1]
        dr  =  dy[2:npts+1]
        rr  = zero(r)

         u  = reshape(  y[npts+2:ny], npts, npde)
        du  = reshape( dy[npts+2:ny], npts, npde)
        ru  = zero(u)

        dudr = zero(u)          # du/dr
        for j = 1:npde
            dudr[:,j] = dur(r,u[:,j])
        end
        dudt = du-dr.*dudr      # du/dt=u_t-dr/dt*u_r

        rhsval = rhs.fun(r,u)

        M    = monitor(r,u[:,1],dudr[:,1])
        gval = g(u[:,1],dudr[:,1])

        res  = zero(y)

        # Sundman transform equations
        res[1] = dt-gval

        # moving mesh equation
        for i = 2:npts-1
            rr[i]   = (dr[i]*dxi^2-gamma*(dr[i-1]+dr[i+1]-2*dr[i]))-gval / epsilon*( (M[i+1]+M[i])*(r[i+1]-r[i])-(M[i]+M[i-1])*(r[i]-r[i-1]) )
        end
        rr[ 1 ] = dr[ 1 ]   # mesh boundary conditions
        rr[end] = dr[end]
        res[2:npts+1]=rr

        # physical equations
        ru = -dudt+gval.*rhsval
        res[npts+2:end]=reshape(ru,npts*npde)

        return res
    end

    return F

end

function meshinit(r0,
                  uinit;
                  monitor = monitorarclength,
                  epsilon = 1e-5,
                  gamma = 2,
                  args...)

    info("Mesh initialization: Start")
    T    = eltype(r0)
    npts = length(r0)
    dxi  = convert(T,1/(npts-1))

    function F(t,r,dr)
        u    = uinit(r)
        M    = monitor(r,u,dur(r,u))
        res  = copy(dr)
        for i = 2:npts-1
            res[i]=epsilon*(dr[i]*dxi^2-gamma*(dr[i-1]+dr[i+1]-2*dr[i]))-( (M[i+1]+M[i])*(r[i+1]-r[i])-(M[i]+M[i-1])*(r[i]-r[i-1]) )
        end
        return res
    end

    tol = min(1,epsilon*1e7)
    err = Inf
    for (t,r,dr) in dasslIterator(F, r0, zero(T); reltol = dxi*tol, abstol = 1e-2*dxi*tol)
        err = norm(dr)*dxi
        if err < eps(T)^(4/5)
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
                   taumax = 30,
                   args...)
    npde = rhs.npde

    # TODO remove the generation of initial data to another function

    # generate initial data
    t0 = zero(T)
    dr = convert(T,(rspan[2]-rspan[1])/(npts-1))
    r0 = [rspan[1]:dr:rspan[2]]
    r0 = meshinit(r0, uinit; args...)
    u0 = zeros(npts,npde)
    u0[:,1]=uinit(r0)
    ur0=zero(u0)
    for j=1:npde
        ur0[:,j]=dur(r0,u0[:,j])
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
    tauout    = Array(T,1)
      tout    = Array(T,1)
      rout    = Array(Array{T,1},1)
      uout    = Array(Array{T,2},1)
     urout    = Array(Array{T,2},1)
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
            ur[:,j]=dur(r,u[:,j])
        end

        push!(tauout, tau)
        push!(tout,   t  )
        push!(rout,   r  )
        push!(uout,   u  )
        push!(urout,  ur )

        if u[2,1] > 1/5 || tau > taumax
            break
        end
    end

    return tauout, tout, rout, uout, urout

end
