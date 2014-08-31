using DASSL
using ArrayViews

include("discretization.jl")

type Equation
    rhs  :: Function
    npde :: Int
end

function arclength(;power=2,A=1,B=1)
    function monitor(r,u,ur)
        M=sqrt(A.+B*abs(ur).^power)
    end
end

# default sundman transform function, suitable for self-similar blow-up in first derivative
function defaultsundman(u,ur)
    1/sqrt(1+ur[1]^2)
end

function extracty(y,dy,npts,npde)
    t   =  y[1]
    dt  = dy[1]
    r   =  y[2:npts+1]
    dr  = dy[2:npts+1]
    u   = reshape(  y[npts+2:end], npts, npde)
    du  = reshape( dy[npts+2:end], npts, npde)
    return t, dt, r, dr, u, du
end

function computederivatives(r,dr,u,du)
    dudr = zero(u)          # du/dr
    for j = 1:size(u,2)
        dudr[:,j] = dur(r,u[:,j])
    end
    dudt = du-dr.*dudr      # du/dt=u_t-dr/dt*u_r
    return dudr, dudt
end

# simple gaussian filter to smoothen the monitor function
function smoothen(M; passes = 1, args...)
    if passes == 0
        return M
    else
        Msmooth = zero(M)
        Msmooth[2:end-1]=(M[1:end-2]+2M[2:end-1]+M[3:end])/4
        Msmooth[1]=(M[1]+M[2])/2
        Msmooth[end]=(M[end-1]+M[end])/2
        return smoothen(Msmooth, passes=passes-1)
    end
end

function movingmeshres(r,dr,M,gamma,gval,epsilon)
    npts = length(dr)
    dxi  = 1/(npts-1)
    resr = zero(dr)          # residua of moving mesh
    for i = 2:npts-1
        resr[i]   = (dr[i]-gamma*(dr[i-1]+dr[i+1]-2*dr[i])/dxi^2)-gval / epsilon*( (M[i+1]+M[i])*(r[i+1]-r[i])/dxi-(M[i]+M[i-1])*(r[i]-r[i-1])/dxi )/dxi
    end
    resr[ 1 ] =  r[ 1 ]   # mesh boundary conditions
    resr[end] = dr[end]
    return resr
end

function physicalres(rhs,r,u,dudt,gval)
    rhsval = rhs(r,u)
    npts, npde = size(dudt)
    resu = -dudt+gval.*rhsval
    resu[1,:] = u[1,:]          # enforce boundary conditions u(0,t)=0
    return reshape(resu,npts*npde)
end

function newF(eqn::Equation;
              epsilon=1e-5,
              monitor=arclength(),
              g=defaultsundman,
              gamma=2,
              args...)

    function F(tau,y,dy)
        npde = eqn.npde

        npts = int((length(y)-1)/(npde+1))

        t, dt, r, dr, u, du = extracty(y,dy,npts,npde)
        dudr, dudt = computederivatives(r,dr,u,du)

        M    = smoothen(monitor(r,u[:,1],dudr[:,1]); args...)
        gval = g(u[:,1],dudr[:,1])

        res  = zero(y)
        res[1]          = dt-gval                                  # Sundman transform equations
        res[2:npts+1]   = movingmeshres(r,dr,M,gamma,gval,epsilon) # moving mesh equations
        res[npts+2:end] = physicalres(eqn.rhs,r,u,dudt,gval)       # physical equations

        return res
    end

    return F

end

function meshinit(r0,
                  uinit;
                  monitor = arclength(),
                  epsilon = 1e-5,
                  gamma = 2,
                  args...)

    info("Mesh initialization: Start")
    T    = eltype(r0)
    npts = length(r0)
    dxi  = convert(T,1/(npts-1))

    function F(t,r,dr)
        u    = uinit(r)
        M    = smoothen(monitor(r,u,dur(r,u)); args...)
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

function wavesolve(eqn   :: Equation,
                   uinit :: Function;
                   npts  = 100,
                   rspan = [0,pi],
                   T     = Float64,
                   taumax = 20,
                   stopcondition = u->(u[2,1] > 1/5),
                   args...)
    npde = eqn.npde

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

    F=newF(eqn; args...)

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

        if stopcondition(u) || tau > taumax
            break
        end
    end

    return tauout, tout, rout, uout, urout

end

include("equations.jl")
