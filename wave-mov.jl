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
        dudt[i,2] = Lu[i]-(d-1)/2*sin(2*u[i,1])./r.^2
    end
    return dudt
end

const myRhs=Rhs(fun,2)

function monitorarclength(r,u,ur)
    M=sqrt(1.+ur.^2)
end

function gnull(u,ur)
    1
end

function newF(rhs::Rhs;eps=1e-5,monitor=monitorarclength,g=gnull)

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

        dudr = dur(r,u)
        dudt = du-dr.*dudr
        rhsval = rhs.fun(r,u)

        M    = monitor(r,view(u,:,1),view(dudr,:,1))
        gval = g(view(u,:,1),view(dudr,:,1))

        # sundman transform equations
        res[1] = dt-gval

        for i = 2:npts-1
            # moving mesh equations
            rr[i]   = 1/dxi^2*(dr[i-1]+dr[i+1]-2*dr[i])-gval / eps / dxi^2*( (M[i+1]+M[i])*(r[i+1]-r[i])-(M[i]+M[i-1])*(r[i]-r[i-1]) )
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

function meshinit(r0,u0;monitor=monitorarclength, eps=1e-5, args...)
    npts = length(r0)
    dxi = 1/(npts-1)

    function F(t,r,dr)
        u    = u0(r)
        dudr = dur(r,u)
        M    = monitor(r,u,dudr)
        res  = copy(dr)
        for i = 2:npts-1
            res[i]=dr[i]-1/ eps / dxi^2*( (M[i+1]+M[i])*(r[i+1]-r[i])-(M[i]+M[i-1])*(r[i]-r[i-1]) )
        end
        return res
    end

    (t,r)=dasslSolve(F,r0,[0.,10]; args...)

    return r[end]
end

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
