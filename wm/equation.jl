type WMEquation <: Equation
    rhs  :: Function
    npde :: Int
    sundman :: Function
    monitor :: Function
    d :: Int
    profile :: Function

    # derivate of the profile at y=0
    dphi0 :: Float64

    function WMEquation(d::Int,version::Int)
        s(r,u,ur)=1/sqrt(ur[1]^2+1)
        m(r,u,ur)=sqrt(1.+abs(ur).^2)
        phi0(y)=2*atan(y/sqrt(d-2))

        this = new()
        this.rhs=rhsWM(d,version)
        this.npde=2
        this.d=d
        this.sundman = s
        this.monitor = m
        this.profile = phi0
        this.dphi0 = 2/sqrt(d-2)
        return this
    end
end

# compute the blow-up time for this particular equation
function blowuptime(res::Results,eqn::WMEquation)
    d  = eqn.d
    t  = res.t[end]
    ur = res.ur[end][1,1]
    T=t+2/sqrt(d-2)/ur # crude estimate on blow-up time
end

function rhsWM(d::Int,version::Int)
    function rhsv1(r,u)
        npts = length(r)
        dudt = zero(u)
        Lu = L(d,r,view(u,:,1),order=2)
        for i=2:npts-1
            dudt[i,1] = u[i,2]
            dudt[i,2] = Lu[i]-(d-1)/2*sin(2*u[i,1])/r[i]^2
        end
        return dudt
    end

    function rhsv2(r,u)
        npts = length(r)
        dudt = zero(u)
        u_rr   = durr(r,u[:,1])
        ubyr   = u[:,1]./r
        ubyr[1]= dur(r[1:3],u[1:3,1])[1]
        ubyr_r = (d-1)*dur(r,ubyr)
        nonln  = (d-1)*(2u[:,1]-sin(2u[:,1]))./r.^2/2

        dudt[:,1] = u[:,2]
        dudt[:,2] = u_rr+ubyr_r+nonln

        return dudt
    end

    function rhsv3(r,u)
        npts = length(r)
        u1   = view(u,:,1)
        Lu   = Lu2(2,r,u1)
        ubyr   = u1./r
        ubyr[1]= dur(r[1:3],u[1:3,1])[1]
        ubyr_r = (d-2)*dur(r,ubyr)
        nonln  = (d-1)*(2u1-sin(2u1))./r.^2/2
        dudt = zero(u)
        dudt[:,1] = u[:,2]
        dudt[:,2] = Lu+ubyr_r+nonln
        dudt[1,:] = 0
        dudt[end,:] = 0
        return dudt
    end

    if     version == 1 return rhsv1
    elseif version == 2 return rhsv2
    elseif version == 3 return rhsv3
    else                return rhsv3
    end
end

function plotmode(res::Results,eqn::WMEquation,s0;xrange=[0:2],args...)
    i0   = indmin(abs(res.s.-s0))
    s    = res.s[i0]
    y    = res.y[i0]
    u    = res.u[i0][:,1]
    ur0  = res.ur[i0][1,1]
    fsol = map(eqn.profile,y)
    Du   = u-fsol
    Dur0 = ur0*exp(-s)-eqn.dphi0
    mode = getmode(eqn)
    p = plot(y,Du,".";title="First mode profile at s=$s0",xlabel="y",xrange=xrange,args...)
    oplot(mode[:,1],mode[:,2]*Dur0,"r")
    return p
end

function getmode(eqn::WMEquation)
    n = 1
    f = open("modes/mode_d$(eqn.d)_n$(n).dat","r")
    mode = readdlm(f)
    close(f)
    return mode
end

function plotconvergencerate(res::Results,eqn::WMEquation;fit=false,args...)

    s = res.s
    p = plot(title="(T-t)u_r(0,t)-ϕ'(0)",ylog=true,xlabel="-log(T-t)";args...)

    num = Curve(res.s,abs(map(z->z[1,1],res.ur).*exp(-s).-eqn.dphi0))

    add(p,num)
    if fit
        (c,lambda)=fitlambda(res,eqn;args...)
        fit = Curve(s,exp(c+s*lambda),color="red")
        setattr(num,label="Numerical data")
        setattr(fit,label="Fit: Ce^λ^s, C=$(round(exp(c),4)), λ=$(round(lambda,4))")
        l = Legend(.1,.2,{num,fit})
        add(p,fit)
        add(p,l)
    end
    return p
end

function fitlambda(res::Results,eqn::WMEquation;filter=s->true,args...)
    ur0 = map(z->z[1,1],res.ur)
    x = -log(res.T-res.t)
    y = log(abs(ur0.*exp(-res.s).-eqn.dphi0))
    xnew = x[map(filter,x)]
    ynew = y[map(filter,x)]
    (a,b)=linreg(xnew,ynew)
end

function summarize(res::Results,eqn::WMEquation,s0)
    tab=Table(2,2)

    # plot the mesh
    # p=plotmesh(tau,r,xrange=[1e-10,pi],yrange=[0,10])
    tab[1,1]=plotmesh(res,eqn,xrange=[1e-8,pi],yrange=[0,20])

    # plot the eigenmode
    tab[1,2]=plotmode(res,eqn,xrange=[0,2])

    # plot the convergence rate to the self-similar solution
    tab[2,1]=plotconvergencerate(res,eqn)

    # plot the self-similar profile
    tab[2,2]=plotu(res,eqn,xlog=true,s0,xrange=[1e-2,100],yrange=[-pi,1.3pi],xrange=[1e-2,1e4])

    return tab
end
