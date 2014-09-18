include("plot.jl")

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
