using ArrayViews

abstract YMEquation <: Equation

include("equation-1.jl")
include("equation-2.jl")

function blowuptime(res::Results,eqn::YMEquation)
    T=res.t[end]+sqrt(eqn.phi0_yy/res.urr[end][1,1])
end

function plotconvergencerate(res::Results,eqn::YMEquation;fit=false,args...)

    s = res.s
    p = plot(title="(T-t)u_r(0,t)-ϕ'(0)",ylog=true,xlabel="-log(T-t)";args...)

    num = Curve(res.s,abs(map(z->z[1,1],res.urr).*exp(-2s).-eqn.phi0_yy))

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

function fitlambda(res::Results,eqn::YMEquation;filter=s->true,args...)
    urr0 = map(z->z[1,1],res.urr)
    x = -log(res.T-res.t)
    y = log(abs(urr0.*exp(-2res.s).-eqn.phi0_yy))
    xnew = x[map(filter,x)]
    ynew = y[map(filter,x)]
    (a,b)=linreg(xnew,ynew)
end

function summarize(res::Results,eqn::YMEquation,s0)
    tab=Table(2,2)

    # plot the mesh
    # p=plotmesh(tau,r,xrange=[1e-10,pi],yrange=[0,10])
    tab[1,1]=plotmesh(res,eqn,xrange=[1e-8,pi],yrange=[0,20])

    # plot mode is not implemented yet
    # plot the eigenmode
    # tab[1,2]=plotmode(res,eqn,xrange=[0,2])

    # plot the convergence rate to the self-similar solution
    tab[2,1]=plotconvergencerate(res,eqn)

    # plot the self-similar profile
    tab[2,2]=plotu(res,eqn,xlog=true,s0,xrange=[1e-2,100],yrange=[0,1.1eqn.a],xrange=[1e-2,1e2])

    return tab
end
