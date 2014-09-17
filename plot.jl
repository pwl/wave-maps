using Winston

function plotmesh(s,r;args...)
    maxpts = 100
    npts   = length(r[1])
    di     = max(1,int(npts/maxpts))

    p = plot(xlog=true,xrange=[r[end][2],r[end][end]];title="Mesh points",xlabel="r",ylabel="-log(T-t)",args...)
    for i = 2:di:npts
        oplot(map(q->q[i],r),s)
    end
    return p
end

function plotmeshu(tau,t,r,u,ur,d;args...)
    T,s  = Ts(tau,t,r,u,ur,d)
    maxpts = 100
    npts   = length(r[1])
    di     = max(1,int(npts/maxpts))

    p = plot(xrange=[0,5];title="Values of u at mesh points",xlabel="r",ylabel="-log(T-t)",args...)
    for i = 2:di:npts
        oplot(map(q->q[i,1],u),s)
    end
    return p
end

function plotxiu(tau,t,r,u,ur,s0,d;args...)
    T,s  = Ts(tau,t,r,u,ur,d)
    i0   = indmin(abs(s.-s0))
    xi   = linspace(0,1,length(r[1]))

    plot(xi,u[i0][:,1],".")
end

function plotmode(r,u,ur,s,s0,d;xrange=[0:2],args...)
    i0   = indmin(abs(s.-s0))
    y    = r[i0].*exp(s[i0])
    fsol = 2*atan(y/sqrt(d-2))
    Du   = u[i0][:,1]-fsol
    Dur0 = ur[i0][1,1]*exp(-s[i0])-2/sqrt(d-2)
    mode = getmode(d,1)
    plot(y,Du,".";title="First mode profile at s=$s0",ylabel="u(y(T-t),-log(T-t))-ϕ(y)",xlabel="y",xrange=xrange,args...)
    oplot(mode[:,1],mode[:,2]*Dur0,"r")
end

function plotu(tau,t,r,u,ur,s0,d;args...)
    T,s  = Ts(tau,t,r,u,ur,d)
    i0   = indmin(abs(s.-s0))
    y    = r[i0].*exp(s[i0])
    fsol = 2*atan(y/sqrt(d-2))
    Du   = u[i0][:,1]
    plot(y[2:end],Du[2:end],".";title="Self-similar profile at s=$s0",xlabel="y",ylabel="u(y(T-t),-log(T-t))",args...)
    oplot(y[2:end],fsol[2:end],"r")
end

function plotconvergencerate(tau,t,r,u,ur,d;fit=false,T=0,args...)
    if T == 0
        T,s  = Ts(tau,t,r,u,ur,d)
    else
        s = -log(T-t)
    end
    p = plot(title="(T-t)u_r(0,t)-ϕ'(0)",ylog=true,xlabel="-log(T-t)";args...)
    num = Curve(s,abs(map(z->z[1,1],ur).*exp(-s).-2/sqrt(d-2)))
    add(p,num)
    if fit
        (c,lambda)=fitlambda(tau,t,r,u,ur,d,T;args...)
        fit = Curve(s,exp(c+s*lambda),color="red")
        setattr(num,label="Numerical data")
        setattr(fit,label="Fit: Ce^λ^s, C=$(round(exp(c),4)), λ=$(round(lambda,4))")
        l = Legend(.1,.2,{num,fit})
        add(p,fit)
        add(p,l)
    end
    return p
end

function getmode(d,n)
    f=open("modes/mode_d$(d)_n$(n).dat","r")
    mode = readdlm(f)
    close(f)
    return mode
end

function summarize(tau,t,r,u,ur,s0,d)
    tab=Table(2,2)

    T,s=Ts(tau,t,r,u,ur,d)

    # plot the mesh
    # p=plotmesh(tau,r,xrange=[1e-10,pi],yrange=[0,10])
    tab[1,1]=plotmesh(s,r,xrange=[1e-8,pi],yrange=[0,20])

    # plot the eigenmode
    tab[1,2]=plotmode(r,u,ur,s,s0,d,xrange=[0,2])

    # plot the convergence rate to the self-similar solution
    tab[2,1]=plotconvergencerate(tau,t,r,u,ur,d)

    # plot the self-similar profile
    tab[2,2]=plotu(tau,t,r,u,ur,s0,d,xlog=true,xrange=[1e-2,100],yrange=[-pi,1.3pi],xrange=[1e-2,1e4])

    return tab
end

function Ts(tau,t,r,u,ur,d)
    T=t[end]+2/sqrt(d-2)/ur[end][1,1]       # crude estimate on blow-up time
    s=-log(abs(T-t))
    return (T,s)
end

function fitlambda(tau,t,r,u,ur,d,T;filter=s->true,args...)
    ur0 = map(z->z[1,1],ur)
    x = -log(T-t)
    y = log(abs(ur0.*exp(-s).-2/sqrt(d-2)))
    xnew = x[map(filter,x)]
    ynew = y[map(filter,x)]
    (a,b)=linreg(xnew,ynew)
end
