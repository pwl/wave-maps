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

function plotmode(r,u,ur,s,s0;d=3,args...)
    i0   = indmin(abs(s.-s0))
    y    = r[i0].*exp(s[i0])
    fsol = 2*atan(y/sqrt(d-2))
    Du   = u[i0][:,1]-fsol
    Dur0 = ur[i0][1,1]*exp(-s[i0])-2/sqrt(d-2)
    println(Dur0)
    mode = getmode(d,1)
    plot(y,Du,"x";title="First mode profile at s=$s0",ylabel="v_1(y)",xlabel="y",args...)
    oplot(mode[:,1],mode[:,2]*Dur0,"r:",xrange=[0:1])
end

function plotu(r,u,ur,s,s0;d=3,args...)
    i0   = indmin(abs(s.-s0))
    y    = r[i0].*exp(s[i0])
    fsol = 2*atan(y/sqrt(d-2))
    Du   = u[i0][:,1]
    plot(y[2:end],Du[2:end];title="Self-similar profile at s=$s0",xlabel="y",ylabel="f(y,s)",args...)
    oplot(y[2:end],fsol[2:end],"r:")
end

function getmode(d,n)
    f=open("modes/mode_d$(d)_n$(n).dat","r")
    mode = readdlm(f)
    close(f)
    return mode
end

function summarize(tau,t,r,u,ur,fn)
    T=t[end]+2/ur[end][1,1]       # crude estimate on blow-up time
    s=-log(T-t)

    # plot the mesh
    # p=plotmesh(tau,r,xrange=[1e-10,pi],yrange=[0,10])
    p1=plotmesh(s,r,xrange=[1e-8,pi],yrange=[0,20])

    # plot the eigenmode
    p2=plotmode(r,u,ur,s,s0,d=d) # @todo: fixed dimension

    # plot the self-similar profile
    p3=plotu(r,u,ur,s,s0,d=d,xlog=true,xrange=[1e-2,100],yrange=[0,pi]) # @todo: fixed dimension

    # plot the blow-up
    p4=plot(s,abs(map(z->z[1,1],ur).*exp(-s).-sqrt(d-2)),title="f'(0,s)-sqrt(d-2)",ylog=true,xlabel="-log(T-t)")

    tab=Table(2,2)
    tab[1,1]=p1
    tab[1,2]=p2
    tab[2,1]=p3
    tab[2,2]=p4
    return tab
end
