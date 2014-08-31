using Winston

function plotmesh(s,r;args...)
    maxpts = 100
    npts   = length(r[1])
    di     = max(1,int(npts/maxpts))

    p = plot(xlog=true,xrange=[r[end][2],r[end][end]];title="Mesh points",args...)
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
    plot(y,Du;title="First mode profile",args...)
    oplot(mode[:,1],mode[:,2]*Dur0,xrange=[0:1])
end

function plotu(r,u,ur,s,s0;d=3,args...)
    i0   = indmin(abs(s.-s0))
    y    = r[i0].*exp(s[i0])
    fsol = 2*atan(y/sqrt(d-2))
    Du   = u[i0][:,1]
    plot(y[2:end],Du[2:end];title="Self-similar profile",args...)
    oplot(y[2:end],fsol[2:end],"r:")
end

function getmode(d,n)
    f=open("modes/mode_d$(d)_n$(n).dat","r")
    mode = readdlm(f)
    close(f)
    return mode
end
