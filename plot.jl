function plotmesh(tau,r)
    maxpts = 100
    npts   = length(r[1])
    di     = max(1,int(npts/maxpts))

    p = plot(xlog=true,xrange=[r[end][2],r[end][end]])
    for i = 2:di:npts
        oplot(map(q->q[i],r),tau)
    end
    display(p)
    return(p)
end
