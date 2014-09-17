using Winston

function plotmesh(res::Results,eqn::Equation;args...)
    maxpts = 100
    di     = max(1,int(res.npts/maxpts))

    p = plot(xlog=true,
             xrange=[res.r[end][2],res.r[end][end]];
             title="Mesh points",
             xlabel="r",
             ylabel="-log(T-t)",
             args...)

    for i = 2:di:res.npts
        oplot(map(q->q[i],res.r),res.s)
    end

    return p
end

# plot the solution profile at the given s0
function plotu(res::Results,eqn::Equation,s0::Real;args...)
    i0   = indmin(abs(res.s.-s0))
    y    = res.r[i0][2:end].*exp(res.s[i0])
    Du   = res.u[i0][2:end,1]

    p = plot(y,Du,".";title="Self-similar profile at s=$s0",xlabel="y",ylabel="u(y,s)",args...)

    # add the asymptotical profile, if available
    if isdefined(eqn,:profile)
        fsol = map(eqn.profile,y)
        oplot(y,fsol,"r")
    end

    return p
end

function plotspan(res::Results,eqn::Equation,plot::Function,s::Real;args...)
    plot(res,eqn,s;args...)
end

function plotspan(res::Results,eqn::Equation,plot::Function,s::Range;args...)
    rows = int(ceil(sqrt(length(s))))
    cols = int(ceil(length(s)/rows))
    t = Table(rows,cols)
    for i = 1:length(s)
        col = mod(i-1,cols)+1
        row = int((i-col)/cols)+1
        t[row,col] = plotspan(res,eqn,plot,s[i];title="s=$(s[i])",args...)
    end
    return t
end

function plotspan(res::Results,eqn::Equation,plot::Function;ns=9,args...)
    s = res.s
    rng = s[1]:((s[end-1]-s[1])/(ns-1)):s[end-1]
    return plotspan(res,eqn,plot,rng;args...)
end
