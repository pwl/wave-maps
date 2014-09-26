type YMEquation2 <: YMEquation
    rhs  :: Function
    npde :: Int
    sundman :: Function
    monitor :: Function
    d :: Int

    profile :: Function
    a
    b
    phi0_yy

    function YMEquation2(d::Int)

        function s(r,u,ur)
            amp = maximum(abs(u))
            # 1/sqrt(amp^2+1)
            1/sqrt(ur[1]^2+1)
        end

        function m(r,u,ur)
            M = zero(ur)
            urr = zero(ur); dur!(urr,r,ur)
            for i = 1:length(M)
                M[i]=sqrt(1+abs(u[i]+r[i]*ur[i])^2.1+abs(2ur[i]+r[i]*urr[i]))
                # M[i]=sqrt(1+abs(ur[i])+sqrt(abs(urr[i])))
            end
            return M
        end

        den = 3d-2sqrt(3d*(d-2))
        a =2(d+4-sqrt(3d*(d-2)))/den
        b =1/3*sqrt(3d*(d-2))*(8-d)/den

        profile(y)=a*y^2/(y^2+b)

        this = new()
        this.rhs = rhsYM2(d)
        this.npde = 2
        this.d = d
        this.sundman = s
        this.monitor = m
        this.profile = profile
        this.a = a
        this.b = b
        this.phi0_yy = 2*a/b
        return this
    end
end

function rhsYM2(d)
    function rhs(r,u)
        npts = length(r)
        dudt = zero(u)
        u1 = view(u,:,1)
        Lu = L(d+2,r,u1,order=2)
        for i=2:npts-1
            dudt[i,1] = u[i,2]
            dudt[i,2] = Lu[i]+u1[i]/r[i]^2*((d-1)-d*(1-r[i]*u1[i])*(2-r[i]*u1[i]))
        end
        return dudt
    end
end

function finalize!(res::Results,eqn::YMEquation2)
    res.npts = length(res.r[1])
    for n = 1:length(res.t)
        for j=1:eqn.npde
            for i=1:res.npts
                r   = res.r[n][i]
                v   = res.u[n][i,j]
                vr  = res.ur[n][i,j]
                vrr = res.urr[n][i,j]
                res.u[n][i,j]   = r*v
                res.ur[n][i,j]  = r*vr+v
                res.urr[n][i,j] = r*vrr+2*vr
            end
        end
    end

    res.T = blowuptime(res,eqn)
    res.s = -log(res.T.-res.t)
    for n = 1:length(res.t)
        push!(res.y,selfsimilary(res.r[n],res.s[n],eqn))
    end
end
