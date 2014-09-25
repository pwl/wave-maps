type YMEquation1 <: YMEquation
    rhs  :: Function
    npde :: Int
    sundman :: Function
    monitor :: Function
    d :: Int

    profile :: Function
    a
    b
    phi0_yy

    function YMEquation1(d::Int)

        function s(r,u,ur)
            urr0=zero(r[1:5])
            dur!(urr0,r[1:5],ur[1:5])
            1/sqrt(abs(urr0[1])+1)
        end

        function m(r,u,ur)
            urr = zero(u); dur!(urr,r,ur)
            M = zero(ur)
            for i = 1:length(M)
                M[i]=sqrt(1+abs(ur[i])^2+abs(urr[i]))
            end
            return M
        end

        den = 3d-2sqrt(3d*(d-2))
        a =2(d+4-sqrt(3d*(d-2)))/den
        b =1/3*sqrt(3d*(d-2))*(8-d)/den

        profile(y)=a*y^2/(y^2+b)

        this = new()
        this.rhs = rhsYM(d)
        this.npde = 2
        this.d = d
        this.sundman = s
        this.monitor = m
        this.profile = profile
        this.a = a
        this.b = b
        this.phi0_yy = 2a/b
        return this
    end
end

function rhsYM(d)
    function rhs(r,u)
        npts = length(r)
        dudt = zero(u)
        u1 = view(u,:,1)
        Lu = L(d,r,u1,order=2)
        for i=2:npts-1
            dudt[i,1] = u[i,2]
            dudt[i,2] = Lu[i]-d*u1[i]*(1-u1[i])*(2-u1[i])/r[i]^2
        end
        return dudt
    end
    return rhs
end
