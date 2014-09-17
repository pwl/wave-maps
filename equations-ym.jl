type YMEquation <: Equation
    rhs  :: Function
    npde :: Int
    sundman :: Function
    monitor :: Function
    d :: Int

    function YMEquation(d::Int)

        function s(r,u,ur)
            urr=dur(r,ur)
            1/sqrt(abs(urr[1])+1)
        end

        function m(r,u,ur)
            urr=dur(r,ur)
            M=sqrt(1.+abs(ur).^2+abs(urr))
        end

        this = new()
        this.rhs=rhsYM(d)
        this.npde=2
        this.d=d
        this.sundman = s
        this.monitor = m
        return this
    end
end

function rhsYM(d)
    function rhs(r,u)
        npts = length(r)
        dudt = zero(u)
        Lu = L(d,r,u[:,1],order=2)
        f(u)=d*u*(1-u)*(2-u)
        for i=2:npts-1
            dudt[i,1] = u[i,2]
            dudt[i,2] = Lu[i]-f(u[i,1])/r[i]^2
        end
        return dudt
    end
    return rhs
end
