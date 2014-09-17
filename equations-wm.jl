type WMEquation <: Equation
    rhs  :: Function
    npde :: Int
    sundman :: Function
    monitor :: Function
    d :: Int

    function WMEquation(d::Int)

        function s(r,u,ur)
            urr=dur(r,ur)
            1/sqrt(ur[1].^2+1)
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

function rhsWM(d::Int,version::Int)
    function rhsv1(r,u)
        npts = length(r)
        dudt = zero(u)
        Lu = L(d,r,u[:,1],order=2)
        for i=2:npts-1
            dudt[i,1] = u[i,2]
            dudt[i,2] = Lu[i]-(d-1)/2*sin(2*u[i,1])/r[i]^2
        end
        return dudt
    end

    function rhsv2(r,u)
        npts = length(r)
        dudt = zero(u)
        u_rr   = durr(r,u[:,1])
        ubyr   = u[:,1]./r
        ubyr[1]= dur(r[1:3],u[1:3,1])[1]
        ubyr_r = (d-1)*dur(r,ubyr)
        nonln  = (d-1)*(2u[:,1]-sin(2u[:,1]))./r.^2 / 2

        dudt[:,1] = u[:,2]
        dudt[:,2] = u_rr+ubyr_r+nonln

        return dudt
    end

    function rhsv3(r,u)
        npts = length(r)
        Lu   = Lu2(2,r,u[:,1])
        ubyr   = u[:,1]./r
        ubyr[1]= dur(r[1:3],u[1:3,1])[1]
        ubyr_r = (d-2)*dur(r,ubyr)
        nonln  = (d-1)*(2u[:,1]-sin(2u[:,1]))./r.^2 / 2
        dudt = zero(u)
        dudt[:,1] = u[:,2]
        dudt[:,2] = Lu+ubyr_r+nonln
        dudt[1,:] = 0
        dudt[end,:] = 0
        return dudt
    end

    if     version == 1 return rhsv1
    elseif version == 2 return rhsv2
    elseif version == 3 return rhsv3
    else                return rhsv3
    end
end
