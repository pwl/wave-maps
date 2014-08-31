const d = 4                          # dimension

function rhsWM(r,u)
    npts = length(r)
    dudt = zero(u)
    Lu = L(d,r,u[:,1],order=1)
    dudt[1,:]=0
    for i=2:npts-1
        dudt[i,1] = u[i,2]
        dudt[i,2] = Lu[i]-(d-1)/2*sin(2*u[i,1])/r[i]^2
        # dudt[i,2] = Lu[i]-(d-1)*u[i,1]/r[i]^2
    end
    return dudt
end

const WM=Equation(rhsWM,2)

# default monitor function
function WMmonitor(r,u,ur)
    M=sqrt(1.+abs(ur).^2*(d-2))
end

function WMsundman(u,ur)
    1/sqrt(ur[1]^2+1)
end
