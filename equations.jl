const d = 4                          # dimension

function rhsWM(r,u)
    npts = length(r)
    dudt = zero(u)
    Lu = L(d,r,u[:,1],order=2)
    for i=2:npts-1
        dudt[i,1] = u[i,2]
        dudt[i,2] = Lu[i]-(d-1)/2*sin(2*u[i,1])/r[i]^2
        # dudt[i,2] = Lu[i]-(d-1)*u[i,1]/r[i]^2
    end
    return dudt
end

function rhsWMv2(r,u)
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

function rhsWMv3(r,u)
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

const WM=Equation(rhsWMv3,2)

# default monitor function
function WMmonitor(r,u,ur)
    M=sqrt(1.+abs(ur).^2)
end

function WMsundman(u,ur)
    1/sqrt((ur[1]*sqrt(d-2)/2)^2+1)
end
