using Winston
using DASSL

npde = 2

function newF(d,rOld,uOld)

    r    = new_mesh(rOld)
    npts = length(r)
    npde = size(uOld,2)
    u0   = zeros(eltype(uOld),npts,npde)
    for i = 1:npde
        u0[:,i] = new_u(rOld,uOld[:,i],r)
    end

    function rhs(t,u)
        du = zero(u)
        du[  :,1] = u[:,2]
        du[  :,2] = L(d,r,u[:,1])-(d-1)/2*sin(2*u[:,1])./r.^2
        return du
    end

    function Fexp(t,y)
        dy  = zeros(npts,npde)
        u   = reshape(y  ,npts,npde)
        du  = reshape(dy ,npts,npde)
        du[  :,1] = u[:,2]
        du[  :,2] = L(d,r,u[:,1])-(d-1)/2*sin(2*u[:,1])./r.^2
        du[  1,:] = 0
        du[end,:] = 0
        return dy
    end

    function Fimp(t,y,dy)
        res = zeros(npts,npde)
        u   = reshape(y  ,npts,npde)
        du  = reshape(dy ,npts,npde)
        res = reshape(res,npts,npde)
        res[  :,1] = -du[:,1] + u[:,2]
        res[  :,2] = -du[:,2] + L(d,r,u[:,1])-(d-1)/2*sin(2*u[:,1])./r.^2
        res[  1,:] = du[  1,:]
        res[end,:] = du[end,:]
        return reshape(res,length(res))
    end

    return (r,u0,Fimp,Fexp)

end

function new_mesh(rOld)
    return [0, rOld[2]/2, rOld[2:end]]
end

function new_u(yOld,uOld,y)
    # linear interpolation does not use y
    du0 = uOld[2]/yOld[2]
    return [0, du0*yOld[2]/2, uOld[2:end]]
end

function wavesolve(d,r0,u0; C = 0.01, rmin = 1e-5, args...)
    status = 0
    T = eltype(r0)

    (r,u,Fi)=newF(d,r0,u0)
    while u[2,1] > C
        (r,u,Fi)=newF(d,r,u)
    end

    tout  = Array(T,1)
    rout  = Array(Array{T,1},1)
    uout  = Array(Array{T,2},1)
    urout = Array(Array{T,2},1)

     tout[1] = 0
     rout[1] = r
     uout[1] = u
    urout[1] = dur(r,u)

    while tout[end] < 10

        r = rout[end]
        u = uout[end]

        if r[2] < rmin
            break
        end

        for (t,y,dy) in dasslIterator(Fi, reshape(u,length(u)), tout[end]; args...)
            print("$t, $(r[2])\r")
            u = reshape(y,length(r),2)
            push!(tout,  t)
            push!(rout,  r)
            push!(uout,  u)
            push!(urout, dur(r,u))

            if u[2,1] > C || tout[end] > 10
                break
            end
        end

        (rout[end],uout[end],Fi)=newF(d,rout[end],uout[end])
    end

    return tout, rout, uout, urout

end

function dur(r,u)
    ur=zero(u)
    npts=length(r)
    npde=size(u,2)
    for j = 1:npde
        uj=u[:,j]
        for i = 1:npts
            if i == 1
                ur[i,j]=(r[2+i]^2*(uj[i]-uj[1+i])+2r[i]*(r[2+i]*(-uj[i]+uj[1+i])+r[1+i]*(uj[i]-uj[2+i]))+r[1+i]^2*(-uj[i]+uj[2+i])+r[i]^2*(-uj[1+i]+uj[2+i]))/((r[i]-r[1+i])*(r[i]-r[2+i])*(r[1+i]-r[2+i]))
            elseif i == npts
                ur[i,j]=(r[i]^2*(-uj[-2+i]+uj[-1+i])+2r[-1+i]*r[i]*(uj[-2+i]-uj[i])+r[-2+i]^2*(uj[-1+i]-uj[i])+r[-1+i]^2*(-uj[-2+i]+uj[i])+2r[-2+i]*r[i]*(-uj[-1+i]+uj[i]))/((r[-2+i]-r[-1+i])*(r[-2+i]-r[i])*(r[-1+i]-r[i]))
            else
                ur[i,j]=(r[1+i]^2*(uj[-1+i]-uj[i])-2r[i]*(r[1+i]*(uj[-1+i]-uj[i])+r[-1+i]*(uj[i]-uj[1+i]))+r[i]^2*(uj[-1+i]-uj[1+i])+r[-1+i]^2*(uj[i]-uj[1+i]))/((r[-1+i]-r[i])*(r[-1+i]-r[1+i])*(r[i]-r[1+i]))
            end
        end
    end

    return ur
end

function L(d,r,u;order=2)
    npts = length(r)
    Lu=zero(u)
    if order == 1
        for i = 1:npts

        end
    elseif order == 2
        for i = 1:npts
            if i == 1
                Lu[i]=0
            elseif i == 2
                Lu[i]=-(r[i]^(-2-d)*(r[i]^d*(2r[i]-r[1+i])*r[1+i]*(-2r[i]*r[1+i]*u[i]+r[1+i]^2*u[i]+r[i]^2*u[1+i])+(r[i]^4*r[1+i]^d*(r[2+i]^2*(u[i]-u[1+i])-2r[1+i]*(r[2+i]*(u[i]-u[1+i])+r[i]*(u[1+i]-u[2+i]))+r[1+i]^2*(u[i]-u[2+i])+r[i]^2*(u[1+i]-u[2+i])))/((r[i]-r[2+i])*(r[1+i]-r[2+i]))))/((r[i]-r[1+i])^2*r[1+i]^2)
            elseif i == npts-1
                Lu[i]=-((r[i]^(1-d)*((r[i-1]^(-1+d)*(r[i-1]-r[1+i])*(-r[i]+r[1+i])^3*(r[i]^2*(u[i-2]-u[i-1])-2r[i-1]*(r[i]*(u[i-2]-u[i-1])+r[i-2]*(u[i-1]-u[i]))+r[i-1]^2*(u[i-2]-u[i])+r[i-2]^2*(u[i-1]-u[i])))/((r[i-2]-r[i-1])*(r[i-2]-r[i])*(r[i-1]-r[i]))+r[i]^(-1+d)*(-1+(r[i]-r[1+i])/(r[i-1]-r[i]))*(r[i-1]-r[1+i])*(r[1+i]^2*(u[i-1]-u[i])-2r[i]*(r[1+i]*(u[i-1]-u[i])+r[i-1]*(u[i]-u[1+i]))+r[i]^2*(u[i-1]-u[1+i])+r[i-1]^2*(u[i]-u[1+i]))+(r[i-1]-r[i])*r[1+i]^(-1+d)*(r[1+i]^2*(-u[i-1]+u[i])+2r[i]*r[1+i]*(u[i-1]-u[1+i])+r[i-1]^2*(u[i]-u[1+i])+r[i]^2*(-u[i-1]+u[1+i])+2r[i-1]*r[1+i]*(-u[i]+u[1+i]))))/((r[i-1]-r[i])*(r[i-1]-r[1+i])^2*(r[i]-r[1+i])^2))
            elseif i == npts
                Lu[i]=(r[i]^(1-d)*(-((r[i-2]^(d-1)*(r[i-1]-r[i])^3*(r[i-1]^2*(u[i-3]-u[i-2])-2r[i-2]*(r[i-1]*(u[i-3]-u[i-2])+r[i-3]*(u[i-2]-u[i-1]))+r[i-2]^2*(u[i-3]-u[i-1])+r[i-3]^2*(u[i-2]-u[i-1])))/((r[i-3]-r[i-2])*(r[i-3]-r[i-1])))+r[i-1]^(-1+d)*(r[i-2]-r[i])*(r[i]^2*(u[i-2]-u[i-1])-2r[i-1]*(r[i]*(u[i-2]-u[i-1])+r[i-2]*(u[i-1]-u[i]))+r[i-1]^2*(u[i-2]-u[i])+r[i-2]^2*(u[i-1]-u[i]))-1/(r[i-2]-r[i])*(r[i-2]-r[i-1])*(r[i-2]+r[i-1]-2r[i])*r[i]^(d-1)*(r[i]^2*(-u[i-2]+u[i-1])+2r[i-1]*r[i]*(u[i-2]-u[i])+r[i-2]^2*(u[i-1]-u[i])+r[i-1]^2*(-u[i-2]+u[i])+2r[i-2]r[i]*(-u[i-1]+u[i]))))/((r[i-2]-r[i-1])^2*(r[i-2]-r[i])*(r[i-1]-r[i])^2)
            else
                Lu[i]=-(1/((r[i]-r[1+i])^2))*r[i]^(1-d)*(-((r[i-1]^(-1+d)*(r[i]-r[1+i])^3*(r[i]^2*(u[i-2]-u[i-1])-2r[i-1]*(r[i]*(u[i-2]-u[i-1])+r[i-2]*(u[i-1]-u[i]))+r[i-1]^2*(u[i-2]-u[i])+r[i-2]^2*(u[i-1]-u[i])))/((r[i-2]-r[i-1])*(r[i-2]-r[i])*(r[i-1]-r[i])^2*(r[i-1]-r[1+i])))-(r[i]^(-1+d)*(r[i-1]-2r[i]+r[1+i])*(r[1+i]^2*(u[i-1]-u[i])-2r[i]*(r[1+i]*(u[i-1]-u[i])+r[i-1]*(u[i]-u[1+i]))+r[i]^2*(u[i-1]-u[1+i])+r[i-1]^2*(u[i]-u[1+i])))/((r[i-1]-r[i])^2*(r[i-1]-r[1+i]))+((-r[i-1]+r[i])r[1+i]^(-1+d)*(r[2+i]^2*(u[i]-u[1+i])-2r[1+i]*(r[2+i]*(u[i]-u[1+i])+r[i]*(u[1+i]-u[2+i]))+r[1+i]^2*(u[i]-u[2+i])+r[i]^2*(u[1+i]-u[2+i])))/((-r[i-1]+r[1+i])*(r[i]-r[2+i])*(r[1+i]-r[2+i])))
            end
        end
    end
    return Lu
end
