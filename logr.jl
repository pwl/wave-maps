const d=4
const (x0,x1)=(-10.,1.)

function initial(npts,u0::Function)
    x=linspace(x0,x1,npts)
    v1=map(u0,exp(x))
    v2=zero(x)
    y0=reshape([v1 v2], 2*npts)
    return y0
end

function rhs(t,y,dy)
    res = zero(y)
    npts=int(length(y)/2)
    dv=reshape(dy,npts,2)
     v=reshape( y,npts,2)
    x=linspace(x0,x1,npts)
    dx=x[2]-x[1]
    vxi=vx(v[:,1])
    vxxi=vxx(v[:,1])
    resv = reshape(res,npts,2)
    resv[:,1]=v[:,2]-dv[:,1]
    resv[:,2]=vxxi+(d-2)*vxi-(d-1)/2*sin(2*v[:,1])-dv[:,2].*exp(2*x)
    resv[1,2]=2*(d-1)*v[1,1]^3-3*(2+d)*exp(-dx)*(coth(dx)-1)*(exp(dx)*v[1,1]-v[2,1])-3*exp(2*x0)*dv[1,2]
    return res
end

function vx(v)
    vxi = zero(v)
    npts = length(v)
    dx = (x1-x0)/(npts-1)
    for i = 1:npts
        if i == 1
            vxi[i]=0
        elseif i == npts
            vxi[i]= (3*v[i]-4*v[i-1]+v[i-2])/(2dx)
        else
            vxi[i]=(v[i+1]-v[i-1])/(2dx)
        end
    end
    return vxi
end

function vxx(v)
    vxxi = zero(v)
    npts = length(v)
    dx = ((x1-x0)/(npts-1))
    vxxi[2:end-1] = (v[1:end-2]-2*v[2:end-1]+v[3:end])/dx^2
    vxxi[1]=0
    vxxi[end]=vxxi[end-1]
    return vxxi
end

function extract(t,y,dy)
    npts = int(length(y[1])/2)
     x = linspace(x0,x1,npts)
     r = exp(x)
     v = map(yi->reshape(yi,npts,2)[:,1],y)
    du = map(vi->exp(-x).*vx(vi),v)
    return r,v,du             # v=u
end
