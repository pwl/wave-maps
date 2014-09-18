# dur and L assume:
# - r[1]=0
# - the function u is symmetric around r=0, i.e. around r[1]

function dur!(ur,r,u)
    # order = 2
    npts=length(r)
    for i = 1:npts
        if i == 1
            Du2=u[i+1]#-u[i]
            Du3=u[i+2]#-u[i]
             h2=r[i+1]#-r[i]
             h3=r[i+2]#-r[i]
            ur[i]=Du2/h2 + (Du3-Du2)/(h2-h3) + Du3/h3
        elseif i == npts
            Du1=u[i-1]-u[i]
            Du2=u[i-2]-u[i]
             h1=r[i-1]-r[i]
             h2=r[i-2]-r[i]
            ur[i]=Du1/h1 + (Du2-Du1)/(h1-h2) + Du2/h2
        else
            Du1=u[i-1]-u[i]
            Du3=u[i+1]-u[i]
             h1=r[i-1]-r[i]
             h3=r[i+1]-r[i]
            ur[i]=Du1/h1 + (Du3-Du1)/(h1-h3) + Du3/h3
        end
    end

    return ur
end

function durr(r,u)
    order = 2
    npts=length(r)
    urr=zero(u)
    Du=zeros(eltype(u),1+2*order)
    h =zeros(eltype(u),1+2*order)
    for i = 1:npts
        if i == 1
            urr[i]=0
        elseif i == npts
            Du[1]=u[i-1]-u[i]
            Du[2]=u[i-2]-u[i]
             h[1]=r[i-1]-r[i]
             h[2]=r[i-2]-r[i]
            urr[i]=(2*Du[1])/(h[1]*(h[1] - h[2])) - (2*Du[2])/((h[1] - h[2])*h[2])
        else
            Du[1]=u[i-1]-u[i]
            Du[3]=u[i+1]-u[i]
             h[1]=r[i-1]-r[i]
             h[3]=r[i+1]-r[i]
            urr[i]=(2*Du[1])/(h[1]*(h[1] - h[3])) - (2*Du[3])/((h[1] - h[3])*h[3])
        end
    end

    return urr
end

# we assume that u is a symmetric function, i.e. u(r)=u(-r)
function L(d,r,u;order=2, args...)
    npts = length(r)
    Lu=zero(u)
    if order == 1
        # this is a first order scheme with error err(Lu)=u'''/3*(h1+h2)
        # if h=h and h1=-h the scheme becomes second order
        for i = 1:npts
            if i == 1
                Lu[i]=NaN
            elseif i == npts
                Du1=u[i-2]-u[i]
                Du2=u[i-1]-u[i]
                 h1=r[i-2]-r[i]
                 h2=r[i-1]-r[i]
                Lu[i]=(Du2*((-1+d)h1-2r[i]))/((h1-h2)h2r[i])+(Du1*((1-d)*h2+2r[i]))/(h1*(h1-h2)*r[i])
            elseif i == 2
                h1=r[i+1]-r[i]
                h2=-r[i]
                D1=u[i+1]-u[i]
                D2=-u[i]
                Lu[i]=(-(D2*(2*r[i]+(1-d)*h1)/h2)+D1*(2*r[i]+(1-d)*h2)/h1)/(r[i]*(h1-h2))
            else
                h1=r[i+1]-r[i]
                h2=r[i-1]-r[i]
                D1=u[i+1]-u[i]
                D2=u[i-1]-u[i]
                @inbounds Lu[i]=(-(D2*(2*r[i]+(1-d)*h1)/h2)+D1*(2*r[i]+(1-d)*h2)/h1)/(r[i]*(h1-h2))
            end
        end
    elseif order == 2
        # this is a second order scheme
        for i = 1:npts
            if i == 1
                Lu[i]=NaN
            elseif i == 2
                Du1=u[i-1]-u[i]
                Du3=u[i+1]-u[i]
                Du4=u[i+2]-u[i]
                Du5=u[i+3]-u[i]
                h1 =r[i-1]-r[i]
                h3 =r[i+1]-r[i]
                h4 =r[i+2]-r[i]
                h5 =r[i+3]-r[i]
                Lu[i]=((-1+d)*Du5*h1*h3*h4)/((h1-h5)*h5*(-h3+h5)*(-h4+h5)*r[i])-(Du1*((-1+d)*h3*h4*h5+2(h3+h4)*(h1-h5)*r[i]))/(h1*(h1-h3)*(h1-h4)*(h1-h5)*r[i])+(Du3*((-1+d)*h1*h4*h5+2(h1+h4)*(h3-h5)*r[i]))/((h1-h3)*h3*(h3-h4)*(h3-h5)*r[i])+(Du4*((-1+d)*h1*h3*h5+2*(h1+h3)*(h4-h5)*r[i]))/((h1-h4)*h4*(-h3+h4)*(h4-h5)*r[i])
            elseif i == npts-1
                Du1=u[i-3]-u[i]
                Du2=u[i-2]-u[i]
                Du3=u[i-1]-u[i]
                Du5=u[i+1]-u[i]
                h1 =r[i-3]-r[i]
                h2 =r[i-2]-r[i]
                h3 =r[i-1]-r[i]
                h5 =r[i+1]-r[i]
                Lu[i]=-(((-1+d)Du1*h2*h3*h5)/(h1*(h1-h2)*(h1-h3)*(h1-h5)r[i]))+(Du5*((-1+d)h1*h2*h3-2(h2+h3)*(h1-h5)r[i]))/((h2-h5)h5*(-h1+h5)*(-h3+h5)r[i])+(Du3*((-1+d)h1*h2*h5-2(h1-h3)*(h2+h5)r[i]))/((h2-h3)h3*(-h1+h3)*(h3-h5)r[i])+(Du2*(-(-1+d)h1*h3*h5+2(h1-h2)*(h3+h5)r[i]))/(h2*(-h1+h2)*(h2-h3)*(h2-h5)r[i])
            elseif i == npts
                Du1=u[i-4]-u[i]
                Du2=u[i-3]-u[i]
                Du3=u[i-2]-u[i]
                Du4=u[i-1]-u[i]
                h1 =r[i-4]-r[i]
                h2 =r[i-3]-r[i]
                h3 =r[i-2]-r[i]
                h4 =r[i-1]-r[i]
                Lu[i]=-(((-1+d)Du1*h2*h3*h4)/(h1*(h1-h2)*(h1-h3)*(h1-h4)r[i]))+(Du4*((-1+d)h1*h2*h3-2(h2+h3)*(h1-h4)r[i]))/((h2-h4)h4*(-h1+h4)*(-h3+h4)r[i])+(Du3*((-1+d)h1*h2*h4-2(h1-h3)*(h2+h4)r[i]))/((h2-h3)h3*(-h1+h3)*(h3-h4)r[i])+(Du2*(-(-1+d)h1*h3*h4+2(h1-h2)*(h3+h4)r[i]))/(h2*(-h1+h2)*(h2-h3)*(h2-h4)r[i])
            else
                Du1=u[i-2]-u[i]
                Du2=u[i-1]-u[i]
                Du4=u[i+1]-u[i]
                Du5=u[i+2]-u[i]
                h1 =r[i-2]-r[i]
                h2 =r[i-1]-r[i]
                h4 =r[i+1]-r[i]
                h5 =r[i+2]-r[i]
                @inbounds Lu[i]=((d-1)*Du5*h1*h2*h4)/((h1-h5)*h5*(-h2+h5)*(-h4+h5)*r[i])-(Du1*((d-1)*h2*h4*h5+2(h2+h4)*(h1-h5)*r[i]))/(h1*(h1-h2)*(h1-h4)*(h1-h5)r[i])+(Du2*((d-1)*h1*h4*h5+2(h1+h4)*(h2-h5)*r[i]))/((h1-h2)*h2*(h2-h4)*(h2-h5)*r[i])+(Du4*((d-1)h1*h2*h5+2(h1+h2)*(h4-h5)r[i]))/((h1-h4)h4*(-h2+h4)*(h4-h5)r[i])
            end
        end
    end
    return Lu
end

# approximates urr+(d-1)(ur/r-u/r^2)
function Lu2(d,r,u)
    npts = length(r)
    Lu=zero(u)
    for i = 1:npts
        if i == 1
            Lu[i]=0
        elseif i == 2
            Lu[i]=(2^d*(-(u[i]/r[i])+((r[i]+r[1+i])^(1-d)*(-r[i]^(-1+d)*u[i]+r[1+i]^(-1+d)*u[1+i]))/(-r[i]+r[1+i])))/r[1+i]
        elseif i == npts
            Lu[i]=(2^d*(-(((r[-2+i]+r[-1+i])^(1-d)*(-r[-2+i]^(-1+d)*u[-2+i]+r[-1+i]^(-1+d)*u[-1+i]))/(-r[-2+i]+r[-1+i]))+((r[-1+i]+r[i])^(1-d)*(-r[-1+i]^(-1+d)*u[-1+i]+r[i]^(-1+d)*u[i]))/(-r[-1+i]+r[i])))/(-r[-2+i]+r[i])
        else
            Lu[i]=(2^d*(-(((r[-1+i]+r[i])^(1-d)*(-r[-1+i]^(-1+d)*u[-1+i]+r[i]^(-1+d)*u[i]))/(-r[-1+i]+r[i]))+((r[i]+r[1+i])^(1-d)*(-r[i]^(-1+d)*u[i]+r[1+i]^(-1+d)*u[1+i]))/(-r[i]+r[1+i])))/(-r[-1+i]+r[1+i])
        end
    end
    return Lu
end
