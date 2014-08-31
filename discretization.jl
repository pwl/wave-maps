# dur and L assume:
# - r[1]=0
# - the function u is symmetric around r=0, i.e. around r[1]

function dur(r,u)
    ur=zero(u)
    npts=length(r)
    for i = 1:npts
        if i == 1
            ur[i]=u[2]/r[2]
        elseif i == 2
            ur[i]=-(-2*r[2]*r[3]*u[2] + r[3]^2*u[2] + r[2]^2*u[3])/(r[2]^2*r[3] - r[2]*r[3]^2)
        elseif i == npts
            ur[i]=(r[i]^2*(-u[-2+i]+u[-1+i])+2r[-1+i]*r[i]*(u[-2+i]-u[i])+r[-2+i]^2*(u[-1+i]-u[i])+r[-1+i]^2*(-u[-2+i]+u[i])+2r[-2+i]*r[i]*(-u[-1+i]+u[i]))/((r[-2+i]-r[-1+i])*(r[-2+i]-r[i])*(r[-1+i]-r[i]))
        else
            ur[i]=(r[1+i]^2*(u[-1+i]-u[i])-2r[i]*(r[1+i]*(u[-1+i]-u[i])+r[-1+i]*(u[i]-u[1+i]))+r[i]^2*(u[-1+i]-u[1+i])+r[-1+i]^2*(u[i]-u[1+i]))/((r[-1+i]-r[i])*(r[-1+i]-r[1+i])*(r[i]-r[1+i]))
        end
    end

    return ur
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
                Lu[i]=(2(u[-2+i]/(r[-2+i]-r[i])+u[-1+i]/(-r[-1+i]+r[i])))/(r[-2+i]-r[-1+i])+((-1+d)*(r[-1+i]-r[i]))/(u[-1+i]-u[i])+(2u[i])/((r[-2+i]-r[i])*(r[-1+i]-r[i]))
            else
                h1=r[i+1]-r[i]
                h2=r[i-1]-r[i]
                D1=u[i+1]-u[i]
                D2=u[i-1]-u[i]
                Lu[i]=-((D2*(2*r[i]+(1-d)*h1))/(r[i]*(h1-h2)*h2))+(D1*(2*r[i]+(1-d)*h2))/(r[i]*h1*(h1-h2))
                # Lu[i]=((h1-h2)*(2r-(d-1)*(h1+h2))*u[r]+h2*(2r+(1-d)*h2)*u[r+h1]-h1*(2r+(1-d)*h1)*u[r+h2])/(r*h1*(h1-h2)*h2)
            end
        end
    elseif order == 2
        # this is a second order scheme
        for i = 1:npts
            if i == 1
                Lu[i]=NaN
            elseif i == 2
                Lu[i]=-(1/((r[2]-r[3])^2*r[3]^2))*r[2]^(-2-d)*(r[2]^d*(2r[2]-r[3])*r[3]*(-2r[2]r[3]u[2]+r[3]^2*u[2]+r[2]^2*u[3])+(r[2]^4*r[3]^d*(r[4]^2*(u[2]-u[3])-2r[3]*(r[4]*(u[2]-u[3])+r[2]*(u[3]-u[4]))+r[3]^2*(u[2]-u[4])+r[2]^2*(u[3]-u[4])))/((r[2]-r[4])*(r[3]-r[4])))
                #Lu[i]=(r[2]^(-2-d)*(r[2]^d*r[3]*(-2r[2]+r[3])*(r[2]-r[4])*(r[3]-r[4])*(-2r[2]*r[3]*u[2]+r[3]^2*u[2]+r[2]^2*u[3])-r[2]^4*r[3]^d*(r[4]^2*(u[2]-u[3])+r[3]^2*(u[2]-u[4])+r[2]^2*(u[3]-u[4])+2r[3]*(r[4]*(-u[2]+u[3])+r[2]*(-u[3]+u[4])))))/((r[2]-r[3])^2*r[3]^2*(r[2]-r[4])*(r[3]-r[4]))
            elseif i == 3
                Lu[i]=-((r[3]^(1-d)*((r[2]^(-2+d)*(r[3]-r[4])^3*(-2r[2]r[3]u[2]+r[3]^2u[2]+r[2]^2u[3]))/r[3]-r[3]^(-1+d)*(r[2]-2r[3]+r[4])*(r[4]^2*(u[2]-u[3])+r[3]^2*(u[2]-u[4])+r[2]^2*(u[3]-u[4])+2r[3]*(r[4]*(-u[2]+u[3])+r[2]*(-u[3]+u[4])))+(1/((r[3]-r[5])*(r[4]-r[5])))*(r[2]-r[3])^3*r[4]^(-1+d)*(r[5]^2*(u[3]-u[4])+r[4]^2*(u[3]-u[5])+r[3]^2*(u[4]-u[5])+2r[4]*(r[5]*(-u[3]+u[4])+r[3]*(-u[4]+u[5])))))/((r[2]-r[3])^2*(r[2]-r[4])*(r[3]-r[4])^2))
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
