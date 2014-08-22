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
                Lu[i]=(r[2]^(-2-d)*(r[2]^d*r[3]*(-2r[2]+r[3])*(r[2]-r[4])*(r[3]-r[4])*(-2r[2]*r[3]*u[2]+r[3]^2*u[2]+r[2]^2*u[3])-r[2]^4*r[3]^d*(r[4]^2*(u[2]-u[3])+r[3]^2*(u[2]-u[4])+r[2]^2*(u[3]-u[4])+2r[3]*(r[4]*(-u[2]+u[3])+r[2]*(-u[3]+u[4])))))/((r[2]-r[3])^2*r[3]^2*(r[2]-r[4])*(r[3]-r[4]))
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
