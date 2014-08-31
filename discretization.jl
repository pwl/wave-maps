# dur and L assume:
# - r[1]=0
# - the function u is symmetric around r=0, i.e. around r[1]

function dur(r,u)
    order = 2
    npts=length(r)
    ur=zero(u)
    Du=zeros(eltype(u),1+2*order)
    h =zeros(eltype(u),1+2*order)
    for i = 1:npts
        if i == 1
            Du[2]=u[i+1]#-u[i]
            Du[3]=u[i+2]#-u[i]
             h[2]=r[i+1]#-r[i]
             h[3]=r[i+2]#-r[i]
            ur[i]=Du[2]/h[2] + (Du[3]-Du[2])/(h[2]-h[3]) + Du[3]/h[3]
        elseif i == 2
            Du[1]=-u[i]
            Du[3]=u[i+1]-u[i]
             h[1]=-r[i]
             h[3]=r[i+1]-r[i]
            ur[i]=Du[1]/h[1] + (Du[3]-Du[1])/(h[1]-h[3]) + Du[3]/h[3]
        elseif i == npts
            Du[1]=u[i-1]-u[i]
            Du[2]=u[i-2]-u[i]
             h[1]=r[i-1]-r[i]
             h[2]=r[i-2]-r[i]
            ur[i]=Du[1]/h[1] + (Du[2]-Du[1])/(h[1]-h[2]) + Du[2]/h[2]
        else
            Du[1]=u[i-1]-u[i]
            Du[3]=u[i+1]-u[i]
             h[1]=r[i-1]-r[i]
             h[3]=r[i+1]-r[i]
            ur[i]=Du[1]/h[1] + (Du[3]-Du[1])/(h[1]-h[3]) + Du[3]/h[3]
        end
    end

    return ur
end

# we assume that u is a symmetric function, i.e. u(r)=u(-r)
function L(d,r,u;order=2, args...)
    npts = length(r)
    Lu=zero(u)
    Du=zeros(eltype(u),1+2*order)
    h =zeros(eltype(u),1+2*order)
    if order == 1
        # this is a first order scheme with error err(Lu)=u'''/3*(h1+h2)
        # if h=h and h1=-h the scheme becomes second order
        for i = 1:npts
            if i == 1
                Lu[i]=NaN
            elseif i == npts
                Du[1]=u[i-2]-u[i]
                Du[2]=u[i-1]-u[i]
                 h[1]=r[i-2]-r[i]
                 h[2]=r[i-1]-r[i]
                Lu[i]=(Du[2]*((-1+d)h[1]-2r[i]))/((h[1]-h[2])h[2]r[i])+(Du[1]*((1-d)*h[2]+2r[i]))/(h[1]*(h[1]-h[2])*r[i])
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
                Lu[i]=(-(D2*(2*r[i]+(1-d)*h1)/h2)+D1*(2*r[i]+(1-d)*h2)/h1)/(r[i]*(h1-h2))
            end
        end
    elseif order == 2
        # this is a second order scheme
        for i = 1:npts
            if i == 1
                Lu[i]=NaN
            elseif i == 2
                Du[1]=u[i-1]-u[i]
                Du[3]=u[i+1]-u[i]
                Du[4]=u[i+2]-u[i]
                Du[5]=u[i+3]-u[i]
                h[1] =r[i-1]-r[i]
                h[3] =r[i+1]-r[i]
                h[4] =r[i+2]-r[i]
                h[5] =r[i+3]-r[i]
                Lu[i]=((-1+d)*Du[5]*h[1]*h[3]*h[4])/((h[1]-h[5])*h[5]*(-h[3]+h[5])*(-h[4]+h[5])*r[i])-(Du[1]*((-1+d)*h[3]*h[4]*h[5]+2(h[3]+h[4])*(h[1]-h[5])*r[i]))/(h[1]*(h[1]-h[3])*(h[1]-h[4])*(h[1]-h[5])*r[i])+(Du[3]*((-1+d)*h[1]*h[4]*h[5]+2(h[1]+h[4])*(h[3]-h[5])*r[i]))/((h[1]-h[3])*h[3]*(h[3]-h[4])*(h[3]-h[5])*r[i])+(Du[4]*((-1+d)*h[1]*h[3]*h[5]+2*(h[1]+h[3])*(h[4]-h[5])*r[i]))/((h[1]-h[4])*h[4]*(-h[3]+h[4])*(h[4]-h[5])*r[i])
            elseif i == npts-1
                Du[1]=u[i-3]-u[i]
                Du[2]=u[i-2]-u[i]
                Du[3]=u[i-1]-u[i]
                Du[5]=u[i+1]-u[i]
                h[1] =r[i-3]-r[i]
                h[2] =r[i-2]-r[i]
                h[3] =r[i-1]-r[i]
                h[5] =r[i+1]-r[i]
                Lu[i]=-(((-1+d)Du[1]h[2]h[3]h[5])/(h[1]*(h[1]-h[2])*(h[1]-h[3])*(h[1]-h[5])r[i]))+(Du[5]*((-1+d)h[1]h[2]h[3]-2(h[2]+h[3])*(h[1]-h[5])r[i]))/((h[2]-h[5])h[5]*(-h[1]+h[5])*(-h[3]+h[5])r[i])+(Du[3]*((-1+d)h[1]h[2]h[5]-2(h[1]-h[3])*(h[2]+h[5])r[i]))/((h[2]-h[3])h[3]*(-h[1]+h[3])*(h[3]-h[5])r[i])+(Du[2]*(-(-1+d)h[1]h[3]h[5]+2(h[1]-h[2])*(h[3]+h[5])r[i]))/(h[2]*(-h[1]+h[2])*(h[2]-h[3])*(h[2]-h[5])r[i])
            elseif i == npts
                Du[1]=u[i-4]-u[i]
                Du[2]=u[i-3]-u[i]
                Du[3]=u[i-2]-u[i]
                Du[4]=u[i-1]-u[i]
                h[1] =r[i-4]-r[i]
                h[2] =r[i-3]-r[i]
                h[3] =r[i-2]-r[i]
                h[4] =r[i-1]-r[i]
                Lu[i]=-(((-1+d)Du[1]h[2]h[3]h[4])/(h[1]*(h[1]-h[2])*(h[1]-h[3])*(h[1]-h[4])r[i]))+(Du[4]*((-1+d)h[1]h[2]h[3]-2(h[2]+h[3])*(h[1]-h[4])r[i]))/((h[2]-h[4])h[4]*(-h[1]+h[4])*(-h[3]+h[4])r[i])+(Du[3]*((-1+d)h[1]h[2]h[4]-2(h[1]-h[3])*(h[2]+h[4])r[i]))/((h[2]-h[3])h[3]*(-h[1]+h[3])*(h[3]-h[4])r[i])+(Du[2]*(-(-1+d)h[1]h[3]h[4]+2(h[1]-h[2])*(h[3]+h[4])r[i]))/(h[2]*(-h[1]+h[2])*(h[2]-h[3])*(h[2]-h[4])r[i])
            else
                Du[1]=u[i-2]-u[i]
                Du[2]=u[i-1]-u[i]
                Du[4]=u[i+1]-u[i]
                Du[5]=u[i+2]-u[i]
                h[1] =r[i-2]-r[i]
                h[2] =r[i-1]-r[i]
                h[4] =r[i+1]-r[i]
                h[5] =r[i+2]-r[i]
                Lu[i]=((d-1)*Du[5]*h[1]*h[2]*h[4])/((h[1]-h[5])*h[5]*(-h[2]+h[5])*(-h[4]+h[5])*r[i])-(Du[1]*((d-1)*h[2]*h[4]*h[5]+2(h[2]+h[4])*(h[1]-h[5])*r[i]))/(h[1]*(h[1]-h[2])*(h[1]-h[4])*(h[1]-h[5])r[i])+(Du[2]*((d-1)*h[1]*h[4]*h[5]+2(h[1]+h[4])*(h[2]-h[5])*r[i]))/((h[1]-h[2])*h[2]*(h[2]-h[4])*(h[2]-h[5])*r[i])+(Du[4]*((d-1)h[1]*h[2]*h[5]+2(h[1]+h[2])*(h[4]-h[5])r[i]))/((h[1]-h[4])h[4]*(-h[2]+h[4])*(h[4]-h[5])r[i])
            end
        end
    end
    return Lu
end
