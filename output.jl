function outputat0line(f,tau,t,ur,T,d)
    s = -log(T-t)
    rate = ur[1]*exp(-s)-2/sqrt(d-2)

    write(f,"$tau $t $s ")
    for i = 1:length(ur)
        write(f, "$(ur[i]) ")
    end
    write(f,"$rate $T\n")
end

function outputat0(tau,t,ur,d,T,fn;every=1)
    f = open(fn,"w")
    write(f,"# [tau] [t] [s] [ur] [urt] [rate] [T]\n")
    for i = 1:every:length(t)
        outputat0line(f,tau[i],t[i],ur[i][1,:],T,d)
    end
    close(f)
end

function outputsolline(tau,t,r,u,ur,T,f)
    y = r/(T-t)
    write(f,"$r $y ")
    for i = 1:length(u)
        uval = u[i]
        urval = ur[i]
        write(f, "$uval $urval ")
    end
    write(f,"$tau $t $T")
    write(f,"\n")
end

function outputsolblock(tau,t,r,u,ur,T,f)
    write(f,"#t= $t\n")
    write(f,"#s= $(-log(T-t))\n")
    for i = 1:length(r)
        outputsolline(tau,t,r[i],u[i,:],ur[i,:],T,f)
    end
    write(f,"\n\n\n")
end

function outputsol(tau,t,r,u,ur,T,fn;every=1)
    f = open(fn,"w")
    write(f,"# [r] [y] [u] [ur] [ut] [urt] [tau] [t] [T]")
    for i = 1:every:length(t)
        outputsolblock(tau[i],t[i],r[i],u[i],ur[i],T,f)
    end
    close(f)
end
