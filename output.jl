# this method saves the values of solution at r=0
function saveat0(res::Results,fn::String)
    open(fn,"w") do f
        write(f,"[t] [s] [tau] [u] [ur] [urr] ... [uN] [urN] [urrN]\n")
        for line = 1:length(res.t)
            saveat0line(f,
                        res.t[line],
                        res.s[line],
                        res.tau[line],
                        vec(res.u[line][1,:]),
                        vec(res.ur[line][1,:]),
                        vec(res.urr[line][1,:]))
        end
    end
end

function saveat0line(f,t,s,tau,u,ur,urr)
    write(f,"$t $s $tau ")
    for i = 1:length(u)
        write(f,"$(u[i]) $(ur[i]) $(urr[i])")
    end
    write(f,"\n")
end
