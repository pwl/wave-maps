include("plot.jl")

const s0 = 4
# const d  = 4

macro testparameter(name,range,fcall)
    println("Testing parameter $(string(name))")
    for val=eval(range)
        println()
        println("$(string(name)) = $val")
        newcall = fcall
        push!(newcall.args,Expr(:kw,name,val))
        @eval (tau,t,r,u,ur)=$newcall
        summarize(tau,t,r,u,ur,"$(string(name))-$(val).png")
    end
end

function summarize(tau,t,r,u,ur,fn)
    T=t[end]+2/ur[end][1,1]       # crude estimate on blow-up time
    s=-log(T-t)

    # plot the mesh
    # p=plotmesh(tau,r,xrange=[1e-10,pi],yrange=[0,10])
    p1=plotmesh(s,r,xrange=[1e-8,pi],yrange=[0,20])

    # plot the eigenmode
    p2=plotmode(r,u,ur,s,s0,d=d) # @todo: fixed dimension

    # plot the self-similar profile
    p3=plotu(r,u,ur,s,s0,d=d,xlog=true,xrange=[1e-2,100],yrange=[0,pi]) # @todo: fixed dimension

    # plot the blow-up
    p4=plot(s,abs(map(z->z[1,1],ur).*exp(-s).-sqrt(d-2)),title="u'(0,t)-sqrt(d-2)",ylog=true)

    tab=Table(2,2)
    tab[1,1]=p1
    tab[1,2]=p2
    tab[2,1]=p3
    tab[2,2]=p4
    savefig(tab,fn, "width", 1024, "height", 1024)
end
