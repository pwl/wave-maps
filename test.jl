include("plot.jl")

macro testparameter(name,range,fcall)
    println("Testing parameter $(string(name))")
    for val=eval(range)
        println()
        println("$(string(name)) = $val")
        newcall = fcall
        push!(newcall.args,Expr(:kw,name,val))
        @eval (tau,t,r,u,ur)=$newcall

        # plot the mesh
        p=plotmesh(tau,r,xrange=[1e-10,pi],yrange=[0,30])
        savefig(p,"mesh-$(string(name))-$(val).png")

        # plot the function
        T=t[end]+2/ur[end][1,1]       # crude estimate on blow-up time
        println(-log(T-t[end]))
        p=plotmode(r,u,ur,-log(T-t),6,d=4) # @todo: fixed dimension
        savefig(p,"mode-$(string(name))-$(val).png")
        p=plotu(r,u,ur,-log(T-t),6,d=4,xlog=true) # @todo: fixed dimension
        savefig(p,"u-$(string(name))-$(val).png")
    end
end
