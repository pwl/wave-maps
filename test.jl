include("plot.jl")

macro testparameter(name,range,fcall)
    println("Testing parameter $(string(name))")
    for val=eval(range)
        println()
        println("$(string(name)) = $val")
        newcall = fcall
        push!(newcall.args,Expr(:kw,name,val))
        @eval (tau,t,r,u,ur)=$newcall
        savefig(plotmesh(tau,r),"$(string(name))-$(val).png")
    end
end
