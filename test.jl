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
        p = summarize(tau,t,r,u,ur,d)
        savefig(p,"$(string(name))-$(val).png", "width", 1024, "height", 1024)
    end
end
