type Results{R<:Real}
    tau :: Array{R,1}
    t   :: Array{R,1}
    r   :: Array{Vector{R},1}
    u   :: Array{Matrix{R},1}
    ur  :: Array{Matrix{R},1}
    urr :: Array{Matrix{R},1}

    ### These fields are filled at the end of the simulation by finalize!
    #
    # number of mesh points
    npts :: Int
    # estimate on the blow-up time
    T   :: R
    # slow time
    s   :: Array{R,1}
    # self-similar spatial variable
    y   :: Array{Vector{R},1}
    #########################################################

    function Results()
        this=new()
        this.tau = Array(R,0)
        this.t   = Array(R,0)
        this.r   = Array(Vector{R},0)
        this.u   = Array(Matrix{R},0)
        this.ur  = Array(Matrix{R},0)
        this.urr = Array(Matrix{R},0)
        this.y   = Array(Matrix{R},0)
        return this
    end
end

function push{T<:Real}(res :: Results{T},
                       tau :: T,
                       t   :: T,
                       r   :: Vector{T},
                       u   :: Matrix{T})
    ur  = zero(u)
    urr = zero(u)

    for j=1:size(u,2)
        dur!(view(ur ,:,j),r,view(u ,:,j))
        dur!(view(urr,:,j),r,view(ur,:,j))
    end

    push!(res.tau, tau )
    push!(res.t,   t   )
    push!(res.r,   r   )
    push!(res.u,   u   )
    push!(res.ur,  ur  )
    push!(res.urr, urr )
end

# crude estimate on blow-up time
function blowuptime(res::Results,eqn::Equation)
    T=res.t[end]
end

# the definition of self-similar variable y
function selfsimilary(r,s,eqn::Equation)
    return exp(s)*r
end

function finalize!(res::Results,eqn::Equation)
    res.T = blowuptime(res,eqn)
    res.s = -log(res.T.-res.t)
    res.npts = length(res.r[1])
    for i = 1:length(res.t)
        push!(res.y,selfsimilary(res.r[i],res.s[i],eqn))
    end
end
