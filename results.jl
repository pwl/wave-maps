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

function setT(res::Results,eqn::Equation,T)
    res.T = T
    res.s = -log(res.T.-res.t)
    for i = 1:length(res.t)
        res.y[i]=selfsimilary(res.r[i],res.s[i],eqn)
    end
end

function finalize!(res::Results,eqn::Equation)
    res.npts = length(res.r[1])
    res.y = copy(res.r)         # allocate the space for res.y
    setT(res,eqn,blowuptime(res,eqn))
end

function save(res::Results,fn::String;
              every=1)          # save only `every` step
    f = open(fn,"w")
    write(f,"# [tau] [t] [r] [y] [u_1] [ur_1] [urr_1] ... [u_N] [ur_N] [urr_N]")
    for line = 1:every:length(res.t)
        outputsolblock(f,res,line)
    end
    close(f)
end

function outputsolblock(f,res,line)
    write(f,"#t= $(res.t[line])\n")
    write(f,"#s= $(res.s[line])\n")
    for i = 1:length(res.r)
        outputsolline(f,
                      res.tau[line],
                      res.t[line],
                      res.r[i][line],
                      res.y[i][line],
                      res.u[i][:,line],
                      res.ur[i][:,line],
                      res.urr[i][:,line])
    end
    write(f,"\n\n\n")
end

function outputsolline(f,tau,t,r,y,u,ur,urr)
    write(f,"$tau $t ")
    write(f,"$r $y ")
    for i = 1:length(u)
        write(f, "$(u[i]) $(ur[i]) $(urr[i]) ")
    end
    write(f,"\n")
end
