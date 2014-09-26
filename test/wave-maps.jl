include("../src/wave-mov.jl")

# set the dimension
dim = 4
# define the equation
eqn=WMEquation(dim);
# define the initial data
u0(r)=-2sin(4r).*exp(-(r.-1).^2)

# compute the numerical solution (to get a reasonable mode profile you
# should consider running this function with npts=800 and
# tolerances=1e-8
res=wavesolve(eqn,u0,npts=200,reltol=1e-6,abstol=1e-6)

# summarize the results at s=6
display(summarize(res,eqn,6))

save(res,"soln.dat",every=10)   # save every 10-th step of the
                                # numerical solution to "soln.dat"
saveat0(res,eqn,"at0.dat",every=100)   # save every 100-th step of the
                                       # derivatives at r=0 to
                                       # "at0.dat"
