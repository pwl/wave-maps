include("wave-mov.jl")

# parameters leading to good results.  maxorder, and tolerances are critical, solution converges with increasing npts
(tau,t,r,u,ur)=wavesolve(W<,r->r+sin(r),taumax=25,npts=50,passes=8,maxorder=3,abstol=1e-6,reltol=1e-6)
T=t[end]+2/ur[end][1,1]
s=-log(T-t)

p=plotu(r,u,ur,s,6,d=4,xlog=true)
display(p)
