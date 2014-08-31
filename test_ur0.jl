include("wave-mov.jl")
(tau,t,r,u,ur)=wavesolve(myRhs,r->r+sin(r),taumax=10,npts=50,gamma=1,passes=0)

T=t[end]+2/ur[end][1,1]

p=plotu(r,u,ur,-log(T-t),6,d=4,xlog=true)
display(p)
