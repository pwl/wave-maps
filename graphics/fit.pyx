# fitting
C=0.01
L1=-0.5
D=20
L2=-1.5
f(x)=C+L1*x+D*exp(L2*x)
flin(x)=C+L1*x
smax = 6
fit [4:smax] flin() withouterrors @at0 u 2:(log(abs($4))) via C,L1
fit [3:4] f() withouterrors @at0 u 2:(log(abs($4))) via D,L2
fit [3:smax] f() withouterrors @at0 u 2:(log(abs($4))) via C,L1,D,L2
