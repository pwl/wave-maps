#   1    2   3   4     5     6     7
# [tau] [t] [s] [ur] [urt] [rate] [T]
at0="wm/d4-npts800/at0.dat"

set output "rate.pdf"
set terminal pdf

# fitting
C=0.01
L1=-0.5
d=20
L2=-1.5
f(x)=C+L1*x+d*exp(L2*x)
flin(x)=C+L1*x
smax = 6
fit [3:smax] flin() withouterrors @at0 u 3:(log(abs($6))) via C,L1
fit [2:3] f() withouterrors @at0 u 3:(log(abs($6))) via d,L2
fit [2:smax] f() withouterrors @at0 u 3:(log(abs($6))) via C,L1,d,L2

# plotting
set nodisplay

set logscale y
set xrange [-.3:]
set yrange [1e-4:1]
set title "Rate of convergence to $\phi_0$"
set title ""
set xlabel "$-\log(T-t)$"
set xlabel "$s$"
set yformat "$10^{%2.0f}$"%(log(y))
set yticks 1e-4,10,1e-1

# plot at0 u 3:(abs($6)) every 100 w p pt 3 ps .4 c Black t "$(T-t)\partial_r u(0,t)-\phi_0'(0)$"
# replot exp(flin(x)) w l lt 2 c BrickRed t "$C_1\cdot (T-t)^{-\lambda_1}$"
# replot 0 w c White t "$C_1=%5.3f$, $\lambda_1=%5.3f$"%(-exp(C),L1)

plot at0 u 3:(abs($6)) every 100 w p pt 3 ps .4 c Black t "$\partial_y U(0,s)-\phi_0'(0)$"
replot exp(flin(x)) w l lt 2 c BrickRed t "$Ce^{\lambda_1 s}$"
replot 0 w c White t "$\lambda_1=%5.3f$, $C=%5.3f$"%(L1,-exp(C))

set display
refresh
