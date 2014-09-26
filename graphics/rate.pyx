#   1    2   3   4     5     6     7
# [tau] [t] [s] [ur] [urt] [rate] [T]
#  1   2    3      4    5   6    7    9
# [t] [s] [tau] [rate] [u] [ur] [ut] [urt]
load "files.pyx"

# set output "rate.pdf"
# set terminal pdf

# fitting
load "fit.pyx"

# plotting
set nodisplay

set logscale y
set xrange [-.3:10]
set yrange [1e-3:100]
set title "Rate of convergence to $\phi_0$"
set title ""
set xlabel "$-\log(T-t)$"
set xlabel "$s$"
# set yformat "$10^{%2.0f}$"%(log(y))
# set yticks 1e-4,10,1e-1

plot at0 u 2:(abs($4)) w p pt 3 ps .4 c Black t "$\partial_y U(0,s)-\phi_0'(0)$"
replot exp(flin(x)) w l lt 2 c BrickRed t "$Ce^{\lambda_1 s}$"
replot 0 w c White t "$\lambda_1=%5.3f$"%(L1)
replot 0 w c White t "$C=%5.3f$"%(-exp(C))

set display
refresh
