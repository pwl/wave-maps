#  1   2   3   4    5     6     7    8   9
# [r] [y] [u] [ur] [ut] [urt] [tau] [t] [T]
at0="wm/d4-npts800/at0.dat"
soln="wm/d4-npts800/soln.dat"
mode="../modes/mode_d4_n1.dat"

# set output "mode.pdf"
# set terminal pdf

d=4
phi0(y)=2*atan(y/sqrt(d-2))
s0=3
cmd="awk 'BEGIN{i=0}; /#s=/{i++;if($2>%f){print i; exit 0}};' %s"%(s0,soln)
ind=os.popen(cmd,"r").read()+0

# fitting
c=-2.7
L1=-0.5
D=5.
L2=-1.5
f(x)=C+L1*x+D*exp(L2*x)
flin(x)=C+L1*x
smax = 6
fit [3:smax] flin() withouterrors @at0 u 3:(log(abs($6))) via C,L1
fit [2:3] f() withouterrors @at0 u 3:(log(abs($6))) via D,L2
fit [2:smax] f() withouterrors @at0 u 3:(log(abs($6))) via C,L1,D,L2

# plotting
set nodisplay

set xrange [0:2]
set nolog
set yrange [*:*]
set yformat auto
set title "First mode"
set title ""
set xlabel "$y$"
set key top left

# plot soln u 2:3 w l t ""
# replot soln u 2:(phi0($2)) w l t ""
plot soln u 2:($3-phi0($2)) i ind w p pt 3 ps .4 c Black every 10 t "$U(y,s)-\phi_0\left(y\right)$"
replot mode u 1:(-$2*exp(flin(s0))) w l lt 2 c BrickRed lw 1.2 t "$Ce^{\lambda_1 s}v_1\left(y\right)$"

set display
refresh
