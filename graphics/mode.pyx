#  1   2    3      4    5   6     7    9    10     11
# [t] [s] [tau] [rate] [u] [ur] [urr] [ut] [utr] [utrr]
#  1   2    3    4   5   6   7     8    9    10     11
# [t] [s] [tau] [r] [y] [u] [ur] [urr] [ut] [utr] [utrr]
load "files.pyx"

# set output "mode.pdf"
# set terminal pdf

d=4
phi0(y)=2*atan(y/sqrt(d-2))
s0=4
cmd="awk 'BEGIN{i=0}; /#s=/{i++;if($2>%f){print i; exit 0}};' %s"%(s0,soln)
ind=os.popen(cmd,"r").read()+0

# fitting
load "fit.pyx"

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
plot soln u 5:($6-phi0($5)) i ind w p pt 3 ps .4 c Black every 5 t "$U(y,s)-\phi_0\left(y\right)$"
replot mode u 1:(-$2*exp(flin(s0))) w l lt 2 c BrickRed lw 1.2 t "$Ce^{\lambda_1 s}v_1\left(y\right)$"

set display
refresh
