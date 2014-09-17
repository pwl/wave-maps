#  1   2   3   4    5     6     7    8   9
# [r] [y] [u] [ur] [ut] [urt] [tau] [t] [T]
soln="wm/d4-npts200-convergence/soln.dat"

set output "convergence.eps"
set terminal eps

subroutine ind(s0)
{
cmd="awk 'BEGIN{i=0}; /#s=/{i++;if($2>%f){print i; exit 0}};' %s"%(s0,soln)
ind=os.popen(cmd,"r").read()+0
return(ind)
}

subroutine marksoff()
{
set noyticks
set xformat ""
unset xlabel
unset ylabel
set nomxticks
}

subroutine markson()
{
set xformat auto
set yformat auto
set xticks
set yticks (0 "$0$",pi "$\pi$")
set xlabel "$\frac{r}{T-t}$"
set ylabel "$u(r,t)$"
}

smin=0
smax=10

set multiplot
sz = 5
set size sz ratio 1
rows = 3
cols = 3

set logscale x
set xrange [1e-2:1e3]
set yrange [-0.5*pi:1.4*pi]
set nokey
set axis x top

set nodisplay

for i = 0 to cols-1 {
for j = 0 to rows-1 {

set origin sz*i,sz*j
if (i==0 && j==rows-1) {markson();} else {marksoff();}

t=(i/cols+(rows-j-1))/rows
s=smin+t*(smax-smin)
sstr="%i"%(ind(s))

power=floor(-s/log(10))
if (power < 0) {labeltext = "$T-t=%2.1f\cdot 10^{%i}$"%(exp(-s)/10**power,power);
} else {labeltext = "$T-t=%2.1f$"%(exp(-s));}
set label 1 labeltext at 0.02,3.5

plot soln u 2:3 w l every ::2 index @sstr
replot 2*atan(x/sqrt(2)) w l lt 2 t ""
replot 0 w lt 2 c Gray, pi w lt 2 c Gray
}
}

set display
refresh
