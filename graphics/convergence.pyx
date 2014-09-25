#  1   2   3   4    5     6     7    8   9
# [r] [y] [u] [ur] [ut] [urt] [tau] [t] [T]
soln="wm/d4-npts200-convergence/soln.dat"

set output "convergence.pdf"
set terminal pdf

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
set xticks
unset xlabel
unset ylabel
set nomxticks
}

subroutine markson()
{
set xformat auto
set yformat auto
set xticks (1e-1 "$10^{-1}$",1 "$1$",10 "$10$", 100 "$10^2$")
set yticks (0 "$0$",pi "$\pi$")
set xlabel "$y$"
set ylabel "$U(s,y)$"
}

smin=0
smax=10

set multiplot
sz = 5
rt = .6
set size sz ratio rt
rows = 3
cols = 2
si = [0.0, 1.3, 2.0, 4.0, 6.3, 10.0]

set logscale x
set xrange [2e-2:1e3]
set yrange [-0.5*pi:1.4*pi]
set nokey

set nodisplay

for i = 0 to cols-1 {
for j = 0 to rows-1 {

set origin sz*i,sz*rt*j
if (i==0 && j==0) {markson();} else {marksoff();}

t=(i/cols+(rows-1)-j)/((cols-1)/cols+(rows-1))
s=smin+t*(smax-smin)
s=si[(i+cols*(rows-1-j))]
sstr="%i"%(ind(s))

labeltext="$s=%2.1f$"%(s)
set label 1 labeltext at 0.2,3.5

plot soln u 2:3 w l every ::1 index @sstr
replot 2*atan(x/sqrt(2)) w l lt 2 t ""
replot 0 w lt 2 c Gray, pi w lt 2 c Gray
}
}

set display
refresh
