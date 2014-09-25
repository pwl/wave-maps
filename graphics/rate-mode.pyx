set output "rate-mode.pdf"
set terminal pdf

set multiplot
sz = 6.5
set size sz ratio 1
set nodisplay

set origin 0,0
load "rate.pyx"

set axis y right
set origin sz+1,0
load "mode.pyx"

unset multiplot
set display
refresh
