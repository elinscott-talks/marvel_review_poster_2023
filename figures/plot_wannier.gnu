unset raxis
unset border
unset ytics 
unset xtics
set sample 1000
set object circle at -9,0 size 0.4 fillcolor rgb "red" fillstyle solid border lt 1 lc rgb "black" lw 8 front
set object circle at -6,0 size 0.4 fillcolor rgb "red" fillstyle solid border lt 1 lc rgb "black" lw 8 front
set object circle at -3,0 size 0.4 fillcolor rgb "red" fillstyle solid border lt 1 lc rgb "black" lw 8 front
set object circle at  0,0 size 0.4 fillcolor rgb "red" fillstyle solid border lt 1 lc rgb "black" lw 8 front
set object circle at  3,0 size 0.4 fillcolor rgb "red" fillstyle solid border lt 1 lc rgb "black" lw 8 front
set object circle at  6,0 size 0.4 fillcolor rgb "red" fillstyle solid border lt 1 lc rgb "black" lw 8 front
set object circle at  9,0 size 0.4 fillcolor rgb "red" fillstyle solid border lt 1 lc rgb "black" lw 8 front
set object circle at 12,0 size 0.4 fillcolor rgb "red" fillstyle solid border lt 1 lc rgb "black" lw 8 front
set object circle at 15,0 size 0.4 fillcolor rgb "red" fillstyle solid border lt 1 lc rgb "black" lw 8 front

set arrow nohead from 0,-0.3 to 0, 0.3 lw 8 
a=3
R0=3
b=0.15
f(x)=exp(-a*x**2)*sin(x)+exp(-a*(x+R0)**2)*sin(x+R0)+exp(-a*(x-R0)**2)*sin(x-R0)+exp(-a*(x+2*R0)**2)*sin(x+2*R0)+exp(-a*(x-2*R0)**2)*sin(x-2*R0)+exp(-a*(x+3*R0)**2)*sin(x+3*R0)+exp(-a*(x-3*R0)**2)*sin(x-3*R0)+exp(-a*(x+4*R0)**2)*sin(x+4*R0)+exp(-a*(x-4*R0)**2)*sin(x-4*R0)+exp(-a*(x+5*R0)**2)*sin(x+5*R0)+exp(-a*(x-5*R0)**2)*sin(x-5*R0)+exp(-a*(x+6*R0)**2)*sin(x+6*R0)+exp(-a*(x-6*R0)**2)*sin(x-6*R0)
pl [-5:16.5]f(x)*cos(b*x) notit w l lw 8 lc rgb "blue" , 0.243*cos(b*x) notit w l lw 8 lc rgb "orange", 0 notit w l lc rgb "black" lw 8 
set term pngcairo enha color size 1200,600
set out "bloch.png"
repl

set term qt 
pl [-5:16.5] exp(-a*(x-R0)**2)*sin(x-R0) notit w l lw 8 lc rgb "blue", 0 notit w l lw 8 lc rgb "black" 

set term pngcairo enha color size 1200,600
set out "wannier.png"
repl
