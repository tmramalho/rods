set term png
set output "vars.png"
set xlabel "Monte carlo steps"
set ylabel "Energy"
plot "energy - T0.2rho0.2.dat" using 0:1, "vars.dat" using 0:2

set output "ac0.png"
set xlabel "Monte carlo steps"
set ylabel "X(t)"
#set logscale y
set xrange [0:40]
f1(x) = a * exp(- x / t);
a = 1
t = 0.5
fit f1(x) "ac0.dat" using 0:1 via t
plot "ac0.dat" using 0:1, a * exp ( -x / t )

set output "ac1.png"
set xlabel "Monte carlo steps"
set ylabel "X(t)"
#set logscale y
set xrange [0:40]
f1(x) = a * exp(- x / t);
a = 1
t = 0.5
fit f1(x) "ac1.dat" using 0:1 via t
plot "ac1.dat" using 0:1, a * exp ( -x / t )

