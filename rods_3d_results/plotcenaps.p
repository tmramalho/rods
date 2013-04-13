set terminal postscript
set xlabel "Temperature"
set ylabel "c"
set output "c.ps"
plot "./mineresult.dat" u 1:2:3 with yerrorbars notitle, "./mineresult.dat" u 1:2 with lines notitle
set ylabel "d1"
set output "d1.ps"
plot "./mineresult.dat" u 1:4:5 with yerrorbars notitle, "./mineresult.dat" u 1:4 with lines notitle
set ylabel "d1f"
set output "d1f.ps"
plot "./mineresult.dat" u 1:6:7 with yerrorbars notitle, "./mineresult.dat" u 1:6 with lines notitle
set ylabel "d2"
set output "d2.ps"
plot "./mineresult.dat" u 1:8:9 with yerrorbars notitle, "./mineresult.dat" u 1:8 with lines notitle
set ylabel "d2f"
set output "d2f.ps"
plot "./mineresult.dat" u 1:10:11 with yerrorbars notitle, "./mineresult.dat" u 1:10 with lines notitle
set ylabel "l"
set output "l.ps"
plot "./mineresult.dat" u 1:12:13 with yerrorbars notitle, "./mineresult.dat" u 1:12 with lines notitle
set ylabel "lf"
set output "lf.ps"
plot "./mineresult.dat" u 1:14:15 with yerrorbars notitle, "./mineresult.dat" u 1:14 with lines notitle
