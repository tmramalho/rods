set terminal png size 1600,1200
set xlabel "temp"
set ylabel "c"
set output "c.png"
plot "./mineresult.dat" u 1:2:3 with yerrorbars, "./mineresult.dat" u 1:2 with lines
set ylabel "d1"
set output "d1.png"
plot "./mineresult.dat" u 1:4:5 with yerrorbars, "./mineresult.dat" u 1:4 with lines
set ylabel "d1f"
set output "d1f.png"
plot "./mineresult.dat" u 1:6:7 with yerrorbars, "./mineresult.dat" u 1:6 with lines
set ylabel "d2"
set output "d2.png"
plot "./mineresult.dat" u 1:8:9 with yerrorbars, "./mineresult.dat" u 1:8 with lines
set ylabel "d2f"
set output "d2f.png"
plot "./mineresult.dat" u 1:10:11 with yerrorbars, "./mineresult.dat" u 1:10 with lines
set ylabel "l"
set output "l.png"
plot "./mineresult.dat" u 1:12:13 with yerrorbars, "./mineresult.dat" u 1:12 with lines
set ylabel "lf"
set output "lf.png"
plot "./mineresult.dat" u 1:14:15 with yerrorbars, "./mineresult.dat" u 1:14 with lines
