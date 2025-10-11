set output "report.pdf"
set terminal pdf
set view equal xy
#set size ratio -1
set key right bottom
#set nokey
#set yrange [0:0.8]
set xlabel "displacement, ux (mm)"
set ylabel "traction (N)"
set style line 1 lc rgb 'red' lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb 'blue' lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 3 lc rgb 'green' lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 4 lc rgb 'cyan' lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 5 lc rgb 'magenta' lt 1 lw 2 pt 5 pi -0.5 ps 0.5

plot "monitoring.no_plastic_heat.log" using ($1):(-$2) with lines title "J2 (backward Euler, no plastic heat) solution", \
     "monitoring.log" using ($1):(-$2) with linespoints title "J2 (backward Euler) solution", \

set key left top
set ylabel "temperature (oC)"
plot "monitoring.log" using ($1):($3) with linespoints title "J2 (backward Euler) solution", \

