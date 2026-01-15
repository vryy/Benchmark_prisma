set term pdf size 9,6
set output "Flux_text.pdf"
set grid
set  key font "CMU-Serif, 18.0"
set  xlabel font " CMU-Serif, 22.0"
set  ylabel font " CMU-Serif, 22.0"
set key spacing 1
set xtics font " CMU-Serif, 21.0"
set ytics font " CMU-Serif, 21.0"
set xlabel "Time [min]" offset -5
set ylabel "Darcy Velocity [cm/min]" offset -1
set key at 120,0.0275
set xrange [0:120]
set xtics 20
set ytics 0.005
set pointsize 0.8
unset grid
plot "../../../plotting/experiments/flux/FluxExp_time.txt" using 1:2 with points pt 5 ps 0.5 lt rgb "gray" title "Exp. (Liakopoulos 1965)", \
     "../../../plotting/Nagel_2phase/flux/FluxNagel_time.txt" using 1:2  with linespoints pt 5 ps 0.4 lw 1 lt rgb "red" title "One phase flow model (Nagel 2010)", \
     "../../../plotting/vp-PFEM-Model/flux/FluxSim_time.grf" using 1:2 with lines lw 1 lt rgb "green" title "v-p PFEM model", \
     "bottom_water_flow.log" using ($1/60):(-$2*60.0/0.01) with lines lw 1 lt rgb "blue" title "partially saturated soils model (lateral permeable)", \
