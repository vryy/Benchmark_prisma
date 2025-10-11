#!/bin/sh

gnuplot -p << EOF
	set grid
	set title 'y-displacement of the flap tip'
	set xlabel 'time [s]'
	set ylabel 'y-displacement [m]'
	plot "tip_deflection.txt" using 1:2 with lines notitle
EOF
