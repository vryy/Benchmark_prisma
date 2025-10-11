#!/bin/sh

gnuplot -p << EOF
	set grid
	set title 'x-displacement of the flap tip'
	set xlabel 'time [s]'
	set ylabel 'x-displacement [m]'
	plot "top_deflection.txt" using 1:2 with lines notitle
EOF
