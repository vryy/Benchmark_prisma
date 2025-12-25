set terminal pdf enhanced color

M = 1.4
N = 1.4
R = 2.0
pc0 = 3.967271e+01
p0 = 2.630831e+01
q0 = 2.534678e+01
pTr = 3.793501e+00
qTr = 1.652658e+01

pk = 5.194767e+00
qk = 1.559380e+01
pck = 3.902063e+01
pTrk = 4.179837e+00
qTrk = 1.649238e+01

set xrange [0:50]
set yrange [0:40]

set output "loading_path.pdf"
set xlabel "p' (Pa)"
set ylabel "q (Pa)"
set parametric
set trange [0.00001:1]

alpha0 = 0.000000e+00; p0 = 2.630831e+01; q0 = 2.534678e+01; pc0 = 3.967271e+01
alpha1 = 6.300000e-01; p1 = 7.766578e+00; q1 = 1.363582e+01; pc1 = 3.967271e+01
alpha2 = 8.631000e-01; p2 = 5.114535e+00; q2 = 1.549614e+01; pc2 = 3.944505e+01
alpha3 = 8.716595e-01; p3 = 5.142005e+00; q3 = 1.552973e+01; pc3 = 3.929748e+01
alpha4 = 8.802189e-01; p4 = 5.169874e+00; q4 = 1.556363e+01; pc4 = 3.915014e+01
alpha5 = 8.856114e-01; p5 = 5.187669e+00; q5 = 1.558519e+01; pc5 = 3.905740e+01
alpha6 = 8.877517e-01; p6 = 5.194767e+00; q6 = 1.559380e+01; pc6 = 3.902063e+01

plot pc0*t, M*pc0*t*((-log(t)/log(R))**(1.0/N)) notitle, \
     pc1*t, M*pc1*t*((-log(t)/log(R))**(1.0/N)) notitle, \
     pc2*t, M*pc2*t*((-log(t)/log(R))**(1.0/N)) notitle, \
     pc3*t, M*pc3*t*((-log(t)/log(R))**(1.0/N)) notitle, \
     pc4*t, M*pc4*t*((-log(t)/log(R))**(1.0/N)) notitle, \
     pc5*t, M*pc5*t*((-log(t)/log(R))**(1.0/N)) notitle, \
     pc6*t, M*pc6*t*((-log(t)/log(R))**(1.0/N)) notitle, \
     t*3e1, M*t*3e1 title "CSL", \
     p0+(p1-p0)*t, q0+(q1-q0)*t ls 1 title "loading path", \
     p1+(p2-p1)*t, q1+(q2-q1)*t ls 1 notitle, \
     p2+(p3-p2)*t, q2+(q3-q2)*t ls 1 notitle, \
     p3+(p4-p3)*t, q3+(q4-q3)*t ls 1 notitle, \
     p4+(p5-p4)*t, q4+(q5-q4)*t ls 1 notitle, \
     p5+(p6-p5)*t, q5+(q6-q5)*t ls 1 notitle
