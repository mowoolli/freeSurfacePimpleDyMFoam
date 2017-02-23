N=3

set xlabel 't [s]'
set ylabel 'heave [m], pitch [deg]'
plot "RBmotion/position.dat" u 1:($4+1.2738038) every N::(N-1) w l lt 1 title "heave"
replot "RBmotion/position.dat" u 1:($6*180/3.141592) every N::(N-1) w l lt 2 title "pitch"
replot "<sed 's/(/ /g' RBmotion/position.dat" u 1:($5*180/3.141592) every N::(N-1) w l lt 3 title "roll"
