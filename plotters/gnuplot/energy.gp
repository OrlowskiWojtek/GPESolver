#set output "fft_energy.png"

set terminal qt persist
file = "../build/data_3_bce/energy.dat"

plot file using 1:7 with lines title "Energy vs Iteration", \
     file using 1:2 with lines title "Kinetic Energy vs Iteration", \
     file using 1:3 with lines title "Potential Energy vs Iteration", \
     file using 1:4 with lines title "Interaction Energy vs Iteration", \
     file using 1:5 with lines title "Dipole-dipole Energy vs Iteration", \
     file using 1:6 with lines title "BMF Energy vs Iteration"

replot
