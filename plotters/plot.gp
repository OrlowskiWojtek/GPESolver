#set terminal png
#set output "fi3_fft_initial_no_short_rewritten.png"

set view map

#plot "../fort_solver/build/fort.1498" index 5 u 1:2:4 w image
#plot "../build/xcut.dat" u 1:4
plot "../build/cut.dat" u 1:2:3 w image
pause -1
