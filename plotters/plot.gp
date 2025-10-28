set terminal png
set output "initial.png"

set view map

plot "../fort_solver/build/fort.1498" i 2 u 1:2:3 w image
