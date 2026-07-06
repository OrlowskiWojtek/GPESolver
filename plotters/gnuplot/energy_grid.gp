set term pdf
set output 'energy_grid.pdf'
set xlabel "n_{xy}" font ",20"
set ylabel  "<E/N> (nK)" font",20"
set size ratio 1
set xtics  font ",18"
set ytics  font ",18"

#set key off
#set label "y" font ",18" at 32,1.2 tc "blue"
#set label "x" font ",18" at 20,1. tc "red"

set colorscheme 
plot "data/energy_grid.dat" using 1:3 w l lw 2 lc rgb "blue" title "n_z=32", \
