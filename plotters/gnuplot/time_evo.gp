set term pdf
set output 'n30000_time.pdf'
set xlabel "t (ms)" font ",20"
set ylabel  "r_{max} ({/Symbol m}m)" font",20"
set size ratio 1
set xrange[:75]
set xtics  font ",18"
set key off
set label "y" font ",18" at 32,-0.2 tc "blue"
set label "x" font ",18" at 20,0.5 tc "red"
set yrange[-2.5:2.5]
set ytics  font ",18"

plot "data/n30000_time_new.dat" u ($1):(-$2/1000):($4*$4/2e7) w p pt 5 ps var lc "red" t "x", \
     "data/n30000_time_new.dat" u ($1):($3/1000):($4*$4/2e7) w p pt 5 ps var lc "blue" t "y" 
