set margins
set term pdf
set output 'eps_140_maxrho.pdf'
set size ratio 1
set key off
set xlabel "N/10^4" font ",16"
set ylabel "{/Symbol r}_{max} ({/Symbol m}m^{-3})" font ",16" textcolor "red"
set y2label "H_{z} ({/Symbol m}m)" font ",16" textcolor "blue"

#set yrange[0:2500]
set y2range[4:10]
#set xrange[0.1:5]

set xtics font",14" 1
set ytics font",14" 400
set y2tics font",14" 2

set arrow  from 2.65,0  to 2.65,1500 nohead dt 2

set label "I" at    1.3,1000 font",16"
set label "II" at   3.9,500 font",16"

plot "data/eps_140_maxrho.dat" u ($1 / 10.):($2*1000) w l lw 2 lc "red", \
     "data/eps_140_height.dat" u ($1 / 10.):($2 / 1e3) w l lw 2 lc "blue" axis x1y2

reset
set term pdf
set output 'eps_140_phase.pdf'
set size ratio 1
set key off
set xlabel "N/10^4" font ",16"
set ylabel "<E/N> (nK)" font ",16"
set xrange[0.1:5.0]
set yrange[-5:40]
set xtics font",14" 1
set ytics font",14" 

set arrow  from 2.65,-5  to 2.65,40 nohead dt 2

set label "I" at    1.3,30 font",16"
set label "II" at   3.5,25 font",16"

set label "total" at 4.1,34 font ",16"
set label "V_{ext}" at 4.1,19.5 font ",16" tc "red"
set label "kinetic" at 3.9,1.9 font ",16" tc "orange"
set label "inter" at 4.1,9.5 font ",16" tc "blue"

plot "data/eps_140_segments.dat"  u ($1 / 10.):(($3)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "orange",\
     "data/eps_140_segments.dat"  u ($1 / 10.):(($5 + $6 + $7)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "blue",\
     "data/eps_140_segments.dat"  u ($1 / 10.):(($4)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "red",\
     "data/eps_140_segments.dat"  u ($1 / 10.):(($2)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "black"


