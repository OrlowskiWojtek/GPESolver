set margins
set term pdf
set output 'eps_145_maxrho.pdf'
set size ratio 1
set key off
set xlabel "N/10^4" font ",16"
set ylabel "{/Symbol r}_{max} ({/Symbol m}m^{-3})" font ",16" textcolor "red"
set y2label "H_{z} ({/Symbol m}m)" font ",16" textcolor "blue"

#set yrange[0:1700]
#set y2range[4:12]
set xrange[0.1:5.4]

set xtics font",14" 1
set ytics font",14" 500
set y2tics font",14" 2

set arrow  from 2.55,0  to 2.55,2000 nohead dt 2
set arrow  from 2.95,0  to 2.95,2000 nohead dt 2

set label "I" at 1.5,600 font",16"
set label "II" at 2.6,600 font",16"
set label "III" at 4.0,600 font",16"

plot "data/eps_145_maxrho.dat" u ($1 / 10.):($2*1000) w l lw 2 lc "red", \
     "data/eps_145_height.dat" u ($1 / 10.):($2 / 1e3) w l lw 2 lc "blue" axis x1y2

reset
set term pdf
set output 'eps_145_phase.pdf'
set size ratio 1
set key off
set xlabel "N/10^4" font ",16"
set ylabel "<E/N> (nK)" font ",16"
set xrange[0.1:5.4]
set yrange[0:40]
set xtics font",14" 1
set ytics font",14" 

set arrow  from 2.55,0  to 2.55,40 nohead dt 2
set arrow  from 2.95,0  to 2.95,40 nohead dt 2

set label "I" at 1.5,30 font",16"
set label "II" at 2.6,33 font",16"
set label "III" at 4.0,35 font",16"

set label "total" at 4.1,28 font ",16"
set label "V_{ext}" at 4.1,19.5 font ",16" tc "red"
set label "kinetic" at 3.9,2.2 font ",16" tc "orange"
set label "inter" at 4.1,9.5 font ",16" tc "blue"
plot "data/eps_145_segments.dat"  u ($1 / 10.):(($3)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "orange",\
     "data/eps_145_segments.dat"  u ($1 / 10.):(($5 + $6 + $7)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "blue",\
     "data/eps_145_segments.dat"  u ($1 / 10.):(($4)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "red",\
     "data/eps_145_segments.dat"  u ($1 / 10.):(($2)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "black"

