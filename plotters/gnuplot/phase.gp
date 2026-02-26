set term pdf
set output 'eps_15_maxrho.pdf'
set size ratio 1
set key off
set xlabel "N/10^4" font ",16"
set ylabel "{/Symbol r}_{max} ({/Symbol m}m^{-3})" font ",16"

set yrange[0:2500]
set xrange[0.1:5]

set xtics font",14" 1
set ytics font",14"

set arrow  from 0.55,0  to 0.55,2500 nohead dt 2
set arrow  from 0.75,0  to 0.75,2500 nohead dt 2
set arrow  from 2.95,0  to 2.95,2500 nohead dt 2
set arrow  from 3.65,0  to 3.65,2500 nohead dt 2

set label "I" at    0.25,2000 font",16"
set label "II" at   0.55,2000 font",16"
set label "III" at  1.5 ,2000 font",16"
set label "IV" at   3.1 ,2000 font",16"
set label "V" at    4.0 ,2000 font",16"

plot "data/eps_15_maxrho.dat" u ($1 / 10.):($2*1000) w l lw 2 lc "red"

reset
set term pdf
set output 'eps_15_phase.pdf'
set size ratio 1
set key off
set xlabel "N/10^4" font ",16"
set ylabel "<E/N> (nK)" font ",16"
set xrange[0.1:5.0]
set yrange[-5:40]
set xtics font",14" 1
set ytics font",14" 

set arrow  from 0.55,-5  to 0.55,40 nohead dt 2
set arrow  from 0.75,-5  to 0.75,40 nohead dt 2
set arrow  from 2.95,-5  to 2.95,40 nohead dt 2
set arrow  from 3.65,-5  to 3.65,40 nohead dt 2

set label "I" at 0.25,30 font",16"
set label "II" at 0.55,30 font",16"
set label "III" at 1.5,30 font",16"
set label "IV" at 3.1,35 font",16"
set label "V" at 4.0,35 font",16"

set label "total" at 4.1,28 font ",16"
set label "V_{ext}" at 4.1,19.5 font ",16" tc "red"
set label "kinetic" at 3.9,9.5 font ",16" tc "orange"
set label "inter" at 4.1,2.2 font ",16" tc "blue"
plot "data/eps_15_segments.dat"  u ($1 / 10.):(($3)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "orange",\
     "data/eps_15_segments.dat"  u ($1 / 10.):(($5 + $6 + $7)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "blue",\
     "data/eps_15_segments.dat"  u ($1 / 10.):(($4)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "red",\
     "data/eps_15_segments.dat"  u ($1 / 10.):(($2)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "black"

