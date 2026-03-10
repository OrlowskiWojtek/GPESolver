set margins
set term pdf
set output 'single_well_maxrho.pdf'
set size ratio 1
set key off
set xlabel "N/10^4" font ",16"
set ylabel "{/Symbol r}_{max} ({/Symbol m}m^{-3})" font ",16" textcolor "red"
set y2label "H_{z} ({/Symbol m}m)" font ",16" textcolor "blue"

set yrange[0:2500]
set y2range[2:12]
set xrange[0.1:5]

set xtics font",14" 1
set ytics font",14" 500
set y2tics font",14" 2

set arrow  from 0.55,0  to 0.55,2500 nohead dt 2
set arrow  from 1.05,0  to 1.05,2500 nohead dt 2
set arrow  from 2.35,0  to 2.35,2500 nohead dt 2
set arrow  from 3.35,0  to 3.35,2500 nohead dt 2

set label "I" at    0.25,1500 font",16"
set label "II" at   0.70,1000 font",16"
set label "III" at  1.5 ,1000 font",16"
set label "IV" at   2.8 ,1000 font",16"
set label "V" at    4.0 ,1000 font",16"

plot "data/single_well_maxrho.dat" u ($1 / 10.):($2*1000) w l lw 2 lc "red", \
     "data/single_well_height.dat" u ($1 / 10.):($2 / 1e3) w l lw 2 lc "blue" axis x1y2

reset
set term pdf
set output 'single_well_phase.pdf'
set size ratio 1
set key off
set xlabel "N/10^4" font ",16"
set ylabel "<E/N> (nK)" font ",16"
set xrange[0.1:5.0]
set yrange[-5:40]
set xtics font",14" 1
set ytics font",14" 

set arrow  from 0.55,-5  to 0.55,40 nohead dt 2
set arrow  from 1.05,-5  to 1.05,40 nohead dt 2
set arrow  from 2.35,-5  to 2.35,40 nohead dt 2
set arrow  from 3.35,-5  to 3.35,40 nohead dt 2

set label "I" at    0.25,30 font",16"
set label "II" at   0.70,30 font",16"
set label "III" at  1.5 ,30 font",16"
set label "IV" at   2.8 ,35 font",16"
set label "V" at    4.0 ,35 font",16"

set label "total" at 4.1,28 font ",16"
set label "V_{ext}" at 4.1,19.5 font ",16" tc "red"
set label "kinetic" at 3.9,9.5 font ",16" tc "orange"
set label "inter" at 4.1,2.2 font ",16" tc "blue"
plot "data/single_well_segments.dat"  u ($1 / 10.):(($3)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "orange",\
     "data/single_well_segments.dat"  u ($1 / 10.):(($5 + $6 + $7)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "blue",\
     "data/single_well_segments.dat"  u ($1 / 10.):(($4)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "red",\
     "data/single_well_segments.dat"  u ($1 / 10.):(($2)*27211.4*11.604*1e9/($1 * 1e3)) w l lw 2 lc "black"

