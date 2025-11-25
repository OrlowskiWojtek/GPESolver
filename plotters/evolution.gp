set terminal gif animate delay 10
set output "fft_evolution.gif"

set view map

set cbrange [0:6e-15]

do for [idx = 0:100] { 
    iter = idx * 100
    plot "../build/cut".iter.".dat" u 1:2:3 w image
}
