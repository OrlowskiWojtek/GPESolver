set terminal gif animate delay 10
set output "fft_evolution.gif"

set view map

do for [idx = 0:300] { 
    iter = idx * 100
    plot "../build/cut".iter.".dat" u 1:2:3 w image
}
