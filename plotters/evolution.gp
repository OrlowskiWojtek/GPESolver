set terminal gif animate
set output "fft_evolution.gif"

set view map

do for [idx = 0:64] {
    iter = idx * 1000
    plot "../build/cut".iter.".dat" u 1:2:3 w image
}
