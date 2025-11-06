set terminal gif animate
set output "fft_evolution.gif"

set view map

do for [idx = 0:30] {
    plot "../build_rewritten/fort.498" index idx u 1:2:3 w image
}
