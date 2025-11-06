set terminal png
set output "fi3_fft_initial_no_short_rewritten.png"

set view map

plot "../build_rewritten/fort.1498" u 1:2:4 w image
