set terminal gif animate
set output "evolution.gif"

set view map

do for [idx = 0:10] {
    plot "../fort_solver/build/fort.498" index idx u 1:2:3 w image
}
