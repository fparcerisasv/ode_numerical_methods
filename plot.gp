set terminal pngcairo size 1200,600
set output "outputs/rkf45_solution.png"
set title "RKF45 Solution: x(t) and y(t)"
set xlabel "Time"
set ylabel "Solution"
plot "outputs/rkf45_output.dat" using 1:2 with lines title "x(t)", \
     "outputs/rkf45_output.dat" using 1:3 with lines title "y(t)"

set output "rkf45_error.png"


set terminal pngcairo size 1200,600
set output "outputs/taylor_solution.png"
set title "taylor Solution: x(t) and y(t)"
set xlabel "Time"
set ylabel "Solution"
plot "outputs/taylor_output.dat" using 1:2 with lines title "x(t)", \
     "outputs/taylor_output.dat" using 1:3 with lines title "y(t)"

set terminal pngcairo size 1200,600
set output "outputs/taylor_map_output.png"
set title "Map"
set xlabel "x"
set ylabel "y"
plot "outputs/taylor_map_output.dat" using 1:2 with lines title "phi(t)"


set terminal pngcairo size 1200,600
set output "outputs/taylor_map_large.png"
set title "Map"
set xlabel "x"
set ylabel "y"
plot "outputs/taylor_large_output.dat" using 1:2 with lines title "phi(t)"

