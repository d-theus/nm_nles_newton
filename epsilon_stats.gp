set title "Зависимость количества итераций от требуемой точности"
set logscale x
unset key
set xlabel "Precision"
set ylabel "Iterations"
set xrange [10e-10:1]
set grid
plot 'epsstats.dat' using 1:2 with lines lw 2
