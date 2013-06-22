set title "Зависимость количества итераций от требуемой точности\nТолько метод наискорейшего спуска"
set logscale xy
unset key
set xlabel "Precision"
set ylabel "Iterations"
set xrange [10e-11:1]
set grid
plot 'epsstats_gd_only.dat' using 1:2 with lines lw 2
