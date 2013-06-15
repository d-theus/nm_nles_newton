set title "Зависимость количества итераций от начального приближения"
set logscale
unset key
set xlabel "epsilon"
set ylabel "i"
set xrange [10e-8:1]
plot 'epsstats.dat' using 1:2 with lines
