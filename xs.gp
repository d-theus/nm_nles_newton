set title "Графическое представление хода алгоритма"
set xrange [0.2:0.8]
set yrange [0.6:2]
unset key
plot sqrt(1-x**2),\
	     x+0.5
replot 'xs.dat' using 1:2 with lines
