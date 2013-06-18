set title "Графическое представление хода алгоритма"
set xrange [0.2:0.8]
set yrange [0.6:2]
plot sqrt(1-x**2) t 'x^2 + y^2 -1 = 0',\
	     x+0.5 t 'x - y + 0.5 = 0'
replot 'xs.dat' using 1:2 with lines t 'approximations'
