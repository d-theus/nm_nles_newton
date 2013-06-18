# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 500, 350 
# set output 'surface1.1.png'
set title "Количество итераций в зависимости от выбора начального приближения"
set isosamples 50
set grid
set xrange [-0.1:0.8]
set xlabel 'X'
set yrange [0.2:1.3]
set ylabel 'Y'
set zrange [-0.01:40]
set ztics 8
f(x,y) = x**2 + y**2 - 1
g(x,y) = x - y + 0.5
set contour base
set cntrparam levels disc 0
unset surface
set style line 1 lc rgb 'blue' lw 2
set style line 2 lc rgb 'orange' lw 2
set style line 3 lc rgb 'blue' lw 2
set style line 4 lc rgb 'blue' lw 2
set style increment userstyle

splot 		f(x,y) with lines t 'x^2 + y^2 - 1 = 0',\
	  	g(x,y) with lines t 'x - y + 0.5 = 0',\
		'p0ok.dat' with vectors linecolor rgb 'red' lw 1.5 t 'Iterations'
set key default
