# set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 500, 350 
# set output 'surface1.1.png'
set title "Количество итераций в зависимости от выбора начального приближения"
set isosamples 50
set grid
set xrange [-0.1:0.8]
set yrange [0.2:1.3]
set zrange [-0.01:40]
f(x,y) = x**2 + y**2 - 1
g(x,y) = x - y + 0.5
unset key
set contour
set cntrparam levels disc 0
unset surface
splot f(x,y), g(x,y), 'p0.dat' with vectors linecolor rgb 'red' lw 2 filled
