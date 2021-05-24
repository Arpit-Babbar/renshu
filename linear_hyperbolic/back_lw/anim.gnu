reset
set key #This is what enables setting legends
set grid
set yran[-1.1:1.1]
set xlabel 'x'
set ylabel 'u'
set title 'Plot of exact and approximate solution'
filename(n) = sprintf("solution_%d.txt",n)
N = 358
do for [i=1:N-1]{ 
   plot filename(i-1) u 1:2 t 'Approximation' w l, filename(i-1) u 1:3 t 'Exact' w l
   #plot filename(i-1) u 1:2 t 'Approximation' w lp lt 2 pt 6, filename(i-1) u 1:3 t 'Exact' w l
   pause 0.2
}

