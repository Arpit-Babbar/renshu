reset
set key #This is what enables setting legends
set grid
set yran[-1.1:1.1]
set xlabel 'x'
set ylabel 'u'
set title 'Plot of exact and approximate solution'
filename(n) = sprintf("solution_%d_ftcs.txt",n)
N = 358
do for [i=1:N-1] {
   plot filename(i) u 1:2 t 'Approximation' w l, filename(i) u 1:3 t 'Exact' w l
   pause 2.0
}

