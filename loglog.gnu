set key #This is what enables setting legends
set logscale y
set logscale x
set grid
set xlabel 'h'
set ylabel 'error'
set title 'Plot of error V/S h'
filename(n) = sprintf("error_vs_h.txt",n)
N = 358
do for [i=1:N-1] {
   plot filename(i) u 1:2 t 'rk3' w l, filename(i) u 1:3 t 'rk3 new' w l
   pause 200.0
}

