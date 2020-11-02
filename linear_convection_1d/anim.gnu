reset
unset key
set grid
set yran[-1.1:1.1]
filename(n) = sprintf("solution_%d_lw.txt",n)
N = 358
do for [i=1:N-1] {
   plot filename(i) w l
   pause 0.2
}
