gnuplot <<_end
set terminal eps enhanced
set style line 1 lc 1 lw 2
set style line 2 lt 2 lc 2 lw 2
set style line 3 lt 2 lc 3 lw 2
set style line 4 lt 2 lc 4 lw 2
set style line 5 lt 2 lc 5 lw 2
set xtics ("1,1" 6, "2,2" 19, "2,3" 32)
set xzeroaxis
set termoption dashed

set title 'FN-DMC Hessian with noise'
set ylabel 'Energy 2nd Derivative'
set xlabel 'Parameters'
set output 'fnHessianNoise.eps'


plot 'summary3_noise' u 3:4 w l ls 1 t 'FN-DMC exact', 'summary3_noise' u 3:5 w l ls 2 t 'FU2000', 'summary3_noise' u 3:6 w l ls 3 t 'VMC', 'summary3_noise' u 3:7 w l ls 4 t 'VMC expr. in DMC', 'summary3_noise' u ($3+13):8 w l ls 1 notitle, 'summary3_noise' u ($3+13):9 w l ls 2 notitle, 'summary3_noise' u ($3+13):10 w l ls 3 notitle, 'summary3_noise' u ($3+13):11 w l ls 4 notitle, 'summary3_noise' u ($3+26):12 w l ls 1 notitle, 'summary3_noise' u ($3+26):13 w l ls 2 notitle, 'summary3_noise' u ($3+26):14 w l ls 3 notitle, 'summary3_noise' u ($3+26):15 w l ls 4 notitle

set output 'fnHessianNoNoise.eps'
set title 'FN-DMC Hessian without noise'
plot 'summary3_noNoise' u 3:4 w l ls 1 t 'FN-DMC exact', 'summary3_noNoise' u 3:5 w l ls 2 t 'FU2000', 'summary3_noNoise' u 3:6 w l ls 3 t 'VMC', 'summary3_noNoise' u 3:7 w l ls 4 t 'VMC expr. in DMC', 'summary3_noNoise' u ($3+13):8 w l ls 1 notitle, 'summary3_noNoise' u ($3+13):9 w l ls 2 notitle, 'summary3_noNoise' u ($3+13):10 w l ls 3 notitle, 'summary3_noNoise' u ($3+13):11 w l ls 4 notitle, 'summary3_noNoise' u ($3+26):12 w l ls 1 notitle, 'summary3_noNoise' u ($3+26):13 w l ls 2 notitle, 'summary3_noNoise' u ($3+26):14 w l ls 3 notitle, 'summary3_noNoise' u ($3+26):15 w l ls 4 notitle

_end
