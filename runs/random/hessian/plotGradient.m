gnuplot <<_end
set terminal eps enhanced
set style line 1 lc 1 lw 2 
set style line 2 lt 2 lc 2 lw 2 
set style line 3 lt 2 lc 3 lw 2  
set style line 4 lt 2 lc 4 lw 2 
set style line 5 lt 2 lc 5 lw 2 
set xtics ("1" 6, "2" 19, "3" 32)
set yrange[-0.003:0.004]
set xzeroaxis
set termoption dashed
set title 'FN-DMC gradient with noise'
set ylabel 'Energy derivative'
set xlabel 'Parameter'
set output 'fnGradientNoise.eps'

plot 'summary2_noise' u 3:4 w l ls 1 t 'FN-DMC exact', 'summary2_noise' u 3:5 w l ls 2 t 'FU2000 approx', 'summary2_noise' u 3:6 w l ls 3 t 'New', 'summary2_noise' u 3:7 w l ls 4 t 'VMC', 'summary2_noise' u 3:8 w l ls 5 t 'VMC expr. in DMC', 'summary2_noise' u ($3+13):9 w l ls 1 notitle , 'summary2_noise' u ($3+13):10 w l ls 2 notitle, 'summary2_noise' u ($3+13):11 w l ls 3 notitle, 'summary2_noise' u ($3+13):12 w l ls 4 notitle, 'summary2_noise' u ($3+13):13 w l ls 5 notitle, 'summary2_noise' u ($3+26):14 w l ls 1 notitle, 'summary2_noise' u ($3+26):15 w l ls 2 notitle, 'summary2_noise' u ($3+26):16 w l ls 3 notitle, 'summary2_noise' u ($3+26):17 w l ls 4 notitle, 'summary2_noise' u ($3+26):18 w l ls 5 notitle 

set output 'fnGradientNoNoise.eps'
set title 'FN-DMC gradient without noise'
plot 'summary2_noNoise' u 3:4 w l ls 1 t 'FN-DMC exact', 'summary2_noNoise' u 3:5 w l ls 2 t 'FU2000 approx', 'summary2_noNoise' u 3:6 w l ls 3 t 'New', 'summary2_noNoise' u 3:7 w l ls 4 t 'VMC', 'summary2_noNoise' u 3:8 w l ls 5 t 'VMC expr. in DMC', 'summary2_noNoise' u ($3+13):9 w l ls 1 notitle , 'summary2_noNoise' u ($3+13):10 w l ls 2 notitle, 'summary2_noNoise' u ($3+13):11 w l ls 3 notitle, 'summary2_noNoise' u ($3+13):12 w l ls 4  notitle, 'summary2_noNoise' u ($3+13):13 w l ls 5 notitle, 'summary2_noNoise' u ($3+26):14 w l ls 1 notitle, 'summary2_noNoise' u ($3+26):15 w l ls 2 notitle, 'summary2_noNoise' u ($3+26):16 w l ls 3 notitle, 'summary2_noNoise' u ($3+26):17 w l ls 4 notitle, 'summary2_noNoise' u ($3+26):18 w l ls 5 notitle 

_end
