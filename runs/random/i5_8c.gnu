# This plot is for 10 different random Hamiltonians and corresponding approximate trial wavefunctions
# This plot demonstrates a few things:
# 1) As expected the FN-DMC gradient has much smaller magnitude than the VMC gradient
# 2) The VMC gradient and the VMC expression for the gradient evaluated for the DMC distribution are very close
# 3) The approximate DMC gradients are accurate, i.e, the zero of the true DMC gradient is much closer to the zero of the approx. DMC gradients than to the zero of the VMC gradient 
# 4) The new approximation for the FN-DMC gradient is more accurate than the FU2000 approximation.
#    In particular, FU2000 sometimes has a 2nd spurious zero.
#set style data linespoints
set title "Parameter derivative of energy for 10 Hamiltonians and guiding/trial wavefunctions" font "Times-Roman,18"
set xlabel "Shifted data sets (wavefn. parameter increases within each set)"
set ylabel "Energy derivative"
set xtics ("1" 5, "2" 15, "3" 25, "4" 35, "5" 45, "6" 55, "7" 65, "8" 75, "9" 85, "10" 95)
set style data lines
#plot [.99:1.01][.9:1.1] '1', '2', '3', '4', '5', '1' u 1:(1) w lines lt -1 notitle, '1' u (1):2 w lines lt -1 notitle
set arrow from graph 0, first 0 to graph 1, first 0 nohead lt 0  # draws horizontal line at y=0
#set arrow from 1, graph 0 to 1, graph 1 nohead lt 0  # draws vertical line at x=0
#plot [0:105] \
#  'gradient8.dat' u (10*($1-1)+$2):8 title 'VMC gradient', \
#  'gradient8.dat' u (10*($1-1)+$2):9 title 'VMC gradient expression in DMC', \
#  'gradient8.dat' u (10*($1-1)+$2):10 title 'DMC gradient', \
#  'gradient8.dat' u (10*($1-1)+$2):11 title 'approx DMC gradient new', \
#  'gradient8.dat' u (10*($1-1)+$2):12 title 'approx DMC gradient FU2000'

plot [0:105][-.08:] \
  'gradient8.dat' u (10*($1-1)+$2):8 title 'VMC gradient', \
  'gradient8.dat' u (10*($1-1)+$2):10 title 'DMC gradient', \
  'gradient8.dat' u (10*($1-1)+$2):11 title 'approx DMC gradient'

set size 1.,1.; set term post eps enhanced color solid "Times-Roman" 24 ; set output 'energy_gradient8c.eps' ; replot
