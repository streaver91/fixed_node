#set style data linespoints
#plot [.99:1.01][.9:1.1] '1', '2', '3', '4', '5', '1' u 1:(1) w lines lt -1 notitle, '1' u (1):2 w lines lt -1 notitle
#set arrow from graph 0, first 1 to graph 1, first 1 nohead lt 0  # draws horizontal line at y=1
#set arrow from 1, graph 0 to 1, graph 1 nohead lt 0  # draws vertical line at x=0
plot [0:] 'summary28' u 2:4 title 'exact 2nd term', \
          'summary28' u 2:5 title '3rd term', \
          'summary28' u 2:6 title 'W_i/W'
          'summary28' u 2:6 title 'W_i/W'
