set terminal postscript eps enhanced color
set output "graph_threads_ofastpar.eps"

set title "-Ofast, non-exact computation, # of threads comparison"

set key bottom right
set key box width 1 height 1

set xlabel "Number of threads"
set xrange[0:13]
set yrange[0:150]
set ytics 10
set ylabel "Performance [GFLOPS]"
set grid ytics

#set format x "%.2t*10e%T"

plot  "data/vec_12_0_4_100000_ofastpar"  using 1:2 title "vec_{12,0,4} 100000 particles"  with linespoints lt 1 pt 11 lw 2 lc rgb "black",\
      "data/vec_12_0_4_400000_ofastpar"  using 1:2 title "vec_{12,0,4} 400000 particles"  with linespoints lt 1 pt 11 lw 2 lc rgb "blue",\
      "data/vec_12_0_4_700000_ofastpar"  using 1:2 title "vec_{12,0,4} 700000 particles"  with linespoints lt 1 pt 11 lw 2 lc rgb "red",\
      "data/vec_12_0_4_1000000_ofastpar" using 1:2 title "vec_{12,0,4} 1000000 particles" with linespoints lt 1 pt 11 lw 2 lc rgb "green",\
      "data/vec_12_0_4_1300000_ofastpar" using 1:2 title "vec_{12,0,4} 1300000 particles" with linespoints lt 1 pt 11 lw 2 lc rgb "pink",\
