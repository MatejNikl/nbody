set terminal postscript eps enhanced color
set output "graph_threads_ofastpar.eps"

set title "-Ofast, non-exact computation, # of threads comparison"

set key bottom right
set key box width 1 height 1

set xlabel "Number of particles"
set xrange[0:4100000]
set yrange[0:150]
set ytics 10
set ylabel "Performance [GFLOPS]"
set grid ytics

#set format x "%.2t*10e%T"

plot  "data/vec_12_0_4_ofastpar" using 1:2 title "vec_{12,0,4} 12 threads" with linespoints lt 1 pt 11 lw 2 lc rgb "black",\
      "data/vec_12_0_4_10_ofastpar"  using 1:2 title "vec_{12,0,4} 10 threads"  with linespoints lt 1 pt 11 lw 2 lc rgb "blue",\
      "data/vec_12_0_4_8_ofastpar"  using 1:2 title "vec_{12,0,4} 8 threads"  with linespoints lt 1 pt 11 lw 2 lc rgb "red",\
      "data/vec_12_0_4_6_ofastpar"  using 1:2 title "vec_{12,0,4} 6 threads"  with linespoints lt 1 pt 11 lw 2 lc rgb "green",\
      "data/vec_12_0_4_4_ofastpar"  using 1:2 title "vec_{12,0,4} 4 threads"  with linespoints lt 1 pt 11 lw 2 lc rgb "pink",\
      "data/vec_12_0_4_2_ofastpar" using 1:2 title "vec_{12,0,4} 2 threads" with linespoints lt 1 pt 11 lw 2 lc rgb "brown",\
      "data/vec_12_0_4_ofast" using 1:2 title "vec_{12,0,4} 1 thread" with linespoints lt 1 pt 11 lw 2 lc rgb "orange",\
