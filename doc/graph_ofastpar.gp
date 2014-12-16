set terminal postscript eps enhanced color
set output "graph_ofastpar.eps"

set title "12 threads, -Ofast, non-exact computation"

set key bottom right
set key box width 1 height 1

set xlabel "Number of particles"
#set xrange[0:1400000]
set yrange[0:150]
set ytics 10
set ylabel "Performance [GFLOPS]"
set grid ytics

#set format x "%.2t*10e%T"

plot  "data/vec_12_0_4_ofastpar" using 1:2 title "vec_{12,0,4}" with linespoints lt 1 pt 11 lw 2 lc rgb "black",\
      "data/vec_16_0_4_ofastpar" using 1:2 title "vec_{16,0,4}" with linespoints lt 1 pt 11 lw 2 lc rgb "brown",\
      "data/vec_8_0_4_ofastpar"  using 1:2 title "vec_{8,0,4}"  with linespoints lt 1 pt 11 lw 2 lc rgb "pink",\
      "data/vec_4_0_4_ofastpar"  using 1:2 title "vec_{4,0,4}"  with linespoints lt 1 pt 11 lw 2 lc rgb "blue",\
      "data/unroll_2_0_ofastpar" using 1:2 title "unroll_{2,0}" with linespoints lt 1 pt 11 lw 2 lc rgb "orange",\
      "data/naive_ofastpar"      using 1:2 title "naive"        with linespoints lt 1 pt 11 lw 2 lc rgb "black",\
      "data/unroll_4_0_ofastpar" using 1:2 title "unroll_{4,0}" with linespoints lt 1 pt 11 lw 2 lc rgb "red",\
      "data/unroll_8_0_ofastpar" using 1:2 title "unroll_{8,0}" with linespoints lt 1 pt 11 lw 2 lc rgb "green",\
      "data/unroll_2_2_ofastpar" using 1:2 title "unroll_{2,2}" with linespoints lt 1 pt 11 lw 2 lc rgb "cyan",\
      "data/unroll_0_4_ofastpar" using 1:2 title "unroll_{0,4}" with linespoints lt 1 pt 11 lw 2 lc rgb "olive",\
      "data/unroll_0_2_ofastpar" using 1:2 title "unroll_{0,2}" with linespoints lt 1 pt 11 lw 2 lc rgb "purple",\
      "data/unroll_0_8_ofastpar" using 1:2 title "unroll_{0,8}" with linespoints lt 1 pt 11 lw 2 lc rgb "yellow",\
