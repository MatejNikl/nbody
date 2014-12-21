set terminal postscript eps enhanced color
set output "graph_paropts_ofastpar.eps"

set title "12 threads, -Ofast, non-exact computation, parallelization options comparison"

set key bottom right
set key box width 1 height 1

set xlabel "Number of particles"
#set xrange[0:1400000]
set yrange[0:150]
set ytics 10
set ylabel "Performance [GFLOPS]"
set grid ytics

#set format x "%.2t*10e%T"

plot  "data/vec_12_0_4_ofastpar"         using 1:2 title "vec_{12,0,4}"         with linespoints lt 1 pt 11 lw 2 lc rgb "black",\
      "data/vec_12_0_4_dynamic_ofastpar" using 1:2 title "vec_{12,0,4,dynamic}" with linespoints lt 1 pt 11 lw 2 lc rgb "brown",\
      "data/vec_12_0_4_guided_ofastpar"  using 1:2 title "vec_{12,0,4,guided}"  with linespoints lt 1 pt 11 lw 2 lc rgb "pink",\
      "data/vec_12_0_4_shared_ofastpar"  using 1:2 title "vec_{12,0,4,shared}"  with linespoints lt 1 pt 11 lw 2 lc rgb "blue",\
      "data/unroll_2_0_ofastpar"         using 1:2 title "unroll_{2,0}"         with linespoints lt 1 pt 11 lw 2 lc rgb "orange",\
      "data/unroll_2_0_jpar_ofastpar"    using 1:2 title "unroll_{2,0,jpar}"    with linespoints lt 1 pt 11 lw 2 lc rgb "red"
