set terminal postscript eps enhanced color
set output "graph_o3par.eps"

set title "12 threads, -O3, non-exact computation"

set key bottom right
set key box width 1 height 1

set xlabel "Number of particles"
#set xrange[0:1400000]
set yrange[0:50]
set ytics 10
set ylabel "Performance [GFLOPS]"
set grid ytics

#set format x "%.2t*10e%T"

plot  "data/vec_12_0_4_o3par" using 1:2 title "vec_{12,0,4}" with linespoints lt 1 pt 11 lw 2 lc 0,\
      "data/vec_16_0_4_o3par" using 1:2 title "vec_{16,0,4}" with linespoints lt 1 pt 11 lw 2 lc 1,\
      "data/vec_8_0_4_o3par"  using 1:2 title "vec_{8,0,4}"  with linespoints lt 1 pt 11 lw 2 lc 2,\
      "data/vec_4_0_4_o3par"  using 1:2 title "vec_{4,0,4}"  with linespoints lt 1 pt 11 lw 2 lc 3,\
      "data/unroll_2_0_o3par" using 1:2 title "unroll_{2,0}" with linespoints lt 1 pt 11 lw 2 lc 4,\
      "data/naive_o3par"      using 1:2 title "naive"        with linespoints lt 1 pt 11 lw 2 lc 5,\
      "data/unroll_4_0_o3par" using 1:2 title "unroll_{4,0}" with linespoints lt 1 pt 11 lw 2 lc 6,\
      "data/unroll_8_0_o3par" using 1:2 title "unroll_{8,0}" with linespoints lt 1 pt 11 lw 2 lc 7,\
      "data/unroll_2_2_o3par" using 1:2 title "unroll_{2,2}" with linespoints lt 1 pt 11 lw 2 lc 8,\
      "data/unroll_0_4_o3par" using 1:2 title "unroll_{0,4}" with linespoints lt 1 pt 11 lw 2 lc 9,\
      "data/unroll_0_2_o3par" using 1:2 title "unroll_{0,2}" with linespoints lt 1 pt 11 lw 2 lc 10,\
      "data/unroll_0_8_o3par" using 1:2 title "unroll_{0,8}" with linespoints lt 1 pt 11 lw 2 lc 11,\
