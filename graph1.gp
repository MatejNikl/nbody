set terminal postscript eps enhanced color
set output "graph1.eps"

set title ""

set key top right
set key box width 1 height 1

set xlabel "Number of particles"
set yrange[0:11]
set ytics 1
set ylabel "Performance [GFLOPS]"
set grid ytics

#set format x "%.2t*10e%T"

plot  "data/naive_o3_gflops.txt" using 1:2 title "naive -O3" with linespoints lt 4 pt 11 lw 2 lc 2,\
      "data/naive_ofast_gflops.txt" using 1:2 title "naive -Ofast" with linespoints lt 4 pt 11 lw 2 lc 3,\
      "data/vec_o3_gflops.txt" using 1:2 title "vect + unroll -O3" with linespoints lt 4 pt 11 lw 2 lc 4,\
      "data/vec_ofast_gflops.txt" using 1:2 title "vect + unroll -Ofast" with linespoints lt 4 pt 11 lw 2 lc 5

