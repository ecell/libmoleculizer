
HISTO_DIR=$1

gnuplot > $HISTO_DIR/dr-curve.png <<EOF
set term png
set yrange [0:*]
set logscale x
set xlabel "Dose number of molecules"
set ylabel "Mean response number of molecules"
plot "$HISTO_DIR/all-stats" using 1:2 notitle with lines
EOF