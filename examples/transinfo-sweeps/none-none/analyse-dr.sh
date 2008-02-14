
DEMO_DIR=$1

for TDDR_FILE in tddrs/*.tddr
do
  # Plot the tddr file.
  plt -n $TDDR_FILE > ${TDDR_FILE%.tddr}.png

  # Make a directory to hold the analysis of this tddr file.
  # For this test, we put these "histo" directories in the
  # tddr directory.
  HISTO_DIR=${TDDR_FILE%.tddr}.histos
  mkdir $HISTO_DIR

  # Take the tails of the columns of the tddr file
  # as samples of the equilibrium distribution,
  # according to the MCMC principle.
  $DEMO_DIR/tail-histos.pl $TDDR_FILE 200 5 $HISTO_DIR

  # Write gnuplot script that will plot the histograms and their
  # gamma approximating distributions.
  #
  # It should be noted that this script is the same everywhere it
  # occurs; it should be written once and for all at the top level.
  $DEMO_DIR/write-gp-doses-script.pl moleculizer-doses \
      $HISTO_DIR/plot-histograms.gp $HISTO_DIR

  # Generate the actual histogram files, one for each column
  # of the tddr file.
  for DATA_FILE in $HISTO_DIR/*.data
  do
    realbin 5 < $DATA_FILE > ${DATA_FILE%.data}.histo
  done

  # Generate the map from doses to response statistics.
  $DEMO_DIR/stats.pl moleculizer-doses $HISTO_DIR >> $HISTO_DIR/all-stats

  # Compute the dose/response transinformation.
  nu -load $DEMO_DIR/transinfo.scm \
      -stats-file $HISTO_DIR/all-stats \
      -transinfo-file $HISTO_DIR/transinfo

  # Generate data files for the gamma density plots done by the gnuplot
  # script written above.
  $DEMO_DIR/do-gamma-densities.pl $HISTO_DIR/all-stats $HISTO_DIR $DEMO_DIR

  # Run the gnuplot script written above.
  gnuplot $HISTO_DIR/plot-histograms.gp > $HISTO_DIR/histograms.png

  # Plot the dose response curve using the means of the response
  # distributions.
  #
  # Adding errorbars conveying the standard deviation
  # of the response distributions might be a nice touch.
  $DEMO_DIR/dr-curve.sh $HISTO_DIR
done