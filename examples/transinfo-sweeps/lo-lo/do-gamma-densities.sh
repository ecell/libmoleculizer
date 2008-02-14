
for stats_file in histos/*.stats
do
  nu -load gamma-densities.scm < $stats_file > ${stats_file%.stats}.gamma
done