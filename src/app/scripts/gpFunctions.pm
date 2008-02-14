
package gpFunctions;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(writeScript writeMultiScript);

use strict;

sub writePlot {
  my $scriptHandle = shift @_;
  my $dataFileName = shift @_;

  # Read the header line out of the data file and parse it
  # into an array of headers.
  open dataFile, $dataFileName
    or die "Couldn't open data file $dataFileName.\n";
  my @headers = split /\s+/, <dataFile>;
  close dataFile;

  # Make a an associative array mapping headers to gnuplot
  # column numbers.
  my %hmap;
  my $colNum = 0;
  while($colNum <= $#headers) {
    my $gnuplotColNum = $colNum + 1;
    $hmap{$headers[$colNum]} = $gnuplotColNum;
    $colNum++;
  }

  # If there are any more arguments, they should be the headers
  # whose columns we want to plot.  Otherwise, all columns are
  # plotted.
  my @plotHeaders;
  if(@_) {
    @plotHeaders = @_;
  }
  else {
    # Omitting the #time header.
    shift @headers;
    @plotHeaders = @headers;
  }

  # Emit the (multiline) gnuplot  plot command.
  print $scriptHandle "plot \\\n";
  while(@plotHeaders) {
    my $headerName = shift @plotHeaders;
    my $gnuplotColNum = $hmap{$headerName};

    if($gnuplotColNum <= 0) {
      die "Couldn't find column header `$headerName'.\n";
    }
    
    print $scriptHandle "\"$dataFileName\" ";
    print $scriptHandle "using 1:$gnuplotColNum ";
    print $scriptHandle "title \"$headerName\"";

    # All lines except the last have to have a comma and
    # escaped newline.
    if(@plotHeaders) {
      print $scriptHandle ",\\\n";
    }
    else {
      print $scriptHandle "\n";
    }
  }
}

# This writes a plot statement to plot the same column from several
# files, instead of several columns from one file.  Subsume these into
# a single function that will do arbitrary file/column combinations???
#
# This is also complicated by not wanting to use the full file name
# as a title for the curve that comes from it, but rather a run number
# of some sort.
sub writeMultiPlot {
  my $scriptHandle = shift @_;
  my $header = shift @_;
  my @dataFileNames = @{shift @_};
  my @titles = @{shift @_};

  # To convert the header into a column number, we use the header
  # line out of the first file given.
  open dataFile, $dataFileNames[0]
    or die "writeMultiPlot: couldn't open data file `$dataFileNames[0]'.\n";
  my @headers = split /\s+/, <dataFile>;
  close dataFile;

  # Which (gnuplot) column is the header in?  Gnuplot numbers columns
  # starting with 1.
  my $gnuplotColNum = 0;
  for(my $colNdx = 0;
      $colNdx <= $#headers;
      $colNdx++)
      {
	if("$header" eq "$headers[$colNdx]")
	    {
	      $gnuplotColNum = $colNdx + 1;
	    }
      }
  if($gnuplotColNum == 0)
      {
	die "writeMultiPlot: couldn't find header `$header' in "
	  . "data file `$dataFileNames[0]'.\n";
      }

  # Emit the gnuplot plot command.
  print $scriptHandle "plot \\\n";
  while(@dataFileNames) {
    my $dataFileName = shift @dataFileNames;
    my $title = shift @titles;
    print $scriptHandle "\"$dataFileName\" ";
    print $scriptHandle "using 1:$gnuplotColNum ";
    print $scriptHandle "title \"$title\"";

    if(@dataFileNames) {
      print $scriptHandle ",\\\n";
    }
    else {
      print $scriptHandle "\n";
    }
  }
}

# Writes the gnuplot "prelude" to stdout.  This gets the general sort
# of plot that I prefer.
sub writeBaseStyle {
  my ($scriptHandle, $title, $yLabel, $legendTitle) = @_;

  print $scriptHandle <<EOF;
set key outside width 1 title "$legendTitle" box
set title "$title"
set xlabel "Simulation time in seconds"
set ylabel "$yLabel"
set data style lines
EOF
}

# Makes the plot semi-log.
sub writeLogStyle {
  my $scriptHandle = shift;
  print $scriptHandle "set logscale y\n";
}

# Causes gnuplot to emit PostScript instead of displaying to the
# screen.
sub writePostScriptStyle {
  my $scriptHandle = shift;
  print $scriptHandle "set terminal postscript color solid \"Times-Roman\"\n";
}

# Causes gnuplot to emit Fig output instead of displaying to the screen.
sub writeFigStyle {
  my $scriptHandle = shift;
  print $scriptHandle "set terminal fig color big\n";
}

# Causes gnuplot to emit png output, instead of displaying to the screen.
sub writePngStyle {
  my $scriptHandle = shift;
  print $scriptHandle "set terminal png\n";
}

# Tells gnuplot to exit this script.
sub writeExit {
  my $scriptHandle = shift;
  print $scriptHandle "exit\n";
}

# Tells gnuplot to reload this script, after a pause.
sub writeReload {
  my $scriptHandle = shift;
  print $scriptHandle "pause 5\n";
  print $scriptHandle "reread\n";
}

# Writes the whole gnuplot script based on options passed in
# through a hash.\
sub writeScript {
  my $scriptHandle = shift;
  my %scriptOptions = @_;
  
  writeBaseStyle($scriptHandle,
		 $scriptOptions{"plotTitle"},
		 $scriptOptions{"yLabel"},
		 $scriptOptions{"legendTitle"});

  if($scriptOptions{"logPlot"}) {
    writeLogStyle($scriptHandle);
  }

  if($scriptOptions{"postScript"}) {
    writePostScriptStyle($scriptHandle);
  }

  if($scriptOptions{"fig"}) {
    writeFigStyle($scriptHandle);
  }

  if($scriptOptions{"png"}) {
    writePngStyle($scriptHandle);
  }

  my @headersToPlot = @{$scriptOptions{"headersToPlot"}};
  writePlot($scriptHandle,
	    $scriptOptions{"dataFileName"},
	    @headersToPlot);

  if($scriptOptions{"reload"}) {
    writeReload($scriptHandle);
  }
  else {
    writeExit($scriptHandle);
  }
}

sub writeMultiScript {
  my $scriptHandle = shift;
  my %multiScriptOptions = @_;

  writeBaseStyle($scriptHandle,
		 $multiScriptOptions{"plotTitle"},
		 $multiScriptOptions{"yLabel"},
		 $multiScriptOptions{"legendTitle"});

  if($multiScriptOptions{"logPlot"}) {
    writeLogStyle($scriptHandle);
  }

  if($multiScriptOptions{"postScript"}) {
    writePostScriptStyle($scriptHandle);
  }

  if($multiScriptOptions{"fig"}) {
    writeFigStyle($scriptHandle);
  }

  if($multiScriptOptions{"png"}) {
    writePngStyle($scriptHandle);
  }

  writeMultiPlot($scriptHandle,
		 $multiScriptOptions{"header"},
		 $multiScriptOptions{"dataFileNames"},
		 $multiScriptOptions{"lineTitles"});

  writeExit($scriptHandle);
}


