use strict;
use warnings;

my $usage = "perl test_extract_fasta_multi_exon.pl <executable-dir> <testfile-dir> <'1' for verbose mode, '0' for quiet mode>\n";
if(scalar(@ARGV) != 3) { die $usage; }

my ($execdir, $testdir, $be_verbose) = @ARGV;

if($be_verbose ne "1" && $be_verbose ne "0") { die $usage; }

$execdir =~ s/\/$//; # remove trailing '/' if one exists
$testdir =~ s/\/$//; # remove trailing '/' if one exists

# definitions of variables
my $fail_ntests       = 11; # number of 'fail' tests,    the tests that are     expected to return a failure (non-zero exit status)
my $nofail_ntests     = 13; # number of 'no fail' tests, the tests that are not expected to return a failure (non-zero exit status)
my @fail_descA        = (); # brief descriptions of each fail test, for informative output
my @nofail_descA      = (); # brief descriptions of each no fail test, for informative output
my @fail_inputA       = (); # input files for each fail test
my @nofail_inputA     = (); # input files for each no fail test
my @nofail_outputA    = (); # output files for each no fail test
my @fail_nerrlinesA   = (); # number of expected stderr output lines for fail tests
my @fail_stdout_okayA = (); # '1' if it's okay for some stdout to be printed for each fail test
my @nofail_noutlinesA = (); # number of expected stdout output lines for no fail tests
my ($f, $n);                # counter over fail and no fail tests, respectively

# check executable files exist
my $extract_fasta = $execdir . "/extract_fasta_multi_exon";
my $id_fasta      = $execdir . "/id_fasta.pl";
if(! -e $extract_fasta) { die "ERROR extract_fasta does not exist in $execdir"; }
if(! -e $id_fasta)      { die "ERROR id_fasta does not exist in $execdir"; }

# Description of tests, included here for reference
# tests that are expected to fail (testing the error checking code in the program)
$fail_descA[0]  = "wrong number of pieces";
$fail_descA[1]  = "wrong order of pieces";
$fail_descA[2]  = "start > end";
$fail_descA[3]  = "too many pieces";
$fail_descA[4]  = "no piece start";
$fail_descA[5]  = "no piece end";
$fail_descA[6]  = "no strand";
$fail_descA[7]  = "not the specified number of pieces (too many)";
$fail_descA[8]  = "not the specified number of pieces (too few)";
$fail_descA[9] = "token after strand";
$fail_descA[10] = "end position exceeds sequence length";
# make sure we have descriptions for all fail tests
for($f = 0; $f < $fail_ntests; $f++) { 
  if((! defined $fail_descA[$f]) || ($fail_descA[$f] eq "")) { 
    die "ERROR (in test script) fail test " . ($f+1) . " does not have a description";
  }
}

# tests that are expected to not fail
$nofail_descA[0]  = "fetching entire sequences";
$nofail_descA[1]  = "fetching single piece sequences";
$nofail_descA[2]  = "fetching a multi-piece sequence on the positive strand";
$nofail_descA[3]  = "fetching a multi-piece sequence on the negative strand";
$nofail_descA[4]  = "overlapping intervals (single piece, same strand)";
$nofail_descA[5]  = "overlapping intervals (single piece, opposite strand)";
$nofail_descA[6]  = "overlapping intervals (multi piece, same strand)";
$nofail_descA[7]  = "overlapping intervals (multi piece, opposite strand)";
$nofail_descA[8]  = "subsumed intervals (single piece, same strand)";
$nofail_descA[9]  = "subsumed intervals (single piece, opposite strand)";
$nofail_descA[10] = "subsumed intervals (multi piece, same strand)";
$nofail_descA[11] = "subsumed intervals (multi piece, opposite strand)";
$nofail_descA[12] = "extra token";
# make sure we have descriptions for all no fail tests
for($n = 0; $n < $nofail_ntests; $n++) { 
  if((! defined $nofail_descA[$n]) || ($nofail_descA[$n] eq "")) { 
    die "ERROR (in test script) nofail test " . ($n+1) . " does not have a description";
  }
}

# check that required inputs/outputs exist
for($f = 0; $f < $fail_ntests; $f++) { 
  $fail_inputA[$f]  = $testdir . "/fail." . ($f+1) . ".in"; # off-by-one
  if(! -e $fail_inputA[$f]) { die "ERROR " . $fail_inputA[$f] . "(input file) does not exist"; }
  $fail_descA[$f] .= " [input file: " . $fail_inputA[$f] . "]";
}
for($n = 0; $n < $nofail_ntests; $n++) { 
  $nofail_inputA[$n]   = $testdir . "/nofail." . ($n+1) . ".in";  # off-by-one
  $nofail_outputA[$n]  = $testdir . "/nofail." . ($n+1) . ".out"; # off-by-one
  if(! -e $nofail_inputA[$n])  { die "ERROR " . $nofail_inputA[$n] . " (input file) does not exist"; } 
  if(! -s $nofail_outputA[$n]) { die "ERROR " . $nofail_outputA[$n] . " (output file) does not exist"; }
  $nofail_descA[$n] .= " [input file: " . $nofail_inputA[$n] . ", output_file: " . $nofail_outputA[$n] . "]";
}

# set default number of stderr lines for fail tests,
# and whether or not it's okay for each fail test to
# create some stdout, or not
for($f = 0; $f < $fail_ntests; $f++) { 
  $fail_nerrlinesA[$f] = 1;
  $fail_stdout_okayA[$f] = 0; # not okay, by default
}
# set exceptions: 
$fail_stdout_okayA[10] = 1; # this test should print some to stdout before failing

$nofail_noutlinesA[0]  = 187;
$nofail_noutlinesA[1]  = 6;
$nofail_noutlinesA[2]  = 5;
$nofail_noutlinesA[3]  = 5;
$nofail_noutlinesA[4]  = 4;
$nofail_noutlinesA[5]  = 4;
$nofail_noutlinesA[6]  = 4;
$nofail_noutlinesA[7]  = 4;
$nofail_noutlinesA[8]  = 4;
$nofail_noutlinesA[9]  = 4;
$nofail_noutlinesA[10] = 4;
$nofail_noutlinesA[11] = 4;
$nofail_noutlinesA[12] = 4;

# Run fail tests (these should fail and give no stdout)
for($f = 0; $f < $fail_ntests; $f++) { 
  my $stderr = "fail." . ($f+1) . ".stderr";
  my $stdout = "fail." . ($f+1) . ".stdout";

  # first create the idfetch input file
  my $idfetch_in = "fail." . ($f+1) . ".idfetch.in";
  my $cmd = "cat $fail_inputA[$f] | awk '{ print \$1 }' | sort | uniq > $idfetch_in";
  RunCommand($cmd, 0, "idfetch input creation", 0);

  # run command and ensure it has expected non-zero exit status (die if it does not)
  $cmd = "idfetch -t 5 -c 1 -G $idfetch_in | $id_fasta | $extract_fasta $fail_inputA[$f] 2> $stderr > $stdout";
  RunCommand($cmd, 1, $fail_descA[$f], $be_verbose); # '1' says: failure is expected

  # check that command had expected number of lines in stderr and 
  # either did or did not have stdout output
  CheckNumLinesInFile($stderr, $fail_nerrlinesA[$f], $fail_descA[$f]);
  if(! $fail_stdout_okayA[$f]) { 
    CheckNumLinesInFile($stdout, 0, $fail_descA[$f]); # 0 says that $stdout should either not exist or exist but with 0 lines
  }

  for my $file ($stderr, $stdout, $idfetch_in) { 
    if(-e $file) { unlink $file; }
  }
}

# Run no fail tests (these should not fail and should give output)
for($n = 0; $n < $nofail_ntests; $n++) { 
  my $stderr = "nofail." . ($n+1) . ".stderr"; # off-by-one, to match input file name
  my $stdout = "nofail." . ($n+1) . ".stdout"; # off-by-one, to match input file name

  # first create the idfetch input file
  my $idfetch_in = "nofail." . ($n+1) . ".idfetch.in";
  my $cmd = "cat $nofail_inputA[$n] | awk '{ print \$1 }' | sort | uniq > $idfetch_in";
  RunCommand($cmd, 0, "idfetch input creation", 0);

  # run command and ensure it has expected non-zero exit status (die if it does not)
  $cmd = "idfetch -t 5 -c 1 -G $idfetch_in | $id_fasta | $extract_fasta $nofail_inputA[$n] 2> $stderr > $stdout";
  RunCommand($cmd, 0, $nofail_descA[$n], $be_verbose); # '0' says: failure is NOT expected

  # check that command had expected stderr and stdout output
  CheckNumLinesInFile($stderr, 0,                      $nofail_descA[$n]); # 0 says that $stderr should either not exist or exist but with 0 lines 
  CheckNumLinesInFile($stdout, $nofail_noutlinesA[$n], $nofail_descA[$n]);

  # make sure output file is what we expect
  CompareFilesWithDiff($stdout, $nofail_outputA[$n], $nofail_descA[$n]);

  for my $file ($stderr, $stdout, $idfetch_in) { 
    if(-e $file) { unlink $file; }
  }
}
print "PASS [$fail_ntests tests failed as expected and $nofail_ntests tests succeeded as expected]\n";

exit 0;  


# Subroutine: RunCommand()
# Args:       $cmd:            command to run, with a "system" command;
#             $expect_failure: '1' if $cmd should fail (return non-zero exit status)
#             $desc:           description to print if command return expected exit status
#             $be_verbose:     '1' to output command before we run it, '0' not to
# Dies:       if $cmd fails and $expect_failure == 0
#             OR
#             if $cmd does not fail and $expect_failure == 1

sub RunCommand {
  my $sub_name = "RunCommand()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cmd, $expect_failure, $desc, $be_verbose) = @_;

  if($be_verbose) { printf("Test: $desc\n"); }   
  #print Running cmd: $cmd (expect failure: $expect_failure)\n", $desc);
  system($cmd);
  if($expect_failure) { 
    if($? == 0) { die "ERROR command with desc $desc did not fail.\nActual command:\n$cmd\n"; }
  }
  else { 
    if($? != 0) { die "ERROR command with desc $desc failed.\nActual command:\n$cmd\n"; }
  }
  return;
}

# Subroutine: CheckNumLinesInFile()
# Args:       $file:       command to run, with a "system" command;
#             $nlines_exp: number of lines expected
#                          if '0',  file can either not exist or exist but be empty. 
#                          if '-1', file should not exist
#             $desc:       description to print if command return expected exit status
# Dies:       if file does not have expected number of lines

sub CheckNumLinesInFile {
  my $sub_name = "CheckNumLinesInFile()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($file, $nlines_exp, $desc) = @_;
  if($nlines_exp < -1) { die "ERROR $sub_name, negative number of expected lines"; }

  my $die_msg = "This means the following test failed: $desc\n";

  if($nlines_exp == -1) { 
    if(-e $file) { die "ERROR file $file exists (it should not exist).\n$die_msg"; }
  }
  elsif($nlines_exp == 0) { 
    if(-s $file) { die "ERROR file $file exists and is non-empty (it should either not exist or be empty).\n$die_msg"; }
  }
  else { # $nlines_exp > 0
    if(! -e $file) { die "ERROR file $file does not exist but it should.\n$die_msg"; }
    my $nlines = `wc -l $file | awk '{ print \$1 }'`;
    chomp $nlines;
    if($? != 0) { die "ERROR problem getting number of lines in file $file.\n$die_msg"; }
    if($nlines != $nlines_exp) { die "ERROR, wrong number of lines ($nlines != $nlines_exp) in file $file.\n$die_msg"; }
  }
  return;
}

# Subroutine: CompareFilesWithDiff()
# Args:       $file1: first file
#             $file2: second file
#             $desc:  description to print if files differ
# Dies:       if 'diff' says $file1 and $file2 are different

sub CompareFilesWithDiff {
  my $sub_name = "CompareFilesWithDiff()";
  my $nargs_exp = 3;

  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($file1, $file2, $desc) = @_;

  my $die_msg = "This means the following test failed: $desc\n";

  my $diff_output = `diff $file1 $file2`;

  if($? != 0 || $diff_output =~ m/\w/) { 
    die "ERROR diff says $file1 and $file2 differ.\n$die_msg" . $diff_output; 
  }
  return;
}
  

