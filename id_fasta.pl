#!/opt/perl-5.16.2/bin/perl
#
# EPN, Wed Feb 18 09:28:21 2015
# This file was originally a single line, here is that line:
# perl -e 'while($line=<>){$l=substr($line,0,1);if($l eq ">"){$l=length($line)-1;$p=substr($line,1,$l);@acc=split(/\ /,$line);@ids=split(/\|/,$acc[0]);$a=@ids;if($a==4){print ">$ids[3] "}else{print ">$ids[1] "};print $p}else{print $line}}'

# Here it is rewritten in multiple lines, with added comments:

while ($line = <> ) {
  if($line =~ m/^\>/) { # starts with '>',  a header line
    chomp $line;
    $remainder = $line;
    $remainder =~ s/^\>//;
    $acc       = $remainder;
    $acc       =~ s/\ +.+$//; # remove everything including first space
    @ids = split(/\|/, $acc);
    if(scalar(@ids) >= 4) { # four '|' delimited tokens in $acc, rewrite ID as the 4th one
      print ">$ids[3] ";
    }
    else { # not >= four (probably 2) '|' delimited tokens in $acc, rewrite ID as the 2nd one
      print ">$ids[1] " 
    }
    print $remainder . "\n"; 
  }
  else { # not a header line, sequence, just print it
    print $line;
  }
}

