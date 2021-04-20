#! /usr/bin/perl

$jobFile = "job.dat";
$machineFile = "machines.dat";

open (IN, "<$jobFile") or die "No job file.\n";
foreach $line (<IN>) {
  $readLine = $line;
}
close (IN);
$jobID = substr($readLine, -8);
print {STDOUT} "$jobID";

@nodes = getJobNodes($jobID);

open (OUT,">$machineFile");
print OUT "test";
foreach $node (@nodes) {
    if ($node =~ /:ppn=(\d+)/) { #node matches :ppn=(\any number of digits)
	my $ppn = $1;
	$node =~ s/:.+//;  #removes : and 1 or more of any character
	for (my $q = 0;$q < $ppn;$q++) {  #while q < ppn; loop through node and write to file
	    print OUT "$node:1\n";
	}
    }
    else {
	print OUT "$node:1\n";
    }
}
close (OUT);

sub getJobNodes {
    my ($jobID) = @_;
    my @nodes;
    my @jobInfoLines = split /\n/,`qstat -n -1 $jobID`;
    my $jobLine = substr $jobInfoLines[-1],87;
    $jobLine =~ s/\s//g;    #removes all whitespace
    $jobLine =~ s/\/\d*\+/\+/g;   #replaces "/any number of digits" with "\+"
    $jobLine =~ s/\/\d*$//;  #removes the first "?"
    @nodes = split /\+/,$jobLine;  #nodes now is a list of jobline split on "\+"
    return @nodes;
}
