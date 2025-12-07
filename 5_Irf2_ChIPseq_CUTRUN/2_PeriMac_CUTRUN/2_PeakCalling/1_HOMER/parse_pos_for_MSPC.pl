#! /usr/bin/perl

use warnings;
use strict;

# a script to extract a gene list from big files

unless (scalar @ARGV == 1) {
	die "\nNot enough command line arguments.\n".
	"Usage : parse_pos_for_MSPC.pl <HOMER pos file>.\n";

}


# file that contains peak names
my $file = shift @ARGV;

#get base name to create new file name
$file =~ /(.+)\.pos/;

# use captured pattern to create output file name
my $out = $1.".cl.bed";

# open pos file for reading
unless (open(FILE,$file)) {
	die "cannot open $file file.\n";
}

#open output file for writing
unless (open(OUTPUT, '>', $out)) {
	die "cannot create $out file. \n";
}

while (my $line = <FILE>) {

	next if ($line =~ m/^#/g);

	# parse tab separated files
	my @dat = get_line_data ($line);
	
	# calculate -log10 p-value for MSPC [p-value vs Local field]
	my $mlog10p = -1 * log10($dat[9]);
	
	# get new bed formatted line
	my $newline = "$dat[1]\t$dat[2]\t$dat[3]\tHOMER_Peak_$.\t$mlog10p\n";
	
	print OUTPUT "$newline";

};

close FILE;
close OUTPUT;
	 
exit;

###########################################################
# SUBROUTINES
###########################################################

###########################################################
# a subroutine that separates fields from a data line and
# returns them in an array

sub get_line_data {

    my $line = $_[0];
    
    chomp $line;  

    my @linedata = split(/\t/, $line);
        
    return @linedata;
}

###########################################################
# a subroutine that calculates log10 values
sub log10 {
    my $n = shift;
    return log($n + 1e-100)/log(10); # add small number
}