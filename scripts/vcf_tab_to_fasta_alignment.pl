#!/usr/bin/perl

# Program to convert output of VCFtools' vcf-to-tab
# to FASTA alignment.

# Sample input file
#	$ head input.vcf.tab
#	chr10	94051	C	./	./	./	./	./	T/T
#	chr10	94056	T	./	./	./	./	./	C/C
#	chr10	94180	G	./	A/A	./	./	./	./


use strict;
use warnings;

my $exclude_het = 0;

my %iupac = (
			'G/' => 'G',
                        'C/' => 'C',
                        'T/' => 'T',
                        'A/' => 'A',
                        './' => '.',

	);

my $input_tab = shift;
chomp $input_tab;

open (TAB, "<$input_tab")
	or die "ERROR: Could not open input file $input_tab.\n";

my $header = <TAB>;

my @col_names = split /\t/, $header;

# Make temporary file with just lines we're going to use
my $temp_tab = $input_tab . "_clean";
open (TEMP, ">$temp_tab")
	or die "ERROR: Could not open temp file $temp_tab.\n";

# Get number of columns
my $num_cols = scalar @col_names;

LINE: foreach my $line (<TAB>) {

	my @data = split /\t/, $line;
	
	# Skip if this is indel (Length of @data will be less than $num_cols)
	if ((scalar @data) < $num_cols) {
		next LINE;
	}
	
	# Skip if any basepairs are actually 2 or more together
	for (my $i = 2; $i < $num_cols; $i++) {
		
		my $bp = $data[$i]; 
		chomp $bp;
		if ($bp =~ /\w{2,}/) {
			next LINE;
		}
	}

	if ($exclude_het) {
		# Exclude heterozygotes. Keep only fixed SNPs
		for (my $i = 2; $i < $num_cols; $i++) {
			
			my $bp = $data[$i]; 
			chomp $bp;
			if ($bp =~ /(\w)\/(\w)/) {
				if ($1 ne $2) {
					next LINE;
				}
			}
		}
	}
	
	# Otherwise write line to pure temporary file
	print TEMP $line;
}
	
close TAB;
close TEMP;

# Now convert cleaned tabular file to FASTA alignment

for (my $i = 3; $i < $num_cols; $i++) {

	my $ind = $col_names[$i];
	chomp $ind;
	
	print ">" . $ind . "\n";
	
	open (TEMP, "<$temp_tab")
		or die "ERROR: Could not open temp file $temp_tab.\n";

	# Count number of bp printed so far in this line
	my $count = 0;
	
	foreach my $line (<TEMP>) {
	
		my @data = split /\t/, $line;
		
		my $nuc = $data[$i];
		chomp $nuc;
		
		# Infer and print basepair. There are a few possibilities 
		
		# If we're reference, just print basepair
		if ($i == 2) {
			print $nuc;
			$count++;
		
		# Missing data
		} elsif ($nuc eq './') {
			print '-';
			$count++;
		
		# Data
		} elsif ($nuc =~ /(\w)\//) {
			my $first = $1;
			
			# Homozygote
				print $first;
				$count++;
			
		}
			
		if ($count == 100) {
			print "\n";
			$count = 0;
		}
	}
	
	close TEMP;
	
	print "\n";
}

exit;
