#!/bin/perl

#########################################################################################
#											#
#	Input: List of Entrez UIDs (integer identifiers, e.g. PMID, GI, Gene ID)	#
#	EFetch Output: Formatted data records (e.g. abstracts, FASTA)			#
#											#
#########################################################################################

use LWP::Simple;
use Term::ProgressBar;
use strict;
use warnings;
use diagnostics;
#use Data::Dumper;

#I/O files
my $input = $ARGV[0];
my $output = $ARGV[1];

#Check required arguments
if (scalar(@ARGV) != 2)
{
    print "Usage: perl http_post.pl <input.list> <output.fasta>\n";
    exit;
}

#Open the input file, in read-only mode
open(my $infilehandle, "<", $input) or die "Error opening $input: $!\n";

#parse list to fetch 200 sequences at the time
my @query;
my $counter = 0;
my $string;

#create array with file
while(my $line = <$infilehandle>)
{
	$counter += 1;
	chomp($line);
	if ($counter < 200)
	{
		$string .= $line . ",";
	}
	else
	{
		$string .= $line;
		push(@query, $string);
		$counter = 0;
		$string = "";
	}
	
}

#remainig IDs
#delete the last comma that has no ID following
chop($string);
push(@query, $string);

#print Dumper(\@query);

#Close input file
close ($infilehandle);

my $db = 'protein';

#open output file for writing
open(OUT, ">>", $output) || die "Can't open file!\n";

#progress bar
print "\n\nDownloading\n\n";
my $nQuery = scalar(@query);
my $progress = Term::ProgressBar->new($nQuery);

#retrieve data in batches of 200
foreach my $entry (@query)
{
	#my $efetch_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&id=$entry&rettype=gp&retmode=text";# GenPept
	my $efetch_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&id=$entry&rettype=fasta&retmode=text";# fasta
	#my $efetch_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&id=$entry";#&rettype=fasta&retmode=text";# ASN.1
	my $efetch_out = get($efetch_url);
	print OUT "$efetch_out";
	$progress->update($_);
}
close OUT;

print ("\n\ndone\n");
