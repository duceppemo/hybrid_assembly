#!/bin/perl

use strict;
use warnings;
use diagnostics;

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Sort::Naturally;
use Sort::Key::Natural qw(natkeysort);


####################
#                  #
#   Option flags   #
#                  #
####################


# Default values
my $in;
my $out;
my $width = 60;
my $to_sort;
my $help;
my $man;

GetOptions(
    'input|i=s' => \$in,  # string
    'output|o=s' => \$out,  # integer
    'width|w=i' => \$width,  # integer
    'sort|s' => \$to_sort,  # flag
    'help|h' => \$help,
    'man' => \$man) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;


##############
#            #
#   Checks   #
#            #
##############


#check if arguments are not empty
if (!(defined($in)) or !(defined($out)))
{
    pod2usage(1);
}

# Check if input file is a valid fasta file ???



########################
#                      #
#   Parse fasta file   #
#                      #
########################


open (my $in_fh, '<', $in) or die "Cannot open $in: $!\n";

my (%fastas, $id);  # data structure to store data
my $rank = 0;  # counter to keep track of original sequence order in fasta file

while (my $line = <$in_fh>)
{
    chomp($line);
    $line =~ s/^\s+|\s+$//g;  # trim remove white space from both ends of a string
    next unless length($line);  # skip empty lines or lines with only whitespace characters

    if ($line =~ /^>/)  # If the header line
    {
        $rank++;
        $id = $line;
    }
    else  # it's a sequence line
    {
        $fastas{$id}{'rank'} = $rank;
        push ( @{ $fastas{$id}{'sequence'} }, uc($line) );  # convert seq to uppercase
    }
}

close($in_fh);

# print Dumper(\%fastas);
# exit;


#########################
#                       #
#   Output fasta file   #
#                       #
#########################


sub format_fasta
{
    #pass arguments
    my ($i, $fh) = @_;

    my $seq = join('', @{ $fastas{$i}{'sequence'} });

    my $length = length($seq);
    my $position = 0;

    print {$fh} ("$i\n");  # header
    while ($position < $length)
    {
      my $part = substr($seq, $position, $width);
      print {$fh} ("$part\n");
      $position = $position + $width;
    }
}

# arrays for sorting naturally the headers
my @order_by_name;
foreach my $id (keys %fastas)  # sort by keys
{
    push (@order_by_name, $id);
}

@order_by_name = nsort(@order_by_name);

# Open output file
open (my $out_fh, '>', $out) or die "Cannot write to $out: $!\n";

if ($to_sort)
{
    # foreach my $rank (sort { $fastas{$a}{'id'} <=> $fastas{$b}{'id'} } keys %fastas)
    foreach my $id (@order_by_name)
    {
        # print ("$id\n");
        format_fasta($id, $out_fh);
    }
}
else
{
    # foreach my $rank (sort { $fastas{$a} <=> $fastas{$b} } keys %fastas)
    foreach my $id (sort { $fastas{$a}{'rank'} <=> $fastas{$b}{'rank'} } keys %fastas)
    {
        # print ("$id\n");
        format_fasta($id, $out_fh);
    }
}

close($out_fh);



__END__
=head1 NAME

formatFasta.pl - Format fasta files: Convert sequence to uppercase and column width to user preference (default 60).

=head1 SYNOPSIS

perl formatFasta.pl [options] -i in.fasta -o out.fasta

=head1 OPTIONS

=over 4

=item B<--help [-h]>

Print this help

=item B<--man>

Print script manual

=item B<--input [-i]>

Input fasta file. Mandatory

=item B<--output [-o]>

Output fasta file. Mandatory

=item B<--sort [-s]>

Sort entries by asccending number according to headers. Optional

=item B<--width [-w]>

Sequence line width. Default 60

=back

=head1 DESCRIPTION

B<formatFasta.pl> will convert sequence to uppercase and will set the sequence line width to user preference.
Default sequence line width is 60 base pairs. Headers are not modified.

=cut

