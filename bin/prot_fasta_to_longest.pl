#!/usr/bin/perl
use warnings;
use strict;
use lib '/opt/conda/lib/perl5/site_perl';

die "Please specify (1)fasta file (2)longest isoform list\n" unless(@ARGV==2);

my $fastafile = $ARGV[0];
my $longest = $ARGV[1];

my $outfile="$fastafile\.largestIsoform.fa";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";

my $outfile2="$fastafile\.summary.txt";
open(my $outhandle2, ">", $outfile2)   or die "Could not open $outfile2 \n";

#Set counts for genes and transcripts. (TBD)
my $uniq_gene_count=0;
my $total_gene_count=0;

#Read in longest iso file and store in a hash:
open(my $inhandle, "<", $longest)   or die "Could not open $longest \n";
my %longest;
my %iso2gene;
while (my $line2=<$inhandle>){
    chomp $line2;
    my @split=split("\t", $line2);
    my $gene=$split[0];
    my $trans=$split[1];
    $trans=~ s/rna-//g;
    $trans=~ s/transcript://g;
    $longest{$trans}="longest";
    $iso2gene{$trans}=$gene;
    print "Trans: $trans equals Gene: $gene\n";
}

use Bio::SeqIO;
my $seqio = Bio::SeqIO->new(-file => "$fastafile", '-format' => 'Fasta');
my %fastadictionary=();
my @headersplit=();



while (my $seq = $seqio->next_seq){ ## selects one sequence at a time
   	## set variables for THIS sequence
    my $id = $seq->id;
    chomp $id;
    print "here $id\n";
    my $string = $seq->seq;
    chomp $string;
    my $new_id;
    my $len=length($string);
    #Remove NCBI rna- lines
    $id=~ s/rna-//g;
    $id=~ s/transcript://g;
    print "what $id\n";
    if ($longest{$id}){
        print $outhandle ">$iso2gene{$id}\n$string\n";
    }
	
	
}


=cut
#print off changes in gene to protein number.

print $outhandle2 "$total_gene_count\t$uniq_gene_count\n";

print "For file : $fastafile\n";
print "Total genes\tUnique proteins\n";
print "$total_gene_count\t$uniq_gene_count\n";

#print "Now print new fasta with one main protein per gene.\n";

foreach my $key ( sort keys %fastadictionary){
	if ($fastadictionary{$key} eq "Sequenceunavailable"){

	}
	else{
		print $outhandle ">$key\n$fastadictionary{$key}\n";
	}
	
}

#print "Finished:  input lines, output lines\n";
