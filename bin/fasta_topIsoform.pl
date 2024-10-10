#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use lib '/opt/conda/lib/perl5/site_perl';

die "Please specify (1)fasta file , (2)top isoform \n" unless(@ARGV==2);

my $fastafile = $ARGV[0];
my $hitlist = $ARGV[1];
my $outfile="$fastafile\.nucl.longest.fa";

open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";

open(my $inhandle,  "<", $hitlist)   or die "Could not open $hitlist \n";



use Bio::SeqIO;
my $seqio = Bio::SeqIO->new(-file => "$fastafile", '-format' => 'Fasta');
my %fastadictionary=();
my @headersplit=();


while (my $seq = $seqio->next_seq){ ## selects one sequence at a time
   	## set variables for THIS sequence
    my $id = $seq->display_id;
	my $string = $seq->seq;
	$id=~ s/rna-//g;
    $id=~ s/transcript://g;
	#print $outhandle ">$id\n$string\n";
	$fastadictionary{$id}="$string";

	
}

my %HIT;
my $len_prots=0;
while ( my $line = <$inhandle> ){
	chomp $line;
	my @sp=split("\t", $line);
	if ($fastadictionary{$sp[1]}){
		#print "yes\";
		print $outhandle ">$sp[0]\n$fastadictionary{$sp[1]}\n";
	}
	else{
		print "No we lack $sp[1]\n";
	}
	$HIT{$sp[1]}="y";
	$len_prots++;
}




