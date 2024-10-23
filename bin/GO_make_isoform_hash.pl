#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with *_gene_alltran_list.txt file and the GO protein file *_Result_All_Combine_GO_format\n\n";


my $in_gofile=`ls *_Result_All_Combine_GO_format`;
chomp $in_gofile;
my @namesplit=split(/\./, $in_gofile);
my $species=$namesplit[0];

my $gene2tran="$species\_gene_alltran_list.txt";

print "My species is $species\n\n";

open(my $filein, "<", $in_gofile)   or die "Could not open $in_gofile\n";

my $out="$species\.transcripts_Combine_GO_format.txt";
open(my $fileout, ">", $out)   or die "Could not open $out\n";


my %GO_hash;

while (my $line=<$filein>){
    chomp $line;
    my @split=split("\t", $line);
    if ($GO_hash{$split[0]}){
    	my $old=$GO_hash{$split[0]};
    	$GO_hash{$split[0]}="$old\,$split[1]";
    }
    else{
    	$GO_hash{$split[0]}=$split[1];
    }
}

open(my $tranin, "<", $gene2tran)   or die "Could not open $gene2tran\n";

while (my $line2=<$tranin>){
    chomp $line2;
    my @split=split("\t", $line2);
    my $gene=$split[0];
    my @trans=split(/\,/, $split[1]);
    foreach my $transcripts (@trans){
    	if ($GO_hash{$gene}){
    		my @go_terms_of_gene=split(/\,/, $GO_hash{$gene});
    		foreach my $GO_term_of_gene (@go_terms_of_gene){
    			print $fileout "$transcripts\t$GO_term_of_gene\n";
    		}
    	}
    }
}

