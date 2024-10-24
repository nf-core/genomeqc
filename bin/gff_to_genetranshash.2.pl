#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with *noquest.gff3 file\n\n";


my $in_gfffile=`ls *noquest.gff3`;
chomp $in_gfffile;
my @namesplit=split(/\./, $in_gfffile);
my $species=$namesplit[0];

print "My species is $species\n\n";

open(my $filein, "<", $in_gfffile)   or die "Could not open $in_gfffile\n";

my $out="$species\_gene_tran_list.txt";
open(my $fileout, ">", $out)   or die "Could not open $out\n";
my $out2="$species\_gene_alltran_list.txt";
open(my $fileout2, ">", $out2)   or die "Could not open $out2\n";
my $out3="$species\_longestisoform.txt";
open(my $fileout3, ">", $out3)   or die "Could not open $out3\n";

my %Gene_tran_hash;
my %longest;

while (my $line=<$filein>){
    chomp $line;
    my @split=split("\t", $line);
    #if line is an mRNA line, then we can check the gene to isoform IDs.
    my $gene;
    my $tran;
    #ignore blank lines, starting with #/
    if ($line =~ /^#/){
        #do nothing
    }
    else{
        #print "$line\n";
        if ($split[2] eq "mRNA"){
            my $length=$split[4]-$split[3];
       
            #Do different split if AUGUSTUS or NCBI
            if ($split[1] eq "AUGUSTUS"){
                #Its an AUGUSTUS GFF
                my @lsplit=split("\;", $split[8]);
                my @genesp=split("\=", $lsplit[1]);
                my @transp=split("\=", $lsplit[0]);
                $gene=$genesp[1];
                $tran=$transp[1];
            }
            elsif($split[1] eq "maker"){
                #Its probably a maker gff with ID and Parent for gene and transcript isoform names:
                my @lsplit=split("\;", $split[8]);
                my %temp_h;
                foreach my $bits (@lsplit){
                    my @spbit=split("\=", $bits);
                    $temp_h{$spbit[0]}=$spbit[1];
                }
                $gene=$temp_h{"Parent"};
                $tran=$temp_h{"ID"};
            }
            elsif($split[1] eq "Genbank"){
                #Its probably a Genbank gff with ID and Parent for gene and transcript isoform names:
                my @lsplit=split("\;", $split[8]);
                my %temp_h;
                foreach my $bits (@lsplit){
                    my @spbit=split("\=", $bits);
                    $temp_h{$spbit[0]}=$spbit[1];
                }
                my $curr_tran=$temp_h{"ID"};
                $gene=$temp_h{"Parent"};
                $tran=$temp_h{"ID"};
            }
            else{
                #Its probably a normal NCBI type:
                my @lsplit=split("\;", $split[8]);
                my %temp_h;
                foreach my $bits (@lsplit){
                    my @spbit=split("\=", $bits);
                    $temp_h{$spbit[0]}=$spbit[1];
                }
                my $fullgene=$temp_h{"Parent"};
                my @fullsp=split("\:", $fullgene);
                my @fullminusdash=split("\-",$fullsp[-1]);
                $gene=$fullminusdash[-1];
                $tran=$temp_h{"transcript_id"};
            }


            print $fileout "$gene\t$tran\n";
            if ($Gene_tran_hash{$gene}){
                my $old=$Gene_tran_hash{$gene};
                $Gene_tran_hash{$gene}="$old\,$tran";
            }
            else{
                $Gene_tran_hash{$gene}=$tran;
            }


            #Add to longest, if longest
            if ($longest{$gene}){
                my @old=split("\t", $longest{$gene});
                if ($length >= $old[1]){
                    $longest{$gene}="$tran\t$length";
                }
            }
            else{
                $longest{$gene}="$tran\t$length";
            }


        }
    }
}

foreach my $key (keys %Gene_tran_hash){
    print $fileout2 "$key\t$Gene_tran_hash{$key}\n";
}

foreach my $gen (keys %longest){
    my @spgen=split("\t", $longest{$gen});
    print $fileout3 "$gen\t$spgen[0]\n";
}

