#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with the GO file of interest\n\n";

#Read in file and infer species from file name
my $in_gofile=`ls *GO_format`;
chomp $in_gofile;
my @namesplit=split(/\./, $in_gofile);
my $species=$namesplit[0];
open(my $filein, "<", $in_gofile)   or die "Could not open $in_gofile\n";

print "My species is $species\n\n";

#Set the output name
my $out="$species\.go_family_expansions.txt";
open(my $fileout, ">", $out)   or die "Could not open $out\n";

my %uniq;
my %go_count;
while (my $line=<$filein>){
    chomp $line;
    my @split=split("\t", $line);

    if ($uniq{$line}){
    	#do nothing, we already read in this line
    }
    else{
    	$go_count{$split[1]}++;
    	$uniq{$line}="EXISTS";
    } 
}

print "Now print out the GO counts\n\n";

foreach my $GO_term (sort keys %go_count){
    print $fileout "$GO_term\t$go_count{$GO_term}\n";
}

