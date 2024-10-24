#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with orthofinder file, plus the GO file of interest\n\n";

my $in_orth="Orthogroups.tsv";
my $in_gofile=`ls *GO_format`;
chomp $in_gofile;
my @namesplit=split(/\./, $in_gofile);
my $species=$namesplit[0];

print "My species is $species\n\n";

open(my $filein, "<", $in_orth)   or die "Could not open $in_orth\n";

my $out="$species\_duplications.txt";
open(my $fileout, ">", $out)   or die "Could not open $out\n";

my $head=<$filein>;
chomp $head;
#print "$head";
my @splithead=split("\t",$head);
my $n=0;
my $sp_col;
foreach my $col (@splithead){
	my @orthnamesp=split(/\./, $col);
	my $species_orth=$orthnamesp[0];
	#print "$species_orth eq $species $n\n";
	if ($species_orth eq $species){
		$sp_col=$n;
	}
	$n++;
}

print "Our species col is $sp_col\n\n";

my $output=0;

while (my $line=<$filein>){
    chomp $line;
    my @split=split("\t", $line);
    my $orthogroup=$split[0];
    my $genes=$split[$sp_col];
    if ($genes){
    	if ($genes =~ m/, /g){
	    	#print "$genes\n";
	    	my @all=split(/, /, $genes);
	    	foreach my $bit (@all){
			#if the gene starts with a number, add a letter before it.
			if ($bit =~ /^\d/){
				$bit="GENE_$bit";
			}
	    		print $fileout "$bit\n";
	    		$output++;
	    	}
	    }
    }
}


if ($output){
	print "Now run the GO enrichment test\n\n";
	`ChopGO_VTS2.pl -i $species\_duplications.txt --GO_file $in_gofile`
}
else{
	print "No genes satisfied the cutoff, will not run orthofinder duplicate go\n";
}
