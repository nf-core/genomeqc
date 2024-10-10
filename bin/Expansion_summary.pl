#!/usr/bin/perl
use warnings;
use strict;


print "Please be in folder with the *expansion.txt files\n\n";

#Find out all the names of GO terms found within our dataset.
my $ALL_go_terms=`cat *.txt | cut -f 1 | sort | uniq`;

#Download the GO iD to name hash:
`download_go_names.R > output`;

#Put go id to name infor into a hash.
my %ID_to_name;
my $input="GO_to_name";
open(my $goin, "<", $input)   or die "Could not open $input\n";
while (my $lineH=<$goin>){
	chomp $lineH;
	my @split=split("\t", $lineH);
	$ID_to_name{$split[1]}=$split[2];
	#print "$split[1] $split[2]\n";
}

#Read in file and infer species from file name
my @in_gofile=`ls *expansions.txt`;

my %species_go_hash;
foreach my $species_files (@in_gofile){
	chomp $species_files;
	my @name=split(/\./, $species_files);
	my $species=$name[0];
	open(my $filein, "<", $species_files)   or die "Could not open $species_files\n";
	while (my $line=<$filein>){
    	chomp $line;
    	my @split=split("\t", $line);
    	$species_go_hash{$species}{$split[0]}=$split[1];
	}
	close $filein;
}

#Set the full output name
my $out="GO_table_counts.tsv";
open(my $fileout, ">", $out)   or die "Could not open $out\n";

#Set the output name
my $out2="GO_table_counts_forCAFE.tsv";
open(my $fileout2, ">", $out2)   or die "Could not open $out2\n";


#Print header:
print $fileout "Desc\tFamily ID";
print $fileout2 "Desc\tFamily ID";
foreach my $species (keys %species_go_hash){
	print $fileout "\t$species";
	print $fileout2 "\t$species";
}
print $fileout "\n";
print $fileout2 "\n";

my %line_to_score;  #This hash we keep to we can print off the top 300 most variable terms.

#Print table contents
foreach my $GO_terms (keys %ID_to_name){
	print $fileout "$ID_to_name{$GO_terms}\t$GO_terms";
	my $count=0;
	my $sum=0;
	my $max=0; #Start the max at zero.
	my $min=100000000000; # We set min super high, so that the next number should be lower than this.
	my @row_info_for_hash=();
	#Now for each species fill in the number of genes associated with each term.
	foreach my $species (keys %species_go_hash){
		if (exists $species_go_hash{$species}{$GO_terms}){
			print $fileout "\t$species_go_hash{$species}{$GO_terms}";
			$count++;
			$sum+=$species_go_hash{$species}{$GO_terms};
			if ($species_go_hash{$species}{$GO_terms} > $max){
				$max=$species_go_hash{$species}{$GO_terms};
			}
			if ($species_go_hash{$species}{$GO_terms} < $min){
				$min=$species_go_hash{$species}{$GO_terms};
			}
			push(@row_info_for_hash, $species_go_hash{$species}{$GO_terms});
		}
		else{
			print $fileout "\t0";
			$count++;
			$min=0;
			push(@row_info_for_hash, "0");
		}
	}

	my $aver=$sum/$count;
	my $range=$max-$min;
	my $score;
	if ($aver == 0){
		$score=0;
	}
	elsif ($range == 0){
		$score=0;
	}
	else{
		$score=$range/$aver;
	}
	
	print $fileout "\t$score\n";

	my $join=join("\t", @row_info_for_hash);
	my $full_line="$ID_to_name{$GO_terms}\t$GO_terms\t$join";
	$line_to_score{$full_line}=$score;
}

#Now print off the top hits for cafe:
my $top300=300;
foreach my $key (sort {$line_to_score{$b} <=> $line_to_score{$a}} keys %line_to_score) {
	if ($top300){
		print $fileout2 "$key\n";
		$top300--;
	}
}

print "Script complete\n\n";

