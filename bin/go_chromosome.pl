#!/usr/bin/perl
use warnings;
use strict;

print "Please be in folder with focal gff3 file and GO hashes\n\n";

my @goes=`ls *Result_All_Combine_GO_format`;
my @in_gfffile=`ls *.gff3`;
my $ortho="Orthogroups.tsv";


#Store orthogroup names hash/
my %orthogroup_hash;
open(my $orthin, "<", $ortho)   or die "Could not open $ortho\n";
my $header=<$orthin>;
chomp $header;
my @colHeadsplit=split("\t", $header);
while (my $lineOrtho=<$orthin>){
    chomp $lineOrtho;
    my $n=1;
    my @colsplit=split("\t", $lineOrtho);
    my $OG= shift(@colsplit);
    foreach my $row (@colsplit){
        my @isoform_sp=split(", ", $row);

        #print "HERE $OG : $n\n";
        my @big_name=split(/\./, $colHeadsplit[$n]);
        my $tidyname=$big_name[0];
        foreach my $iso (@isoform_sp){
            $orthogroup_hash{$tidyname}{$iso}=$OG;
        }
        $n++;
    }
}



my @jobs;   #To store pairs of gff and go files to run

#Check we have the right files for each species:
foreach my $gofile (@goes){
    chomp $gofile;
    #print "$gofile\n";
    my @sp_gp =split(/\./, $gofile);
    my $match=0;
    foreach my $gfffile (@in_gfffile){
        chomp $gfffile;
        my @sp_gff=split(/\./, $gfffile);
        #print "$gfffile cat $sp_gp[0]  and $sp_gff[0]\n";
        

        if ($sp_gp[0] eq $sp_gff[0]){
            $match=1;
            push (@jobs, "$gofile $gfffile");
        }
    }
    if ($match){
        #print "$sp_gp[0] is matched\n";
    }
    else{
        print "ERROR: $sp_gp[0] does not have an equivalent gff file\n"
    }
}




#Now run through the jobs and prepare the input files.
my %Gene_tran_hash;

foreach my $species (@jobs){
    my @sp=split(/\ /, $species);
    my $go=$sp[0];
    my $gff=$sp[1];
    chomp $gff;
    
    my @go_split =split(/\./, $go);
    my $species_name=$go_split[0];

    my $out="$species_name\.go_r_file.txt";
    open(my $fileout, ">", $out)   or die "Could not open $out\n";

    my $out2="$species_name\.go_r_file.noDuplicates.txt";
    open(my $fileout2, ">", $out2)   or die "Could not open $out2\n";

    open(my $filein, "<", $gff)   or die "Could not open $gff\n";
    while (my $line=<$filein>){
        chomp $line;
        my @split=split("\t", $line);
        #if line is an mRNA line, then we can check the gene to isoform IDs.
        my $gene;
        my $tran;
        my $scaffold=$split[0];

        #ignore blank lines, starting with #/
        if ($line =~ /^#/){
            #do nothing
        }
        else{
            #print "$line\n";
            if ($split[2] eq "mRNA"){
           
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


                print $fileout  "$gene\t$scaffold\n";
                #if the species gene exists in the orthofinder table, then print it to the OG run:
                if ($orthogroup_hash{$species_name}{$gene}){
                    print $fileout2 "$orthogroup_hash{$species_name}{$gene}\t$scaffold\n";
                }
                
                if ($Gene_tran_hash{$gene}){
                    my $old=$Gene_tran_hash{$gene};
                    $Gene_tran_hash{$gene}="$old\,$tran";
                }
                else{
                    $Gene_tran_hash{$gene}=$tran;
                }
            }
        }
    }

    #Now make a Orthogroup background file:

    open(my $filego, "<", $go)   or die "Could not open $go\n";

    my $out3="$species_name\.go_r_file.noDuplicates_BK.txt";
    open(my $fileout3, ">", $out3)   or die "Could not open $out3\n";

    my %exist_hit;
    while (my $linego=<$filego>){
        chomp $linego;
        my @splitgo=split("\t", $linego);

        if ($orthogroup_hash{$species_name}{$splitgo[0]}){
            if ($exist_hit{"$orthogroup_hash{$species_name}{$splitgo[0]}\t$splitgo[1]"}){
                #do nothing; hit exists.
            }
            else{
                print $fileout3 "$orthogroup_hash{$species_name}{$splitgo[0]}\t$splitgo[1]\n";
                $exist_hit{"$orthogroup_hash{$species_name}{$splitgo[0]}\t$splitgo[1]"}="YES";
            }
        }
    }


    #close handles, sort duplicated go list, and run ChopGO.

    close $fileout;
    close $filein;
    print "Now running go chromosome analysis for $species_name\n";

    #Sort the duplicate file:
    `sort $species_name\.go_r_file.noDuplicates.txt | uniq > $species_name\.go_r_file.noDuplicates.uniq.txt`;

    my $out4="$species_name\.go_r_file.noDuplicates.noSameScaffold.txt";
    open(my $fileout4, ">", $out4)   or die "Could not open $out4\n";
    open(my $filein4, "<", "$species_name\.go_r_file.noDuplicates.uniq.txt")   or die "Could not open $species_name\.go_r_file.noDuplicates.uniq.txt\n";
    my %hit;
    my %bad;
    while (my $line4=<$filein4>){
        chomp $line4;
        my @split=split("\t", $line4);
        if ($hit{$split[0]}){
            #delete $hit{$split[0]};
            $bad{$split[0]}="duplicate_scaffold";
        }
        $hit{$split[0]}=$split[1];
    }

    foreach my $key ( keys %hit ){
        if ($bad{$key}){
            #don't use it
        }
        else{
            print $fileout4 "$key\t$hit{$key}\n";
        }
    }

    close $fileout4;
    close $filein4;

    #Remove background duplicated orthogroups.
    my $out5="$species_name\.go_r_file.noDuplicates_BK.noDuplicates.txt";
    open(my $fileout5, ">", $out5)   or die "Could not open $out5\n";
    open(my $filein5, "<", "$species_name\.go_r_file.noDuplicates_BK.txt")   or die "Could not open $species_name\.go_r_file.noDuplicates_BK.txt\n";
    while (my $line5=<$filein5>){
        chomp $line5;
        my @split=split("\t", $line5);
        if ($bad{$split[0]}){
            #do nothing
        }
        else{
            print $fileout5 "$line5\n";
        }
    }
    close $fileout5;
    close $filein5;

    `mkdir Unfiltered_Go_$species_name`;
    `ChopGO_ChromoGoatee.pl -i $species_name\.go_r_file.txt --GO_file $go -sp $species_name`;
    `mv *_res.tab Unfiltered_Go_$species_name`;
    `mv *_res.tab.pdf Unfiltered_Go_$species_name`;
    `mkdir Filtered_dup_Go_$species_name`;
    `ChopGO_ChromoGoatee.pl -i $species_name\.go_r_file.noDuplicates.noSameScaffold.txt --GO_file $species_name\.go_r_file.noDuplicates_BK.noDuplicates.txt -sp $species_name`;
    `mv *_res.tab Filtered_dup_Go_$species_name`;
    `mv *_res.tab.pdf Filtered_dup_Go_$species_name`;
}

print "\n\nScript finished\n";
