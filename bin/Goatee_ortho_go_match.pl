#!/usr/bin/perl
#loopfold.sortmerna.pl
use warnings;
use strict;

die "Please specify (1) Orthofinder.csv (2) File focal species fasta\n" unless(@ARGV==2);

my $ORTHOFINDER = $ARGV[0];
open(my $filein, "<", $ORTHOFINDER)   or die "Could not open $ORTHOFINDER\n";

my $FOCAL = $ARGV[1];
#open(my $focal, "<", $FOCAL)   or die "Could not open $FOCAL\n";
if ($FOCAL =~ m/gz$/){
  my @sl=split(".gz", $FOCAL);
  my $NEWFOCAL = "$sl[0]";
  `zcat $FOCAL > $NEWFOCAL`;
  $FOCAL = $NEWFOCAL;
  print "ungzipping\n";
}
else{
  print "nah\n";
}

my $outfile="Result_All_Data";
open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";
my $outfile2="$FOCAL\_Result_All_Combine_GO_format";
open(my $outhandle2, ">", $outfile2)   or die "Could not open $outfile2 \n";
my $outfile3="OG_GO_format.tsv";
open(my $outhandle3, ">", $outfile3)   or die "Could not open $outfile3 \n";

#
print "\n#####################################################\n";
print "Starting script to join orthofinder results to GO IDs\n\n";
print "Step 1: Work out header column order and make species to GO hashes\n\n";


#Read the header off the orthofinder file to match up the species with each column (then later, read the GO files and make a hash)
my $header=<$filein>;
chomp $header;
my @split_head=split("\t", $header);
my $orthogroup_name=shift @split_head; #Remove orthogroup column
print "Header of Orthofinder file:\n$header\n\n";
my $colposition=0;
my %col_to_sp_store;
my %sp_store_to_col;

my %Nested_species_go_store;

foreach my $col (@split_head){
  chomp $col;
  $col =~ s/\r\n|\n|\R|\r//g;
  #my @split_dot=split("\.", $col);
  my $exp_gofile="$col\.go.txt";
  print "Col of Orthofinder table: $col\n";

  print "TESTED: $exp_gofile\n";
  if (-e "$exp_gofile"){
    print "3. Match to GO database\n";
    $col_to_sp_store{$colposition}=$exp_gofile;
    $sp_store_to_col{$exp_gofile}=$colposition;
    #Read in the go txt files and make them into a hash of species gene ID to a list of all GO terms associated (e.g. ENS23523525 -> "GO:000032342,GO:0000248214")
    my $GO="$exp_gofile";
    open(my $goin, "<", $GO)   or die "Could not open $GO\n";
    my $ignoreheader=<$goin>;
    while (my $line=<$goin>){
      chomp $line;
      $line =~ s/\r\n|\n|\R|\r//g;
      my @sp=split("\t", $line);
      my $gene=$sp[1];
      my $go_id=$sp[2];
      #print "$line\n";
      if ($go_id){
        if ($go_id =~ m/GO/){
          
          #print "<-$go_id\->\n";
          if ($Nested_species_go_store{$exp_gofile}{$gene}){
            my $old=$Nested_species_go_store{$exp_gofile}{$gene};
            $Nested_species_go_store{$exp_gofile}{$gene}="$old\t$go_id";
          }
          else{
            $Nested_species_go_store{$exp_gofile}{$gene}=$go_id;
          }
        }
      }
    }
  }
  #If we find the focal column:
  elsif($FOCAL eq "$col\.fa" || $FOCAL eq "$col\.fasta"){
    print "1. This happened $FOCAL eq $col\.fa || $FOCAL eq $col\.fasta \n";
    $col_to_sp_store{$colposition}=$FOCAL;
  }
  #Else, means we probably messed up or a typo:
  else{
    print "2. Did not find a match for $col, expected to be: $exp_gofile\nMake sure this is the behaviour you expect, else, maybe you didn't put the correctly titled file in the folder with the GO or protein files.\n";
  }
  $colposition++;
}




print "\n\nRead in the Orthofinder to GO information ( step 1 )\n\n";
print "Step 2: Matching focal genes to Orthofinder orthologs and stripping GO assiugnments\n\n";

#Now read through the rest of the orthofinder file to match each Orthogroups genes for each species to its corresponding lines in their GO files. So we can get for say orthogroup 000001: Apis has two genes, and they have 30 GO annotations from their ensembl ortho data.
my %focal_to_orthofinder_store;
my %target_orthoid_to_GOterms;
while (my $line=<$filein>){
  chomp $filein;
  #print "$line\n";
  my @split=split("\t", $line);
  my $orthogroup_name=shift @split; #Remove orthogroup column
  my $current_col=0;
  
  foreach my $sp_line (@split){
    #print "$sp_line\n";
    my $species=$col_to_sp_store{$current_col};
    if ($species){
      
      my @genes=split(", ",$sp_line);


      #If the focal species, just add a key from focal gene ID to Orthogroup ID:
      #print "HERE: $FOCAL eq $species\n";
      if ($FOCAL eq "$species"){
        #print "yes\n";
        foreach my $gene (@genes){
          $focal_to_orthofinder_store{$gene}=$orthogroup_name;
        }
      }
      #else for the rest of the species we find the gene go hash entries for each gene:
      else{
        #print "No: $FOCAL eq $species\n";
        foreach my $gene (@genes){
          if ($Nested_species_go_store{$species}{$gene}){
            if ($target_orthoid_to_GOterms{$orthogroup_name}{$species}){
              my $old=$target_orthoid_to_GOterms{$orthogroup_name}{$species};
              $target_orthoid_to_GOterms{$orthogroup_name}{$species}="$old\t$Nested_species_go_store{$species}{$gene}";
            }
            else{
              $target_orthoid_to_GOterms{$orthogroup_name}{$species}=$Nested_species_go_store{$species}{$gene};
            }
          }
        }
      }
    }
    $current_col++;
  }
}


print "\n\nLinked orthofinder IDs to GO terms ( step 2 )\n\nStep 3:Now get the focal IDs to GO terms per species\n\n";

#Then pull everything together.

my %done; #This is to make sure we dont get duplicates in the result;

foreach my $focal_genes ( keys %focal_to_orthofinder_store ){

  my $OG=$focal_to_orthofinder_store{$focal_genes};
  
  #print "$OG\t$focal_genes\t$focal_to_orthofinder_store{$focal_genes}\n";

  
  foreach my $species ( keys %{ $target_orthoid_to_GOterms{$OG} } ){
    #print "here something $species\n";
    my $goterms=$target_orthoid_to_GOterms{$OG}{$species};
    print $outhandle "$focal_genes\t$OG\t$species\t$goterms\n";
    my @goterms_sep=split("\t", $goterms);
    foreach my $terms (@goterms_sep){
      my $combtest="$focal_genes\t$terms";
      if ($done{$combtest}){
        #do nothing, already in final doc.
      }
      else{
        if($focal_genes =~ /^\d/){
          $focal_genes="GENE_$focal_genes";
        }
        print $outhandle2 "$focal_genes\t$terms\n";
        print $outhandle3 "$OG\t$terms\n";
        $done{$combtest}="DONE";
      }
    }
  }
}

close $filein;
close $outhandle;
close $outhandle2;
close $outhandle3;
print "Finished\n";

