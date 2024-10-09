#!/usr/bin/perl
use warnings;
use strict;
use Scalar::Util qw(looks_like_number);

print "Please be in folder with N0.tsv, Base/Gamma_change.tab, Base/Gamma_branch_probabilities.tab  and the Go folder\n";

my $pval = $ARGV[0];
my $type = $ARGV[1];
my $go_max_plot = $ARGV[2];

#Set up output name
my $outname1="CAFE_summary.txt";
open(my $out1, ">", $outname1)   or die "Could not open $outname1\n";

#Read in OG hash
my $go="OG_GO_format.tsv";

#Read in Orthofinder file
my $file=`ls N0.tsv`;
chomp $file;
open(my $filein, "<", $file)   or die "Could not open $file\n";
my $header=<$filein>;
#get rid of mac and and endings
$header =~ s/\n//g;
$header =~ s/\r//g;
chomp $header;
my @head_N0=split("\t", $header);

#print "header : $header\n";

my %HOG_TO_OG;
my %GENE_DATA;
my %GENE_DATA_OG;
my %OG_DATA;
#Make a hash for the whole line per HOG
my %HOG_to_line;

#Create a background OG list for each species, for GO. See in next section.
my %Background_OGs;
print "Create Background OG list\n";

while (my $line = <$filein>){
    chomp $line;
    $line =~ s/\n//g;
    $line =~ s/\r//g;
    my @splitl=split("\t", $line);
    my $hog=$splitl[0];
    my $og=$splitl[1];
    $HOG_TO_OG{$hog}=$og;
    my $n=0;
    #add hog to line:
    $HOG_to_line{$hog}=$line;
    foreach my $col (@splitl){

        #Save -> Species, HOG, and genes
        #This will print any gene that is associated with each OG, which may be in a HOG not changing. 
        $GENE_DATA{$head_N0[$n]}{$hog}=$col;

        if ($GENE_DATA_OG{$head_N0[$n]}{$og}){
            my $old_h=$GENE_DATA_OG{$head_N0[$n]}{$og};
            $GENE_DATA_OG{$head_N0[$n]}{$og}="$old_h\, $col";
            #print "IOt happens  $head_N0[$n]   $og   $old_h  GENEs  $col \n";
        }
        else{
            $GENE_DATA_OG{$head_N0[$n]}{$og}=$col;
        }

        $OG_DATA{$head_N0[$n]}{$hog}=$og;

        #Create a background for each species, to show which OGs are present in each species, and hence were tested in each species (Background for GO).
        if ($col){
            if ($Background_OGs{$head_N0[$n]}){
                my $olde=$Background_OGs{$head_N0[$n]};
                $Background_OGs{$head_N0[$n]}="$olde\n$og";
            }
            else{
                $Background_OGs{$head_N0[$n]}=$og;
            }
        }
        $n++;
        #print "<$col>";
    }
    #my $len=scalar(@splitl);
    #print "$hog $n\n";
}

#Read in Pvalue and count file
my %PVALUE_DATA;
my %COUNT_DATA;
my $file2=`ls Out_cafe/Base_branch_probabilities.tab`;
my $file3=`ls Out_cafe/Base_change.tab`;
chomp $file2;
chomp $file3;
open(my $filein2, "<", $file2)   or die "Could not open $file2\n";
open(my $filein3, "<", $file3)   or die "Could not open $file3\n";
my $header2=<$filein2>;
my $header3=<$filein3>;
chomp $header2;
chomp $header3;

my @head_pval=split("\t", $header2);
my @head_pval_simp;
foreach my $colh (@head_pval){
    if ($colh =~ m/^</g){
        #print "$colh\n";
        my @spo=split("\>", $colh);
        my @spo2=split("\<", $spo[0]);
        #print "$spo2[1]\n";
        push (@head_pval_simp, $spo2[1]);
    }
    else{
        my @spo=split("\<", $colh);
        push (@head_pval_simp, $spo[0]);
    }
}

print "Reading in Out_cafe/Base_branch_probabilities.tab\n";
while (my $line2 = <$filein2>){
    chomp $line2;
    my @splitl=split("\t", $line2);
    my $hog=$splitl[0];
    my $n=0;
    foreach my $col (@splitl){
        #Save -> Species, HOG, and genes
        my $sp=$head_pval_simp[$n];
        $PVALUE_DATA{$sp}{$hog}=$col;
        $n++;
    }
}

print "Reading in Out_cafe/Base_change.tab\n";
while (my $line_count = <$filein3>){
    chomp $line_count;
    my @splitc=split("\t", $line_count);
    my $hog=$splitc[0];
    my $m=0;
    foreach my $col (@splitc){
        my $sp=$head_pval_simp[$m];
        $COUNT_DATA{$sp}{$hog}=$col;
        $m++;
    }
}

my %SPECIES_TOTAL;
my %SPECIES_EXPANSION;
my %SPECIES_CONTRACTION;

my %SPECIES_EXPANSION_SUM;
my %SPECIES_CONTRACTION_SUM;

print "Calculate significant expanded and contracted gene families\n";

foreach my $species (keys %PVALUE_DATA){


    #set up output gene list files 
    my $outnagenes5="$species\.sig.neg.genes.txt";
    open(my $outgene5, ">", $outnagenes5)   or die "Could not open $outnagenes5\n";
    my $outnagenes6="$species\.sig.pos.genes.txt";
    open(my $outgene6, ">", $outnagenes6)   or die "Could not open $outnagenes6\n";

    print $outgene6 "$header\n";
    print $outgene5 "$header\n";

    #print "SP $species\n";
    foreach my $hog (keys %{$PVALUE_DATA{$species}}){
        my $pval=$PVALUE_DATA{$species}{$hog};
        #print "HOGGY $hog\n";
        if (looks_like_number($pval) ){
            if($pval <= 0.05){
                my $direction=$COUNT_DATA{$species}{$hog};
                #print "HERE $direction\n";
                $SPECIES_TOTAL{$species}++;
                if ($direction eq "\+0" || $direction eq "0"){
                    #Do nothing
                    print "HEY, why is this true:\n$species $hog $PVALUE_DATA{$species}{$hog} has a direction $direction\n";
                }
                elsif($direction =~ m/\-/g){
                    $SPECIES_CONTRACTION{$species}++;
                    $SPECIES_CONTRACTION_SUM{$species}+=$direction;

                    #print the hog lines of contractions
                    print $outgene5 "$HOG_to_line{$hog}\n";
                }
                else {
                    $SPECIES_EXPANSION{$species}++;
                    $SPECIES_EXPANSION_SUM{$species}+=$direction;

                    #print the hog lines of expansions
                    print $outgene6 "$HOG_to_line{$hog}\n";
                }
            }
        }
    }
}



#Count the number of HOGs with + or - change, not regarding pvalue.
my %SPECIES_TOTAL_C;
my %SPECIES_EXPANSION_C;
my %SPECIES_CONTRACTION_C;
my %SPECIES_C_OGS;

print "Run through the Count data\n";

#Run through the COUNT DATA and count 
foreach my $species3 (keys %COUNT_DATA){
    #print "HEEER: $species3\n";

    #set up output gene list files 
    my $outnagenes3="$species3\.neg.genes.txt";
    open(my $outgene3, ">", $outnagenes3)   or die "Could not open $outnagenes3\n";
    my $outnagenes4="$species3\.pos.genes.txt";
    open(my $outgene4, ">", $outnagenes4)   or die "Could not open $outnagenes4\n";


    foreach my $hog (keys %{$COUNT_DATA{$species3}}){
        if ($COUNT_DATA{$species3}{$hog} eq "\+0" || $COUNT_DATA{$species3}{$hog} eq "0"){
            #Do nothing
            #print "NO $species3 $hog $COUNT_DATA{$species3}{$hog}\n";
        }
        elsif($COUNT_DATA{$species3}{$hog} =~ m/\-/g){
        
            # Just a print line to check things are working.->   print "OORR $species3 $hog $COUNT_DATA{$species3}{$hog}  $HOG_TO_OG{$hog}\n";
            $SPECIES_CONTRACTION_C{$species3}++;
            $SPECIES_TOTAL_C{$species3}++;
            my $og=$HOG_TO_OG{$hog};



            
            if ($SPECIES_C_OGS{$species3}{'neg'}){
                my $old=$SPECIES_C_OGS{$species3}{'neg'};
                $SPECIES_C_OGS{$species3}{'neg'}="$old\n$og";
                my $genes_related_to_OG=$GENE_DATA_OG{$species3}{$og};

                #print "genes: $genes_related_to_OG $species3 $og\n";
                if ($genes_related_to_OG){
                    print $outgene3 "$og\t$genes_related_to_OG\n";
                }
                else{
                    print $outgene3 "$og\tNo_genes\n";
                }
                
            }
            else{
                $SPECIES_C_OGS{$species3}{'neg'}="$og";

                #print "genes: $genes_related_to_OG $species3 $og\n";
                my $genes_related_to_OG=$GENE_DATA_OG{$species3}{$og};
                if ($genes_related_to_OG){
                    print $outgene3 "$og\t$genes_related_to_OG\n";
                }
                else{
                    print $outgene3 "$og\tNo_genes\n";
                }
            }
        }
        else {
            #print "X   $species3 $hog $COUNT_DATA{$species3}{$hog}\n";
            $SPECIES_EXPANSION_C{$species3}++;
            $SPECIES_TOTAL_C{$species3}++;
            my $og=$HOG_TO_OG{$hog};
            
            if ($SPECIES_C_OGS{$species3}{'pos'}){
                my $old=$SPECIES_C_OGS{$species3}{'pos'};
                $SPECIES_C_OGS{$species3}{'pos'}="$old\n$og";

                #print "genes: $genes_related_to_OG $species3 $og\n";
                my $genes_related_to_OG=$GENE_DATA_OG{$species3}{$og};
                if ($genes_related_to_OG){
                    print $outgene4 "$og\t$genes_related_to_OG\n";
                }
                else{
                    print $outgene4 "$og\tNo_genes\n";
                }
            }
            else{
                $SPECIES_C_OGS{$species3}{'pos'}="$og";

                #print "genes: $genes_related_to_OG $species3 $og\n";
                my $genes_related_to_OG=$GENE_DATA_OG{$species3}{$og};
                if ($genes_related_to_OG){
                    print $outgene4 "$og\t$genes_related_to_OG\n";
                }
                else{
                    print $outgene4 "$og\tNo_genes\n";
                }
            }
        }
    }
}

print "Print GO_lists\n";

foreach my $species4 (keys %SPECIES_C_OGS){
    #print "$species4 $SPECIES_C_OGS{$species4}\n";
    if ($species4 eq "#FamilyID" || $species4 eq "FamilyID"){
        #Do nothing, this is not a useful line
    }
    else{

        if ($SPECIES_C_OGS{$species4}{'pos'}){
            my $outna2="$species4\.pos.txt";
            open(my $out2, ">", $outna2)   or die "Could not open $outna2\n";
            print $out2 "$SPECIES_C_OGS{$species4}{'pos'}\n";
            close $out2;
        }
        
        if ($SPECIES_C_OGS{$species4}{'neg'}){
            my $outna3="$species4\.neg.txt";
            open(my $out3, ">", $outna3)   or die "Could not open $outna3\n";
            print $out3 "$SPECIES_C_OGS{$species4}{'neg'}\n";
            close $out3; 
        }

    }

}

print "Print Background GO lists\n";

foreach my $species6 (keys %Background_OGs){
    chomp $species6;

    if ($Background_OGs{$species6}){
        if ($species6 =~ m/Gene Tree Parent Clade/){
            #ignore it, this is the header.
        }
        else{
            #print "ITS HERE: $species6 \n";
            my $out_back="$species6\.BK.txt";
            open(my $outb, ">", $out_back)   or die "Could not open $out_back\n";
            print $outb "$Background_OGs{$species6}\n";
            `sort $species6\.BK.txt | uniq > $species6\.BK.txt.uniq`;
        }

    }
    
}




print "Print a summary table:\n";

print $out1 "Total_HOGs_significant\tExpansion_HOGs_significant\tContraction_HOGs_significant\tExpansion_genes_significant\tContraction_genes_significant\tExpansion_HOGs_total\tContraction_HOGs_total\n";

foreach my $species2 (sort keys %SPECIES_TOTAL){
    print $out1 "$species2\t$SPECIES_TOTAL{$species2}";
    if ($SPECIES_EXPANSION{$species2}){
        print $out1 "\t$SPECIES_EXPANSION{$species2}";
    }
    else{
        print $out1 "\t0";
    }
    if ($SPECIES_CONTRACTION{$species2}){
        print $out1 "\t$SPECIES_CONTRACTION{$species2}";
    }
    else{
        print $out1 "\t0";
    }


    if ($SPECIES_EXPANSION_SUM{$species2}){
        print $out1 "\t$SPECIES_EXPANSION_SUM{$species2}";
    }
    else{
        print $out1 "\t0";
    }
    if ($SPECIES_CONTRACTION_SUM{$species2}){
        print $out1 "\t$SPECIES_CONTRACTION_SUM{$species2}";
    }
    else{
        print $out1 "\t0";
    }

    #Not significant expansions and contractions:
    if ($SPECIES_EXPANSION_C{$species2}){
        print $out1 "\t$SPECIES_EXPANSION_C{$species2}";
    }
    else{
        print $out1 "\t0";
    }
    if ($SPECIES_CONTRACTION_C{$species2}){
        print $out1 "\t$SPECIES_CONTRACTION_C{$species2}";
    }
    else{
        print $out1 "\t0";
    }
    print $out1 "\n";
}



#Now run Chopgo

foreach my $species5 (keys %SPECIES_TOTAL){
    print "Running GO on $species5\n";
    if (looks_like_number($species5)){
    	#`mv $species5\.pos.txt Node_$species5\.pos.txt`;
    	#`mv $species5\.neg.txt Node_$species5\.neg.txt`;

        if ( -e "$species5\.pos.txt" ){
            `cut -f 1 $species5\.pos.txt > Node_$species5\.pos.txt`;
            my $line_count = `wc -l Node_$species5\.pos.txt | awk '{print \$1}'`;
            if ( $line_count >= 10){ 
                print "ChopGO_VTS2.pl -i Node_$species5\.pos.txt --GO_file $go -bg OG_GO_format.tsv -pval $pval -pval_type $type -max_plot $go_max_plot\n";
                `ChopGO_VTS2.pl -i Node_$species5\.pos.txt --GO_file $go -bg OG_GO_format.tsv -pval $pval -pval_type $type -max_plot $go_max_plot`;
            }
            if ( $line_count <= 100){ 
                print "WARNING, GO likely to fail on $line_count genes (Run $species5)\n";
            }

        }
    	if ( -e "$species5\.neg.txt" ){
            `cut -f 1 $species5\.neg.txt > Node_$species5\.neg.txt`;
            my $line_count = `wc -l Node_$species5\.neg.txt | awk '{print \$1}'`;
            if ( $line_count >= 10){ 
                print "ChopGO_VTS2.pl -i Node_$species5\.neg.txt --GO_file $go -bg OG_GO_format.tsv -pval $pval -pval_type $type -max_plot $go_max_plot\n";
                `ChopGO_VTS2.pl -i Node_$species5\.neg.txt --GO_file $go -bg OG_GO_format.tsv -pval $pval -pval_type $type -max_plot $go_max_plot`;
            }
            if ( $line_count <= 100){ 
                print "WARNING, GO likely to fail on $line_count genes (Run $species5)\n";
            }
        }

    }
    else{
        if ( -e "$species5\.pos.txt" ){
            my $line_count = `wc -l $species5\.pos.txt | awk '{print \$1}'`;
            if ( $line_count <= 100){ 
                print "WARNING, GO likely to fail on $line_count genes (Run $species5)\n";
            }
            if ( $line_count >= 10){ 
                `ChopGO_VTS2.pl -i      $species5\.pos.txt --GO_file $go -bg $species5\.BK.txt.uniq  -pval $pval -pval_type $type -max_plot $go_max_plot`;
            }
        }
        if( -e "$species5\.neg.txt" ){
            my $line_count = `wc -l $species5\.neg.txt | awk '{print \$1}'`;
            if ( $line_count <= 100){ 
                print "WARNING, GO likely to fail on $line_count genes (Run $species5)\n";
            }
            if ( $line_count >= 10){ 
            `ChopGO_VTS2.pl -i      $species5\.neg.txt --GO_file $go -bg $species5\.BK.txt.uniq  -pval $pval -pval_type $type -max_plot $go_max_plot`;
            }
        }
    }
}


close $out1;

print "\nFinished Running\n\n";



