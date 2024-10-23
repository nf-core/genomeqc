#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use Cwd qw(abs_path cwd);
use File::Path qw(make_path);

# INITIALIZE PATH
my $binPath = abs_path($0);
$0 =~ s/^.*\///;
$binPath =~ s/\/$0$//; # vast-tools/bin

my $GO_file;
my $input;
my $background;
my $species;
my $helpFlag=0;
my $path_to_DB="$binPath/../VASTDB";
my $plot = 1;
my $pval_cutoff=0.05;
my $num_cutoff=10;
my $plot_only=0;
my $plot_open=0;
my $install_topGO;
my $meth="none";
my $sort_enrich=0;
my $min_genes=2;
my $max_genes=100000;
my $infer=0;

GetOptions(	  "help" => \$helpFlag,
		  'i=s' => \$input,
		  "bg=s" => \$background,
		  'sp=s' => \$species,
		  "db=s" => \$path_to_DB,
		  "GO_file=s" => \$GO_file,
		  "plot" => \$plot,
		  "pval=s" => \$pval_cutoff,
		  "pval_type=s" => \$meth,
		  "max_plot=s" => \$num_cutoff,
		  "plot_only" => \$plot_only,
		  "open_plot" => \$plot_open,
		  "install_topGO" => \$install_topGO,
		  "filt_enrich=s" => \$sort_enrich,
		  "min_genes_req=s" => \$min_genes,
		  "max_genes_req=s" => \$max_genes,
		  "allow_inferred" => \$infer,
    );



my $EXIT_STATUS=0;
if (defined $input && ((defined $species && defined $path_to_DB) || (defined $GO_file))){
    print "Starting...\n\n";
}
else{
    $EXIT_STATUS=1;
}

if ($helpFlag or $EXIT_STATUS){
    die "
usage: ChopGO.pl -i <INPUT_list> (-sp <SPECIES_CODE> -db </path/to/DB> OR --GO_file Custom_GO_file) [options]

compulsory:
    -i               Input list of genes (compulsory). Sep by \'\\n\'. Optional second column for subdivisions of list.  

    -sp              Species must be (Hsa,Bla,Mmu code). Check availability in DB.
OR
    --GO_file        Custom file with GeneID\tGO_term associations

optional:
    -bg              Background must be a list of genes sep by \\n.
    -db              Database, give path to DB of GO annotations. Default (./VASTDB)
    -install_topGO   Install topGO in R.

plotting options:

    -plot            Plot (BOOLEAN, optional), will create a histogram/s. Default=0=(OFF).
    -pval            Choose pval cutoff for histogram plots. Default=0.05.
    -pval_type	     Choose correction (holm, hochberg, hommel, bonferroni, BH, BY, fdr or none). Default=none.
    -filt_enrich     Filter by fold enrichment. Default=0.
    -max_plot        Max number of results to plot for each histogram (Bp,MF,CC). Default=10.
    -min_genes_req   Choose minimum number of genes containing each GO term. Default=2.
    -max_genes_req   Choose maximum number of genes containing each GO term. Default=100000.
    -allow_inferred  Allow inferred parent (GO, not in original file) into final enrichment. Default=0=OFF.
    -plot_only       Plot Only (BOOLEAN, optional), don't do GO enrichments. Default=0=OFF.
    -open_plot       Open plots once made (BOOLEAN, optional). Default=0=OFF.

Pre-requisites:
    R 

*** Questions: 
    Chris Wyatt (chris.wyatt\@crg.eu)
    Manuel Irimia (mirimia\@gmail.com)

"
}


my $path_to_GOs;
#my inferred file path
if ($infer){}
elsif ($GO_file){
    $path_to_GOs=$GO_file;
    $infer=$path_to_GOs;
}
else {
    $path_to_GOs="$path_to_DB\/$species\/FILES/GO_FILE_$species";
    $infer=$path_to_GOs;
}

if (-e $path_to_GOs){
	#good for you
}
else{
	print "

GO file does not exist:\n
= $path_to_GOs\n
This may be because there is no GO file in the DB you are pointing at.\n
If this is the case, you might do the following: \n
1,  Go to https://github.com/vastgroup/vast-tools.
2,  Download the GO_file of the species of interest.

Alternatively, you may download your own annotation, e.g.:

1,  Go to http://www.ensembl.org/index.html.\n
2,  Go to BioMart.\n
3,  Choose database Ensembl Genes (latest version).\n
4,  Choose species.\n 
5,  Click Attributes.\n
6,  Under Gene, select \"Gene ID\" only.\n
7,  Under External, select : \"GO Term Accession\" and \"GO Term Name\".\n
8,  Click results, Then Go (to save the file to your computer).\n
9,  Finally, make sure file is tab separated, and it goes GENE ID (e.g. ENG000..), GO (e.g.GO:000001), GO name (e.g. Liver related). IN this order.\n
10, Save file as GO_FILE_species (e.g. GO_FILE_Hsa), in the DB folder under the same three letter species name.\n
	"
}


#if you just want to plot
if ($plot_only == 0){

#Files to be printed, filled in later.
my @ALL_made_files;

my @methods=("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none");



### DATABASE CHOICE, print to user
print "GO = $species GO\nDb = $path_to_DB\n\n" if !defined $GO_file;
print "GO file = $GO_file\n\n" if defined $GO_file;

###BACKGROUND###
my $outfile="BACKGROUND\.forR";
my %Background_hash;
if($background){
    ## FIND BACKGROUND GENES
    open(my $IN_b, "<", $background)   or die "Could not open $background \n";
    while (my $line=<$IN_b>){
	$line=~s/\r//g;
	chomp $line;
	my @input_here=split("\t", $line);
	if (scalar @input_here > 1.5){
	    $Background_hash{$input_here[0]}="HIT";
	}
	else{
	    $Background_hash{$line}="HIT";
	}
	
    }
    my $n_back = keys %Background_hash;
    
    ###MAKE BACKGROUND FOR R;
    my $Table_b1;
    if ($GO_file){
	$Table_b1 = $GO_file;
    }
    else {
	$Table_b1 = "$path_to_DB\/$species\/FILES/GO_FILE_$species";
    }
    open(my $IN, "<", $Table_b1)   or die "Could not open $Table_b1 \n";
    open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";
    my %Gene_Go_Hash;
    my %tot_genes;
    while (my $line=<$IN>){
	$line=~s/\r//g;
	chomp $line;
	my @linesplit= split("\t", $line);
	my $gene=$linesplit[0];
	$tot_genes{$gene}="HIT";
	if ($Background_hash{$gene}){
	    my $GO=$linesplit[1];
	    if ($GO){
		if ($Gene_Go_Hash{$gene}){
		    my $old=$Gene_Go_Hash{$gene};
		    #print  "2+ $gene  $old\t$GO\n";
		    $Gene_Go_Hash{$gene}="$old\",\"$GO";
		}
		else{
		    $Gene_Go_Hash{$gene}=$GO;
		    #print "1st $gene $GO\n";
		}
	    }
	}
    }
    my $n_all=keys %tot_genes;
    print "BACKGROUND\nOf $n_all total Genes with GO annotation, $n_back were chosen for the background (based on user -b list)\n\n";
    
    print $outhandle "Chop.gene2GO<- list()\n";
    foreach my $key ( keys %Gene_Go_Hash ){
	print $outhandle "Chop.gene2GO\$$key <- c(\"$Gene_Go_Hash{$key}\")\n";
    }
}
else {
    #all genes as background. No background chosen by user.
    ###MAKE BACKGROUND FOR R;
    my $Table_b2;
    if ($GO_file){
	$Table_b2 = $GO_file;
    }
    else {
        $Table_b2 = "$path_to_DB\/$species\/FILES/GO_FILE_$species";
    }
    open(my $IN, "<", $Table_b2)   or die "Could not open $Table_b2 \n";
    open(my $outhandle, ">", $outfile)   or die "Could not open $outfile \n";
    my %Gene_Go_Hash;
    while (my $line=<$IN>){
	$line=~s/\r//g;
	chomp $line;
	#$line=~s/-/_/g;
	#$line=~s/:/_/g;
	my @linesplit= split("\t", $line);
	my $gene=$linesplit[0];
	$gene=~s/-/_/g;
	$gene=~s/:/_/g;
	my $GO=$linesplit[1];
	if ($GO){
	    if ($Gene_Go_Hash{$gene}){
		my $old=$Gene_Go_Hash{$gene};
		#print  "2+ $gene  $old\t$GO\n";
		$Gene_Go_Hash{$gene}="$old\",\"$GO";
	    }
	    else{
		$Gene_Go_Hash{$gene}=$GO;
		#print "1st $gene $GO\n";
	    }
	}
    }
    
    print $outhandle "Chop.gene2GO<- list()\n";
    foreach my $key ( keys %Gene_Go_Hash ){
	print $outhandle "Chop.gene2GO\$$key <- c(\"$Gene_Go_Hash{$key}\")\n";
    }
}


###QUERY INPUT###

my $Table_i = $input;
open(my $IN, "<", $Table_i)   or die "Could not open $Table_i \n";

if ($Table_i=~ m/\//){
    my @split=split ("\/", $Table_i);
    $Table_i=$split[-1];
    $Table_i=~ s/\-/\_/g;
    
}

if ($Table_i =~ m/\-/){
    print "Name of input file contains \"-\"\'s , replacing with \"\_\"\'s, for processing in R\n";
    $Table_i=~ s/\-/\_/g;
}

my $outfile2="$Table_i\.converted";
my $outfile3="$Table_i\.R_GO_subs";

open(my $outhandle2, ">", $outfile2)   or die "Could not open $outfile2 \n";
open(my $outhandle3, ">", $outfile3)   or die "Could not open $outfile3 \n";

print "\nStep 2\n\n";

my %chromo_score;

my %Gene_Go_Hash;
while (my $line=<$IN>){
    $line =~ s/\r//g;
    $line =~ s/\"//g;
    chomp $line;
    #print "$line\n";
    my @linesplit= split("\t", $line);
    my $gene=$linesplit[0];
    $gene=~s/-/_/g;
    $gene=~s/:/_/g;


    if ($linesplit[1]){
		my $GROUP=$linesplit[1];
		if ($Gene_Go_Hash{$GROUP}){
		    my $old=$Gene_Go_Hash{$GROUP};
		    $Gene_Go_Hash{$GROUP}="$old\",\"$gene";
		    $chromo_score{$GROUP}++;
		}
		else{
		    $Gene_Go_Hash{$GROUP}=$gene;
		    $chromo_score{$GROUP}=1;
		}
    }
    else{
		my $GROUP="$Table_i";
		if ($Gene_Go_Hash{$GROUP}){
		    my $old=$Gene_Go_Hash{$GROUP};
		    $Gene_Go_Hash{$GROUP}="$old\",\"$gene";
		}
		else{
		    $Gene_Go_Hash{$GROUP}=$gene;
		}
    }
}	


#Go through chromosomes and if any contain 20 or less genes, ignore them, as topGO may fail.
if (%chromo_score){
	foreach my $key ( keys %chromo_score ){
		my $val=$chromo_score{$key};
		if ($val <=20){
			#Remove entry from hash, as its too small really for GO.
			delete($Gene_Go_Hash{$key});
		}
	}
}

#Read through multi GROUP file and ensure no lists are less than 20 genes (which won't often run through GO).



### RUN R COMPARISONS ###

print $outhandle2 "Chop.WGCNA2Gene<- list()\n";
foreach my $key ( keys %Gene_Go_Hash ){
    
    print $outhandle2 "Chop.WGCNA2Gene\$$key <- c(\"$Gene_Go_Hash{$key}\")\n";
    
    #print GO runs
    print $outhandle3 "selGenes<-Chop.WGCNA2Gene\$",$key,"\n";
    print $outhandle3 "inGenes <- factor(as.integer(names(Chop.gene2GO) %in% selGenes))\n";
    print $outhandle3 "names(inGenes) <- names(Chop.gene2GO)\n";
    
    #BP
    print $outhandle3 "GOdata <- new(\"topGOdata\", ontology=\"BP\", allGenes=inGenes, annot=annFUN.gene2GO, gene2GO=Chop.gene2GO)\n";
    print $outhandle3 "resultFisher <- runTest(GOdata, algorithm = \"classic\", statistic = \"fisher\")\n";
    print $outhandle3 "allRes_BP <- GenTable(GOdata, classicFisher = resultFisher,orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 50)\n";
    foreach my $meths (@methods){
	print $outhandle3 "allRes_BP\$",$meths,"<-p.adjust(allRes_BP\$classicFisher, method = \"",$meths,"\")\n";
    }
    print $outhandle3 "allRes_BP\$FoldChange<-allRes_BP\$Significant/allRes_BP\$Expected\n";
    print $outhandle3 "allRes_BP\$ontology<-\"BP\"\n";
    
    #MF
    print $outhandle3 "GOdata <- new(\"topGOdata\", ontology=\"MF\", allGenes=inGenes, annot=annFUN.gene2GO, gene2GO=Chop.gene2GO)\n";
    print $outhandle3 "resultFisher <- runTest(GOdata, algorithm = \"classic\", statistic = \"fisher\")\n";
    print $outhandle3 "allRes_MF <- GenTable(GOdata, classicFisher = resultFisher,orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 50)\n";
    foreach my $meths (@methods){
	print $outhandle3 "allRes_MF\$",$meths,"<-p.adjust(allRes_MF\$classicFisher, method = \"",$meths,"\")\n";
    }
    print $outhandle3 "allRes_MF\$FoldChange<-allRes_MF\$Significant/allRes_MF\$Expected\n";
    print $outhandle3 "allRes_MF\$ontology<-\"MF\"\n";
    
    #CC
    print $outhandle3 "GOdata <- new(\"topGOdata\", ontology=\"CC\", allGenes=inGenes, annot=annFUN.gene2GO, gene2GO=Chop.gene2GO)\n";
    print $outhandle3 "resultFisher <- runTest(GOdata, algorithm = \"classic\", statistic = \"fisher\")\n";
    print $outhandle3 "allRes_CC <- GenTable(GOdata, classicFisher = resultFisher,orderBy = \"classicFisher\", ranksOf = \"classicFisher\", topNodes = 50)\n";
    foreach my $meths (@methods){
	print $outhandle3 "allRes_CC\$",$meths,"<-p.adjust(allRes_CC\$classicFisher, method = \"",$meths,"\")\n";
    }
    print $outhandle3 "allRes_CC\$FoldChange<-allRes_CC\$Significant/allRes_CC\$Expected\n";
    print $outhandle3 "allRes_CC\$ontology<-\"CC\"\n";
    
    #ALL MERGE
    print $outhandle3 "ALL_res_",$key,"<-rbind(allRes_BP, allRes_MF, allRes_CC)\n";
    print $outhandle3 "x.sub <- subset(ALL_res_",$key,", none < ",$pval_cutoff,")\n";
    print $outhandle3 "f.sub<-x.sub[sort.list(x.sub\$none),]\n";
    print $outhandle3 "e.sub <- subset(f.sub, FoldChange > ",$sort_enrich,")\n";
    print $outhandle3 "write.table(e.sub, \"$species\_$key\_res.tab\", sep=\"\\t\", quote=FALSE, eol=\"\\n\", row.names=F)   \n";
    my $out_name="$species\_$key\_res.tab";
    push (@ALL_made_files, $out_name);
}

print "Finished compiling scripts, now running R... (can take ~40 seconds for each query)\n\n";

### RUN MAIN CODE ###
my $outfile4="R_CODE_TO_RUN";
open(my $outhandle4, ">", $outfile4)   or die "Could not open $outfile4 \n";
if (defined $install_topGO){
    print $outhandle4 "source(\"http://bioconductor.org/biocLite.R\")\n";
    print $outhandle4 "biocLite(\"topGO\")\n";
}

print $outhandle4 "suppressWarnings(suppressMessages(library (topGO)))\n";
print $outhandle4 "source (\"$outfile\")\nsource (\"$outfile2\")\nsource (\"$outfile3\")\n";
print $outhandle4 "save.image(\"Image.Rdata\")\nquit(save = \"no\")\n";

`R --vanilla < R_CODE_TO_RUN >output.ofthis.test`;

print "Start Plotting\n\n";

### RUN HITOGRAM CODE ###

if (defined $plot){
    my @files=@ALL_made_files;
    foreach my $results (@files){
	chomp $results;
	my $in_here="$binPath\/MakeHist_fromChartGO-jcvi.pl";
	if ($in_here=~ m/ /){
	    $in_here=~ s/ /\\ /g;
	}
	if ($in_here=~ m/\(/){
	    $in_here=~ s/\(/\\(/g;
	}
	if ($in_here=~ m/\)/){
	    $in_here=~ s/\)/\\)/g;
	}
	#print "HERE: $in_here\n";
	`$in_here $results $pval_cutoff $num_cutoff $plot_open $meth $sort_enrich $min_genes $max_genes $infer`;
    }
}

#if you just want to plot. Miss all the other steps.
}
else {
    my $Table_i = $input;
    open(my $IN, "<", $Table_i)   or die "Could not open $Table_i \n";
    my %Gene_Go_Hash;
    while (my $line=<$IN>){
		$line =~ s/\r//g;
		$line =~ s/\"//g;
		chomp $line;
		#print "$line\n";
		my @linesplit= split("\t", $line);
		my $gene=$linesplit[0];
		if ($linesplit[1]){
		    my $GROUP=$linesplit[1];
		    if ($Gene_Go_Hash{$GROUP}){
				my $old=$Gene_Go_Hash{$GROUP};
				$Gene_Go_Hash{$GROUP}="$old\",\"$gene";
		    }
		    else{
				$Gene_Go_Hash{$GROUP}=$gene;
		    }
		}
		else{
		    my $GROUP="$Table_i";
		    if ($Gene_Go_Hash{$GROUP}){
				my $old=$Gene_Go_Hash{$GROUP};
				$Gene_Go_Hash{$GROUP}="$old\",\"$gene";
		    }
		    else{
				$Gene_Go_Hash{$GROUP}=$gene;
		    }
		}
    }	
    
    my @ALL_made_files;
    foreach my $key ( keys %Gene_Go_Hash ){
		my $out_name="$species\_$key\_res.tab.tab";
		push (@ALL_made_files, $out_name);
    }
    
    
    foreach my $results (@ALL_made_files){
	chomp $results;
	my $in_here="$binPath\/MakeHist_fromChartGO-jcvi.pl";
	if ($in_here=~ m/ /){
	    $in_here=~ s/ /\\ /g;
	}
	if ($in_here=~ m/\(/){
	    $in_here=~ s/\(/\\(/g;
	}
	if ($in_here=~ m/\)/){
	    $in_here=~ s/\)/\\)/g;
	}
	#print "HERE: $in_here\n";
	print "$in_here $results $pval_cutoff $num_cutoff $plot_open $meth $sort_enrich $min_genes $max_genes $infer\n";
	`perl $in_here $results $pval_cutoff $num_cutoff $plot_open $meth $sort_enrich $min_genes $max_genes $infer`;
    }
}

print  "Script completed\n\n";

#TIDY UP , temporary files
$input=~s/\-/\_/g;
#`rm -f ACTUAL_R_CODE BACKGROUND.forR R_CODE_TO_RUN output.ofthis.test $input\.converted $input\.R_GO_subs tmp_file Image.Rdata`;


