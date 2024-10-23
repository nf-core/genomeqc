#!/usr/bin/perl

die "\nUsage: MakeHist_fromChartGO_v2.pl topGO_Chart_output [min_p-value max_number_categories]\n\n".
    "    >> Min Raw p-value for a GO category to be plotted (default 0.01)\n".
    "    >> Maximum number of GO categories to be plotted per Type (BP, MF, CC) (default 10)\n".
    "    >> Open or not (default 0)\n".
    "\n" if !$ARGV[0];

$min_p=$ARGV[1];
$num_results=$ARGV[2];
$open_me=$ARGV[3];
$correction=$ARGV[4]; #none or bonferroni.
$enrich=$ARGV[5]; #0
$min_gene=$ARGV[6]; #6
$max_gene=$ARGV[7]; #100000
$inferred=$ARGV[8]; #FILES/GO_FILE

### Set colors:
$col_BP = "red";
$col_MF = "forestgreen";
$col_CC = "blue3";

### open input and output files
if (-e $ARGV[0]){
    open (I, $ARGV[0]) || die "Can't open the input file \"$ARGV[0]\"\n";
    $root=$ARGV[0];
    $out_file="$root\_$correction\_pval\-$correction\-$min_p.ALL.tab";
    open (O, ">$out_file") || die "Can't open the output file";

    print O "TYPE\tGO\tTERM\tLogP-$correction\n";


    open (I_2, $inferred);
    my %hits;
    while ($line=<I_2>){
        #print "Stuff $line\n";
        chomp $line;
        @t_2=split(/\t/, $line);
        ($go_2)=$t_2[1];
        $hits{$go_2}="HITS";
    }


    $header=<I>;
    chomp $header;
    @split_head=split("\t", $header);
    my $n_h=0;
    my $pval_choice_col;
    foreach my $cols_head (@split_head){
        #print "head: $cols_head $n_h\n";
        if ($cols_head eq $correction){
            $pval_choice_col=$n_h;
            #print "$cols_head eq $correction = $pval_choice_col\n";
        }
        $n_h++;
    }
    print "type $pval_choice_col\n";

    while (<I>){
        #print "Stuff\n";
        if ($_ =~m/GO:0016645/){print "Here\n"};
        my $save=$_;
        chomp;
        #if ($save =~m/GO:0016645/){print "$save $save $save\n"};
        @t=split(/\t/);
        ($go)=$t[0];
        print "dfasfef $go" if ($save =~m/GO:0016645/);
        if ($inferred){
            print "dfasfef $go" if ($save =~m/GO:0016645/);
            if ($hits{$go}){
                #if ($save =~m/GO:0016645/){print "Here2\n"};
                ($term)=$t[1];
                if ($term=~ m/\'/){
                    #print "YES $term\n";
                    $term=~ s/\'/p/g;
                }
                ($type)=$t[-1];
                ($fold_enrich)=$t[-2];
                if ($fold_enrich=~ m/inf/){
                    $fold_enrich="100";
                }
                ($annotated)=$t[2];

                if ($annotated <= $min_gene){
                    #Do nothing
                    #if ($save =~m/GO:0016645/){print "Here2\n"};
                }
                else{
                    #if ($save =~m/GO:0016645/){print "Here2\n"};
                    if ($annotated >= $max_gene){
                       #Do nothing
                    }
                    else{
                        $p;
                        $e;
                        $log_p;
                        if ($fold_enrich >= $enrich){
                            #$e=$t[-2];
                            $p=$t[$pval_choice_col];
                            #print "$p \n";
                            if ($p ==0){
                                $log_p=$p;
                            }
                            elsif ($p =~ m /\< 1e-30/){
                                $log_p="1e-30";
                            }
                            else{
                                $log_p=sprintf("%.2f",-log($p)/log(10));
                            }
                            
                            if ($p<$min_p){
                                $terms{$go}=$term;
                                push(@{$order{$type}},"$log_p=$go");
                            }
                            else{
                                #print "NO $p\n";
                            }
                        }
                    }
                }
            }
            else{
                #NADA
            }
        }
        else{
            print "catastophe $go" if ($save =~m/GO:0016645/);
            ($term)=$t[1];
            ($type)=$t[-1];
            ($annotated)=$t[2];
            if ($annotated <= $min_gene){
                #Do nothing
                #print "$go" if ($save =~m/GO:0016645/);
            }
            else{
                #print "dfasfef $go" if ($save =~m/GO:0016645/);
                if ($annotated >= $max_gene){
                   #Do nothing
                }
                else{
                    $p;
                    $e;
                    $log_p;
                    if ($fold_enrich >= $enrich){
                        #$e=$t[-2];
                        $p=$t[$pval_choice_col];
                        $log_p=sprintf("%.2f",-log($p)/log(10));
                        if ($p<$min_p){
                            $terms{$go}=$term;
                            push(@{$order{$type}},"$log_p=$go");
                        }
                        else{
                            #print "NO $p\n";
                        }
                    }
                }
            }
        }
        
    }



    @cols;
    $Good_to_Go=0;

    #print gap for the key
    push (@cols,$col_MF);
    print O "BLANK\t\t\t\n";
    push (@cols,$col_MF);
    print O "BLANK\t\t\t\n";

    @t=sort{$b<=>$a}(@{$order{CC}});
    @t=reverse @t;
    my $len_1=scalar(@t);
    foreach $t (@t){
        if (($len_1 <= $num_results) || !$num_results){
        ($p,$go)=$t=~/(.+?)\=(.+)/;
        $a++;
        print O "CC\t$go\t$terms{$go}\t$p\n";
        push (@cols,$col_CC);
        $Good_to_Go=1;
        }
        $len_1--;
    }
    push (@cols,$col_CC);
    print O "BLANK\t\t\t\n";
    #print "CC $ARGV[0] $Good_to_Go\n";

    @t=sort{$b<=>$a}(@{$order{MF}});
    @t=reverse @t;
    my $len_2=scalar(@t);
    foreach $t (@t){
        if (($len_2 <= $num_results) || !$num_results){
        ($p,$go)=$t=~/(.+?)\=(.+)/;
        $a++;
        print O "MF\t$go\t$terms{$go}\t$p\n";
        push (@cols,$col_MF);
        $Good_to_Go=1;
        }
        $len_2--;
    }
    push (@cols,$col_MF);
    print O "BLANK\t\t\t\n";
    #print "MF $ARGV[0] $Good_to_Go\n";

    @t=sort{$b<=>$a}(@{$order{BP}});
    @t=reverse @t;
    my $len_3=scalar(@t);
    foreach $t (@t){
        if (($len_3 <= $num_results) || !$num_results){
        ($p,$go)=$t=~/(.+?)\=(.+)/;
        $a++;
        print O "BP\t$go\t$terms{$go}\t$p\n";
        push (@cols,$col_BP);
        $Good_to_Go=1;
        }   
        $len_3--;
    }

    #print "BP $ARGV[0] $Good_to_Go\n";

    my $length_result=scalar(@cols);
    my $number_of_bars_total=$num_results*3;
    my $fract=$length_result/$number_of_bars_total;



    $plot_file="TopGO_Pval_barplot_"."$root.pdf";

    my $RCOMMAND="ACTUAL_R_CODE";
    open(my $out_Rcmds,   "> $RCOMMAND")  or die "error opening $RCOMMAND. $!";
    my $height=2+(8*$fract);
    print $out_Rcmds "pdf (\"$plot_file\", width=10, height=$height)\n";
    print $out_Rcmds "par(las=2)\n";
    print $out_Rcmds "BP<-read.table(\"$out_file\", h\=T,sep=\"\t\")\n";
    print $out_Rcmds "par(mar=c(5,25,4,10))\n";
    my @split_col=split("_", $ARGV[0]);

    my $palette_here=join ("\"\,\"", @cols);
    #my $Types=join(/","/, @cols);
    print $out_Rcmds "barplot(BP\$LogP.$type_p, names.arg = BP\$TERM, horiz=TRUE, xlab = \"log pvalues (correction=",$correction,")\", main= \"$root\", cex.main=2,  width = 1, col=c(\"$palette_here\"))\n";

    print $out_Rcmds "legend(\"bottomright\", title=\"GO Categories\", c(\"BP\",\"MF\",\"CC\"), ".
        "fill=c(\"$col_BP\", \"$col_MF\", \"$col_CC\"), horiz=TRUE, cex=0.8)\n";
    print $out_Rcmds "dev.off()\n";

    #my $length_out=`grep [MF|CC|BF] $out_file`;

    #print "DONE\n";
    if (!$Good_to_Go){
        print STDERR "FAILED, TopGO_Pval_barplot_"."$root.pdf.\nNo significant hits with these thresholds. Try lowering pval or removing correction\n";
    }
    else{
        `R --vanilla <ACTUAL_R_CODE> output.ofthis.test`;
        print STDERR "Plotted TopGO_Pval_barplot_"."$root.pdf\n";
        if ($open_me){
            `open $plot_file`;
        }
    }


    #Tidy up
    `rm -f $out_file`;

}
