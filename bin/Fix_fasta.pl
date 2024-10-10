#!/usr/bin/perl
use strict;
use warnings;

my $line;
while ($line = <ARGV>) {
	chomp $line;
	my $line2= <ARGV>;
	chomp $line2;
    if ($line2 =~ m/Sequence unavailable/){
    	#do nothing
    }
    else{
    	print ">$line\n$line2\n";
    }
}
