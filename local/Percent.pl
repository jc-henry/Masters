#!/usr/bin/perl 
use strict;
use warnings;

if($#ARGV !=2){
print "usage: Percent [input directory] [output directory] [percentage]\n";
exit;
}
my $input=$ARGV[0];
my $output=$ARGV[1];
my $percentage=$ARGV[2];

my @seqs;
my @names;
opendir(DH,$input) or die("Could not open input directory: $input\n");
my @inputfiles=readdir(DH);
closedir(DH);

opendir(DH,$output) or die("Could not open output directory: $output\n");
closedir(DH);

foreach my $inputfile(@inputfiles){
next if($inputfile=~/^\.$/);
next if($inputfile=~/^\.\.$/);
if(-e "$output/$inputfile"){
	print("File $inputfile does already exist in $output. Skipping\n");
	next;
}
@seqs=();
@names=();
my $whole="";
my $max=0;
open(Doc,"$input/$inputfile");
while(<Doc>){
my $line=$_;
if($line !~ />/){
	push(@seqs,$line);
	if(length($line)>$max){
		$max=length($line);
	}
}
else{
	push(@names,$line);
}
}
close Doc;
my $cut=0.01*$percentage*$max;
for(my $j=0;$j<@seqs;$j++){
	if(length($seqs[$j])>=$cut){
		$whole=$whole.$names[$j].$seqs[$j];
	}
}
my $file="$output/$inputfile";
open(FILE,">> $file");
print FILE $whole;
close FILE;
}

