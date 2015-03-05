#perl /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Scripts/design_oligos.pl /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Results/SNPs_in_3utrs_MDA_LM2_motif_filtered_kmers_and_structural_seeds_final_report.txt /Users/hani/Life/Projects/Databases/Transcript_info.txt /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Results/SNPs_in_3utrs_MDA_LM2_motif_filtered_kmers_and_structural_seeds_final_oligos.txt
use lib "/Users/hani/Life/Projects/Perl/Tools/";

use FileHandle;
use strict ;
use Getopt::Long ;
use Data::Dumper ;
use PostScript::Simple;
use Sets ;
use Table ;
use Subs ;
use Fasta ;

my $viennafile = shift @ARGV ;

open I, "< $viennafile" or die ;
while (my $l=<I>){
	$l =~ s/\s+$// ;
	$l =~ s/^>// ;
	my $n = $l ;
	my @a = split(/_/, $l) ;
	my $ref = $a[0]."_".$a[1] ;
	
	my $l=<I> ;
	$l =~ s/\s+$// ;
	my $seq = $l ;
	$seq =~ s/U/T/g ;

	my $l=<I> ;
	$l =~ s/\s+$// ;
	my $ss = $l ;
	$ss =~ s/\s.*$// ;
	
	my $l=<I> ;
	$l =~ /= (.*) kcal/ ;
	my $fe = $1 ;	
	print $n, "\t", $ref, "\t", $seq, "\t", $ss, "\t", $fe, "\n" ;
	my $l=<I> ;	
}
close(I) ;
