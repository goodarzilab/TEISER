#perl /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Scripts/design_oligos.pl /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Results/SNPs_in_3utrs_MDA_LM2_motif_filtered_kmers_and_structural_seeds_final_report.txt /Users/hani/Life/Projects/Databases/Transcript_info.txt /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Results/SNPs_in_3utrs_MDA_LM2_motif_filtered_kmers_and_structural_seeds_final_oligos.txt
use lib "/Users/hani/Life/Projects/Perl/Tools/";

use FileHandle;
use strict ;
use Getopt::Long ;
use Data::Dumper ;

my $viennafile = shift @ARGV ;

open I, "< $viennafile" or die ;
while (my $l=<I>){
	$l =~ s/\s+$// ;
	$l =~ s/^>// ;
	my $n = $l ;
	my @a = split(/_/, $l) ;
	#my $ref = $a[0]."_".$a[1]."_".$a[2]."_".$a[3] ;
	#my $ref = $a[0]."_".$a[1] ;
	#my $ref = $a[0]."_".$a[1] ;
	my $ref = $a[0] ;
	my $mseq = $a[1] ;
	my $mss = $a[2] ;

	my $lcnt=0 ;
	for (my $i=0 ; $i<length($mss) ; $i++){
		$lcnt++ if (substr($mss, $i, 1) eq ".") ;
		last if (substr($mss, $i, 1) eq "(") ;
	}
	my $rcnt=0 ;
	for (my $i=length($mss)-1 ; $i>=0 ; $i--){
		$rcnt++ if (substr($mss, $i, 1) eq ".") ;
		last if (substr($mss, $i, 1) eq ")") ;
	}
	$mseq = substr($mseq, $lcnt, length($mseq)-$rcnt-$lcnt) ;
	$mss = substr($mss, $lcnt, length($mss)-$rcnt-$lcnt) ;
	$mseq = &seq_to_re($mseq) ;
	
	my $l=<I> ;
	$l =~ s/\s+$// ;
	my $seq = $l ;
	$seq =~ s/U/T/g ;

	my $matseq = substr($seq, 10+$lcnt, length($mss)) ;
	#print Dumper $matseq, $mseq ; <STDIN> ;

	my $l=<I> ;
	$l =~ s/\s+$// ;
	my $ss = $l ;
	$ss =~ s/\s.*$// ;
	my $matss = substr($ss, 10+$lcnt, length($mss)) ;

	if ($matss eq $mss){
		my $l=<I> ;
		$l =~ /= (.*) kcal/ ;
		my $fe = $1 ;	
		print $n, "\t", $ref, "\t", $seq, "\t", $ss, "\t", $fe, "\n" ;
		my $l=<I> ;	
	}else{
		my $l=<I> ;
		$l =~ /= (.*) kcal/ ;
		my $fe = $1 ;	
		my $l=<I> ;	
	}
}
close(I) ;

sub seq_to_re {
  my $S = shift @_ ;
  $S =~ s/N/./g ;
  $S =~ s/Y/[TC]/g ;
  $S =~ s/R/[AG]/g ;
  $S =~ s/K/[TG]/g ;
  $S =~ s/M/[AC]/g ;
  $S =~ s/S/[GC]/g ;
  $S =~ s/W/[AT]/g ;
  $S =~ s/B/[GTC]/g ;
  $S =~ s/D/[GAT]/g ;
  $S =~ s/H/[ACT]/g ;
  $S =~ s/V/[GCA]/g ;

  return $S ;
}
