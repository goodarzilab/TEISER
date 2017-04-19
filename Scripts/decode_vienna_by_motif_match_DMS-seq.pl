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
use Stat ;

my $dmsinvivofile = "/Users/hani/Life/Projects/Applications/TEISERv1.0/TEISER_Data/human/GSM1297493_K562vivo300DMS.txt" ;
#my $dmsinvivofile = "/Users/hani/Life/Projects/Applications/TEISERv1.0/TEISER_Data/human/GSM1297494_K562vivo400DMS.txt" ;
my $dmsdenaturefile = "/Users/hani/Life/Projects/Applications/TEISERv1.0/TEISER_Data/human/GSM1297495_K562denatured.txt" ;
my $dms;
open I, "< $dmsinvivofile" or die ;
while (my $l=<I>){
	$l =~ s/\s+$// ;
	my ($ref, $pos, $cnt) = split(/\t/, $l) ;
	$ref =~ s/\S+ref\|// ;
	$ref =~ s/\.\d+\|\s+// ;
	$pos =~ s/\s+$// ;
	$dms->{$ref}{$pos} = $cnt ;
}
my $den;
open I, "< $dmsdenaturefile" or die ;
while (my $l=<I>){
	$l =~ s/\s+$// ;
	my ($ref, $pos, $cnt) = split(/\t/, $l) ;
	$ref =~ s/\S+ref\|// ;
	$ref =~ s/\.\d+\|\s+// ;
	$pos =~ s/\s+$// ;
	$den->{$ref}{$pos} = $cnt ;
}

my $reportfile = shift @ARGV ;
my $fa = Fasta->new;   # open seqs
$fa->setFile($reportfile);

while (my $a_ref = $fa->nextSeq()) {
	my ($l, $seq) = @$a_ref;
	$l =~ s/\s+$// ;
	$l =~ s/^>// ;
	my $n = $l ;
	my @a = split(/_/, $l) ;
	my $ref = $a[0]."_".$a[1] ;
	next if (! defined $dms->{$ref}) ;
	#my $ref = $a[0] ;
	my $mseq = $a[2] ;
	my $mss = $a[3] ;
	my $pos = $a[5]+1-10 ;
	
	$seq =~ s/U/T/g ;

	
	my @iA ; my @iC ;
	my @dA ; my @dC ;
	for (my $i=0 ; $i<length($seq) ; $i++) {
		my $char = substr($seq, $i, 1) ;
		my $loc = $pos+$i ;
		if ((defined $dms->{$ref}{$loc}) and ($char eq 'A' or $char eq 'C')) {
			#print $l, "\t", $ref, "\t", $seq, "\t", $pos, "\t", $loc, "\t", $char, "\t", $dms->{$ref}{$loc}, "\n" ; <STDIN> ;
			if ($char eq 'C'){
				$iC[$i]=$dms->{$ref}{$loc} ;
				$dC[$i]=$den->{$ref}{$loc} ;
			}else{
				$iA[$i]=$dms->{$ref}{$loc} ;
				$dA[$i]=$den->{$ref}{$loc} ;
			}
		}else{
			$iA[$i]=0 ;
			$iC[$i]=0 ;
			$dA[$i]=0 ;
			$dC[$i]=0 ;
		}
	}
	my @iSig ;
	my @dSig ;
	my $total_a = Sets::maxInArray(\@iA) ;
	my $total_c = Sets::maxInArray(\@iC) ;
	my $total_cnt = $total_a + $total_c ;
	next if $total_a==0 ;
	next if $total_c==0 ;
	for (my $i=0 ; $i<length($seq) ; $i++) {
		$iA[$i]/=$total_a ;
		$iC[$i]/=$total_c ;
		$iSig[$i] = $iA[$i]+$iC[$i] ;
	}
	my $total_a = Sets::maxInArray(\@dA) ;
	my $total_c = Sets::maxInArray(\@dC) ;
	next if $total_a==0 ;
	next if $total_c==0 ;
	for (my $i=0 ; $i<length($seq) ; $i++) {
		$dA[$i]/=$total_a ;
		$dC[$i]/=$total_c ;
		$dSig[$i] = $dA[$i]+$dC[$i] ;
	}

	my ($r,$z,$p) = Stat::pearson(\@iSig, \@dSig) ;
	my $iG = Stat::calc_gini(@iSig) ;
	my $dG = Stat::calc_gini(@dSig) ;
	my $deltaG = $iG-$dG ;
	#print Dumper $r, $deltaG;

	#if ($r<0.6 and $deltaG>0.1){
	print $l, "\t", $ref, "\t", $seq, "\t", $r, "\t", $deltaG, "\t", $total_cnt, "\n" ;
	#}
	
	
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
