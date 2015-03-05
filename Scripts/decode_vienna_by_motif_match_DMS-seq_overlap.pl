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
use Stat ;

my $dmsinvivofile = "/Users/hani/Life/Projects/Applications/TEISERv1.0/TEISER_Data/human/GSM1297493_K562vivo300DMS.txt" ;
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

my $viennafile = shift @ARGV ;

open I, "< $viennafile" or die ;
while (my $l=<I>){
	$l =~ s/\s+$// ;
	$l =~ s/^>// ;
	my ($n, $ref, $seq, $ss, $fold) = split(/\t/,$l) ;
	$seq =~ s/T/U/gi ;
	my @a = split(/_/, $n) ;
	my $ref = $a[0]."_".$a[1] ;
	next if (! defined $dms->{$ref}) ;
	#my $ref = $a[0] ;
	my $mseq = $a[2] ;
	my $mss = $a[3] ;
	my $pos = $a[5]+1-10 ;
	my $lcnt=0 ;
		
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
	my $iS ; my $iL ; 
	my $dS ; my $dL ; 
	my $total_a = Sets::maxInArray(\@iA) ;
	my $total_c = Sets::maxInArray(\@iC) ;
	my $total_cnt = $total_a + $total_c ;
	next if $total_a==0 ;
	next if $total_c==0 ;
	for (my $i=0 ; $i<length($seq) ; $i++) {
		$iA[$i]/=$total_a ;
		$iC[$i]/=$total_c ;
		$iSig[$i] = $iA[$i]+$iC[$i] ;
		my $char = substr($ss, $i, 1) ;
		if ($char eq "."){
			$iL+= $iSig[$i];
		}else{
			$iS+= $iSig[$i];
		}
	}
	my $total_a = Sets::maxInArray(\@dA) ;
	my $total_c = Sets::maxInArray(\@dC) ;
	next if $total_a==0 ;
	next if $total_c==0 ;
	for (my $i=0 ; $i<length($seq) ; $i++) {
		$dA[$i]/=$total_a ;
		$dC[$i]/=$total_c ;
		$dSig[$i] = $dA[$i]+$dC[$i] ;
		my $char = substr($ss, $i, 1) ;
		if ($char eq "."){
			$dL+= $dSig[$i];
		}else{
			$dS+= $dSig[$i];
		}
	}

	my ($r,$z,$p) = Stat::pearson(\@iSig, \@dSig) ;
	my $iG = Stat::calc_gini(@iSig) ;
	my $dG = Stat::calc_gini(@dSig) ;
	my $deltaG = $iG-$dG ;

	my $iR = ($iL+1e-5)/($iS+1e-5) ;
	my $dR = ($dL+1e-5)/($dS+1e-5) ;
	
	#if ($r<0.8 and $deltaG>-0.01 and $total_cnt>=10 and $dR<$iR){
	if ($r<0.8 and $deltaG>-0.01 and $total_cnt>=10 and $dR<$iR){
		print $n, "\t", $ref, "\t", $seq, "\t", $ss, "\t", $r, "\t", $deltaG, "\t", $iR, "\t", $dR,"\t", $total_cnt, "\n" ;
		print join("\t", split('',$seq)), "\n" ;
		print join("\t", split('',$ss)), "\n" ;
		print join("\t", @iSig), "\n" ;
		print join("\t", @dSig), "\n" ;
		print "\"", join(";", @iSig), "\"\n" ;
		my $val = join(";", @iSig) ;
		print ("java -Xmx8g -cp /Users/hani/Life/Projects/Applications/VARNAv3-9.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN $seq -structureDBN \"$ss\" -o \"/Users/hani/Desktop/$n.eps\" -baseOutline \"#000000\" -bpStyle \"simple\" -bp \"#000000\" -rotation -30 -colorMap \"$val\" -colorMapStyle \"0:#FFFFFF;0.0001:#4c4cff;0.4:#FFFF00;1:#FF0000\" -bpStyle line -bp \"#000000\"") ;
		system("java -Xmx8g -cp /Users/hani/Life/Projects/Applications/VARNAv3-9.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN $seq -structureDBN \"$ss\" -o \"/Users/hani/Desktop/$n.eps\" -baseOutline \"#000000\" -bpStyle \"simple\" -bp \"#000000\" -rotation -30 -colorMap \"$val\" -colorMapStyle \"0:#FFFFFF;0.0001:#4c4cff;0.4:#FFFF00;1:#FF0000\" -bpStyle line -bp \"#000000\"") ;

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
