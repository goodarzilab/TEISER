#!/usr/bin/perl

use lib "$ENV{FIREDIR}/SCRIPTS";
if ((!$ENV{FIREDIR}) || ($ENV{FIREDIR} eq '')) {
  die "Please set the FIREDIR environment variable.\n";
}

#
#  INPUT : one set of aa sequences, one list of genomes
#

use MyBlast;
use Fasta;
use Sets;
use Sequence;
use Table;
use strict;


#
#  get a new sequence object to retrieve BLAST seqs
#


my $mb = MyBlast->new;
$mb->setVerbose(0);
$mb->setBlastProgram("blastp");
$mb->setEvalueThreshold((defined($ARGV[1])?$ARGV[1]:"1e-10"));
$mb->setMegablast(1);
$mb->setNbProcessors(2);
#$mb->setFilter(0);

#
#  traverse all the sequences in file
#

if (! -e "$ARGV[0].pin") {
  
  #system("formatdb -i $ARGV[0] -p T -o T") == 0 or die "Please format $ARGV[0] using formatdb -i $ARGV[0] -p T -o T\n";
  
}

my $s_tmpstore1 = Sets::getTempFile("/tmp/tmp1.seq");

my $fa = Fasta->new;
$fa->setFile($ARGV[0]);


my $h_ref_done = undef;
#if (defined($ARGV[1])) {
#  my $ta = Table->new;
#  $ta->loadFile($ARGV[1]);
#  $h_ref_done = $ta->getIndex(0);
#}



while (my $a_ref = $fa->nextSeq()) {
  
  my ($n, $s) = @$a_ref;

  if (defined($h_ref_done) && defined($h_ref_done->{ $n })) {
    next;
  }


  my @ss  = split /N+/, $s;
  my $max_length = 0;
  foreach my $sss (@ss) {
    if (length($sss) > $max_length) {
      $max_length = length($sss);
    }
  }
  
  if ($max_length < 28) {
    print "$n\n";
    next;
  }
    

  #my @ss  = split //, $s;
  #my @ssn = grep /N/, @ss;
  #if (scalar(@ssn) > 0.9 * scalar(@ss)) {
  #  print "$n\n";
  #  next;
  #}

  open SEQ, ">$s_tmpstore1";
  print SEQ ">$n\n$s\n\n";
  close SEQ;
  
  
  # use BLAST to align the current sequence against the reference sequence 
  $mb->setQueryDatabase($s_tmpstore1, "T");
  
  # set the database
  $mb->setDatabaseDatabase($ARGV[0], "T");
  
  my $a_hits = $mb->blastallMultiple;
  
  print "$n";
  
  foreach my $h (@$a_hits) {
    if ($n ne $h->{HIT_NAME}) {
      print "\t" . $h->{HIT_NAME};
      #foreach my $s ( @{ $h->{HSPS} } ) {
      #	print "  $s->{ALIGNLEN}\n";
      #      }
    }
  }

  print "\n";

}


unlink $s_tmpstore1;
