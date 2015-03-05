#!/usr/bin/perl

my $teiserdir ;
BEGIN{
    if ((!$ENV{TEISERDIR}) || ($ENV{TEISERDIR} eq '')) {
	$teiserdir="./" ;
	print "The TEISERDIR environment variable is not set. It is set to default.\n";
    }
    else{
	print "The TEISERDIR environment variable is".$ENV{TEISERDIR}."\n" ;
	$teiserdir = $ENV{TEISERDIR};
    }
}

my $programdir = $teiserdir."/Programs" ;
my $scriptdir  = $teiserdir."/Scripts" ;

my $pwd        = `pwd`; $pwd =~ s/\n//;
use lib "$teiserdir/Scripts";

use PBS ;
use strict;
use Sets;
use Table;
use Getopt::Long;
use Data::Dumper;

my $dirname = shift @ARGV ;
my $expfile = shift @ARGV ;

for (my $i=0 ; $i<19 ; $i++){
  my $pbs = PBS->new ;
  
  my $seedfile = "$dirname/filtered.%02d.bin"
  my $dataoutfile = "$dirname/filtered.%02d.txt"

  $pbs->setScriptName("$seedfile.script");
  $pbs->addCmd("export TEISERDIR=$teiserdir") ;
  $pbs->addCmd("cd $pwd") ;

  $pbs->addCmd("perl $programdir/mi_mi_seed.pl -seedfile $seedfile -dataoutfile $dataoutfile -expfile $expfile -quantized 1") ;
  my $teiser_jobid ;
  if ($submit==0){
    $pbs->execute ;
  }elsif ($submit==1){
    my $teiser_jobid = $pbs->submit ;
    push(@dep_jobs, $teiser_jobid);
    print "Submitted job $teiser_jobid.\n";
  }
}
