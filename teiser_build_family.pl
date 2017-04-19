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
use File::Basename;

if (@ARGV == 0) {
  die "Usage: perl teiser_build_family.pl --expfile=FILE --index=INT --fastafile=FILE --submit=0/1 --type=UP/DN --minr=0.5\n";
}

my @arg_copy = @ARGV;

Getopt::Long::Configure("pass_through");

my $expfile           = undef ;
my $index             = undef ;
my $fastafile         = undef ;
my $type              = "DN" ;
my $minr              = 0.5 ;

my $quantized     = 1 ;
my $exptype       = "discrete" ;

my $submit        = 0 ;

GetOptions ( 'expfile=s'          => \$expfile,
	     'index=s'            => \$index,
	     'fastafile=s'        => \$fastafile,
	     'submit=s'           => \$submit,
       'type=s'             => \$type,
       'minr=s'             => \$minr) ;

##################Phase 1####################
my $pwd     = `pwd`; $pwd =~ s/\n//;

my $fn = Sets::filename($expfile);
my $maindir = $expfile."_META" ;
my $teiserdir = $expfile."_META/".$fn."_TEISER/".$type."/" ;

my $qfile = $expfile."_META/".$fn."_TEISER/".$type."/".$fn.".q" ;
my $optimfile = $expfile."_META/".$fn."_TEISER/".$type."/".$fn.".optim.bin" ;
my $indexfile = $expfile."_META/".$fn."_TEISER/".$type."/".$fn.".optim.index" ;
my $motiffile = $expfile."_META/".$fn."_TEISER/".$type."/".$fn.".optim.i".$index.".bin" ;
my $familyfile = $expfile."_META/".$fn."_TEISER/".$type."/".$fn.".optim.i".$index.".family.bin" ;
my $seedstxt = $expfile."_META/".$fn."_TEISER/".$type."/".$fn.".seeds.txt" ;
my $seedsfile = $expfile."_META/".$fn."_TEISER/".$type."/".$fn.".seeds.bin" ;
my $tmpdir  = $expfile."_META/".$fn."_TEISER/".$type."/tmp" ;
print("mkdir $tmpdir\n") ;
system("mkdir $tmpdir") ;

print("$programdir/get_motif_at_index -motiffile $optimfile -indexfile $indexfile -outfile $motiffile\n");
system("$programdir/get_motif_at_index -motiffile $optimfile -indexfile $indexfile -outfile $motiffile");

my $globtxt = $fn ;
$globtxt =~ s/\.txt/_s*.$type.txt_TEISER\/*.seed.bin/ ;
print ("ls -al $maindir/$globtxt | awk '{print \$9}' > $seedstxt\n");
system("ls -al $maindir/$globtxt | awk '{print \$9}' > $seedstxt");

print ("$programdir/combine_motifs_into_file -motiffile $seedstxt -outfile $seedsfile\n");
system ("$programdir/combine_motifs_into_file -motiffile $seedstxt -outfile $seedsfile");

my $splitseedsfile = $tmpdir."/$fn";
print ("$programdir/split_seeds_into_files -motiffile $seedsfile -outfile $splitseedsfile -count_per_file 5000\n");
system ("$programdir/split_seeds_into_files -motiffile $seedsfile -outfile $splitseedsfile -count_per_file 5000");

my @splitfiles = glob ("$splitseedsfile.split*.bin");

my @splitoutfiles ;
my @dep_jobs ;
for my $sfile (@splitfiles){
    my $pbs = PBS->new ;
    my $scriptname = $sfile;
    $scriptname =~ s/\.bin$/.script/ ;
    $pbs->setScriptName("$scriptname");
    $pbs->addCmd("export TEISERDIR=$teiserdir") ;
    $pbs->addCmd("cd $pwd") ;

    my $soutfile = $sfile;
    $soutfile =~ s/\.bin$/.condInfo/ ;
    push(@splitoutfiles, $soutfile."_0.seeds.bin") ;

    $pbs->addCmd("$programdir/mi_cond_inf_ratio -expfile $qfile -seedfile $sfile -motiffile $motiffile -quantized 1 -rna_fastafile $fastafile -motifoutfile $soutfile -minr $minr -dG_t -1") ;
    my $teiser_jobid ;
    if ($submit==0){
      $pbs->execute ;
    }elsif ($submit==1){
      my $teiser_jobid = $pbs->submit ;
      push(@dep_jobs, $teiser_jobid);
      print "Submitted job $teiser_jobid.\n";
    }
}

my $pbs = PBS->new ;
$pbs->setScriptName("$tmpdir/combine.script");
$pbs->addCmd("export TEISERDIR=$teiserdir") ;
$pbs->addCmd("cd $pwd") ;
$pbs->addCmd("$programdir/combine_motifs_from_files ".scalar(@splitfiles)." ".join(' ',@splitoutfiles)." $familyfile") ;

if ($submit==0) {
  $pbs->execute ;
} elsif ($submit==1) {
  foreach my $id (@dep_jobs) {
    $pbs->addDepJob($id);
  }
  my $jobid = $pbs->submit ;
  print "Submitted job $jobid.\n";
}
