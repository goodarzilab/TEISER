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

if (@ARGV == 0) {
  die "Usage: perl teiser.pl --expfile=FILE --exptype=TXT --species=SP\n";
}

my @arg_copy = @ARGV;

Getopt::Long::Configure("pass_through");

my $expfile           = undef ;
my $phase             = undef ;
my $fold_cv           = 3 ;

GetOptions ( 'expfile=s'          => \$expfile,
	           'phase=i'            => $phase,) ;


#############phase I: make train/test set (3-fold cross-validation)#############
if ($phase == 1){
  my $fn = Sets::filename($expfile);
  my $maindir = $expfile."_XV" ;
  mkdir $maindir if (! -d $maindir) ;
  system("cp $expfile $maindir");
  my $ta = Table->new ;
  $t_expfile = "$maindir/$fn" ;
  my @files ;
  foreach my ($i=0 ; $i<$fold_cv ; $i++){
    my $trainfile = $t_expfile ;
    $trainfile =~ s/\.\S\S\S$/_train_set_$i.txt/ ;
    my $testfile = $t_expfile ;
    $testfile =~ s/\.\S\S\S$/_test_set_$i.txt/ ;
  }
  $ta->loadFile($expfile) ;
}

if ($exptype eq "discrete") {
  $quantized = 1 ;
} elsif ($exptype eq "continuous"){
  $quantized = 0 ;
}

if (defined($species)) {
  #  read species file
    my $species_data = readSpeciesData($species, $teiserdir);
    $seedfiles         = $species_data->{"seedfiles"} if (! (defined($seedfiles)));
    $fastafile_up      = $species_data->{"fastafile_up"} if (! (defined($fastafile_up)));
    $fastafile_dn      = $species_data->{"fastafile_dn"} if (! (defined($fastafile_dn)));
    $fastafile_ort_up  = $species_data->{"fastafile_ort_up"} if (! (defined($fastafile_ort_up)));
    $fastafile_ort_dn  = $species_data->{"fastafile_ort_dn"} if (! (defined($fastafile_ort_dn))); 
}

##################Phase 1####################
my $pwd     = `pwd`; $pwd =~ s/\n//;


my $try ;
$try->{UP} = $do_up ;
$try->{DN} = $do_dn ;
my $hfasta ;
$hfasta->{UP}    = $fastafile_up ;
$hfasta->{DN}    = $fastafile_dn ;
$hfasta->{ORTUP} = $fastafile_ort_up ;
$hfasta->{ORTDN} = $fastafile_ort_dn ;


my @up_dn_jobs ;
foreach my $do (keys %$try){
  next if ($try->{$do} == 0) ;

  my $fastafile     = $hfasta->{$do} ;
  my $fastafile_ort = $hfasta->{"ORT$do"} ;

  my @dep_jobs ;

  my $dir = Sets::dirname($expfile) ;
  my $fn = Sets::filename($expfile) ;
  my $ext = "";
  if ($fn =~ /\./){
      my @e = split(/\./, $fn) ;
      $ext = pop(@e) ;
      $fn = join(".", @e) ;
  }

  my $metadir = $maindir."/$fn.$ext\_TEISER" ;
  if ($ext eq ""){
      $metadir = $maindir."/$fn\_TEISER" ;
  }
  mkdir $metadir if (! -d $metadir) ;
  $metadir .= "/$do";
  mkdir $metadir if (! -d $metadir) ;

  my $ta = Table->new ;
  $ta->loadFile($seedfiles) ;
  my $seedfns = $ta->getArray ;

  for (my $i=0 ; $i<@$seedfns ; $i++){
    my $seed = $teiserdir."/".$seedfns->[$i][0] ;

    my @arg = @arg_copy;
  
    my $t_expfile ;
    if ($ext eq ""){
	system("cp $expfile $maindir/$fn\_s$i.$do");
	$t_expfile = "$maindir/$fn\_s$i.$do";
        print "Submitting file $t_expfile\n" ;
    }else{
	system("cp $expfile $maindir/$fn\_s$i.$do.$ext");
	$t_expfile = "$maindir/$fn\_s$i.$do.$ext";
	print "Submitting file $t_expfile\n" ;
    }

    my $pbs = PBS->new ;
    $pbs->setScriptName("$t_expfile.script");
    if ($platform eq 'tcluster') {
      $pbs->addCmd("setenv TEISERDIR $teiserdir");	      
    } else {
      $pbs->addCmd("export TEISERDIR=$teiserdir") ;
    }
    $pbs->addCmd("cd $pwd") ;

    for (my $i=0 ; $i<@arg ; $i++) {
      if ($arg[$i] =~ /--submit/) {
	$arg[$i] = "--submit=0" ;
      }
      if ($arg[$i] =~ /--expfile/) {
	$arg[$i] = "--expfile=$t_expfile" ;
      }
      if ($arg[$i] =~ /--fastafile_rna/) {
	$arg[$i] = "--fastafile_rna=$fastafile" ;
      }
    }

    my $args = join(" ", @arg);
    if (! ($args =~ /--fastafile_rna/)) {
      $args .= " --fastafile_rna=$fastafile" ;
    }

    $pbs->addCmd("perl $teiserdir/teiser.pl $args --seedfile=$seed --domireport=0 --dopagerun=0 --domotifint=0 --dodrawmatrix=0 --dodrawmimatrix=0 --dodrawpagematrix=0 --dodrawmotifs=0 --doconservation=0 --type=$do") ;
    my $teiser_jobid ;
    if ($submit==0){
      $pbs->execute ;
    }elsif ($submit==1){
      my $teiser_jobid = $pbs->submit ;
      push(@dep_jobs, $teiser_jobid);
      print "Submitted job $teiser_jobid.\n";
    }
  }

  ##################Phase 2##################
  system("cp $expfile $metadir");

  my $pbs = PBS->new ;
  $pbs->setScriptName("$metadir/teiser_nondis_$do.script");

  $pbs->setPlatform($platform) if (defined($platform));

  if ($platform eq 'tcluster') {
    $pbs->addCmd("setenv TEISERDIR $teiserdir");	      
  } else {
    $pbs->addCmd("export TEISERDIR=$teiserdir") ;
  }


  $pbs->addCmd("cd $pwd") ;
  for (my $i=0 ; $i<@$seedfns ; $i++){
      if ($ext eq ""){
	  $pbs->addCmd("rm $maindir/$fn\_s$i.$do") ;
      }else{
	  $pbs->addCmd("rm $maindir/$fn\_s$i.$do.$ext") ;
      }
  }

  my $extname = "$do.$ext" ;
  $extname = "$do" if ($ext eq "") ;
  my $motiffileout = "$metadir/$fn.$ext.seed.bin" ;
  $motiffileout = "$metadir/$fn.seed.bin" if ($ext eq "") ;
  $pbs->addCmd("echo \"combining motifs.\""); 
  my %PARAMS = ("maindir"        => $maindir,
		"metadir"         => $metadir,
		"filename"        => $fn,
		"extname"         => $extname,
		"suffix"          => "_TEISER",
		"motifoutfile"    => $motiffileout,
		"nseedfiles"      => scalar(@$seedfns),
	        "mode"            => 1, ) ;

  my $cmd = &get_cmd_combine(\%PARAMS);
  $pbs->addCmd($cmd);
  
  my @arg = @arg_copy;
  for (my $i=0 ; $i<@arg ; $i++) {
    if ($arg[$i] =~ /--submit/) {
      $arg[$i] = "--submit=0" ;
    }
    if ($arg[$i] =~ /--expfile/) {
      $arg[$i] = "--expfile=$maindir/$fn.$ext" ;
      $arg[$i] = "--expfile=$maindir/$fn" if ($ext eq "") ;
    }
    if ($arg[$i] =~ /--seedfile/) {
      $arg[$i] = "--seedfile=$metadir/$fn.$ext.seed.bin" ;
      $arg[$i] = "--seedfile=$metadir/$fn.seed.bin" if ($ext eq "") ;
    }
    if ($arg[$i] =~ /--fastafile_rna/) {
      $arg[$i] = "--fastafile_rna=$fastafile" ;
    }
    if ($arg[$i] =~ /--fastafile_ort/) {
      $arg[$i] = "--fastafile_ort=$fastafile_ort" ;
    }
  }
  my $args = join(" ", @arg);
  if (! ($args =~ /--seedfile/)) {
      if ($ext eq ""){
	  $args .= " --seedfile=$metadir/$fn.seed.bin" ;
      }else{
	  $args .= " --seedfile=$metadir/$fn.$ext.seed.bin" ;
      }
  }
  if (! ($args =~ /--fastafile_rna/)) {
    $args .= " --fastafile_rna=$fastafile" ;
  }
  if (! ($args =~ /--fastafile_ort/)) {
    $args .= " --fastafile_ort=$fastafile_ort" ;
  }
  $pbs->addCmd("perl $teiserdir/teiser.pl $args --domifind=0 --suffix=$do --type=$do") ;

  if ($submit==0) {
    $pbs->execute ;
  } elsif ($submit==1) {
    foreach my $id (@dep_jobs) {
      $pbs->addDepJob($id);
    }

    my $jobid = $pbs->submit ;
    push(@up_dn_jobs, $jobid) ;
    print "Submitted job $jobid.\n";
  }
}

my $dir = Sets::dirname($expfile) ;
my $fn = Sets::filename($expfile) ;
my $ext = "";
if ($fn =~ /\./){
    my @e = split(/\./, $fn) ;
    $ext = pop(@e) ;
    $fn = join(".", @e) ;
}

my $metadir = $maindir."/$fn.$ext\_TEISER/UP_DN" ;
$metadir = $maindir."/$fn\_TEISER/UP_DN" if ($ext eq "") ;
mkdir $metadir if (! -d $metadir) ;

mkdir "$metadir/Motifs" if (! -d "$metadir/Motifs") ;

system("cp $expfile $metadir");

my $pbs = PBS->new ;
$pbs->setScriptName("$metadir/teiser_nondis.script");

$pbs->setPlatform($platform) if (defined($platform));

if ($platform eq 'tcluster') {
  $pbs->addCmd("setenv TEISERDIR $teiserdir");	      
} else {
  $pbs->addCmd("export TEISERDIR=$teiserdir") ;
}
$pbs->addCmd("cd $pwd") ;


$pbs->addCmd("echo \"combining results.\""); 
my %PARAMS = ("motifdir_up"    => "$maindir/$fn.$ext\_TEISER/UP/Motifs",
	      "motifdir_dn"     => "$maindir/$fn.$ext\_TEISER/DN/Motifs",
	      "summaryfile_up"  => "$maindir/$fn.$ext\_TEISER/UP/$fn.$ext.summary",
	      "summaryfile_dn"  => "$maindir/$fn.$ext\_TEISER/DN/$fn.$ext.summary",
	      "matrixfile"      => "$metadir/$fn.$ext.mimatrix",
	      "summaryfile"     => "$metadir/$fn.$ext.summary", );
if ($ext eq ""){
    %PARAMS = ("motifdir_up"    => "$maindir/$fn\_TEISER/UP/Motifs",
              "motifdir_dn"     => "$maindir/$fn\_TEISER/DN/Motifs",
              "summaryfile_up"  => "$maindir/$fn\_TEISER/UP/$fn.summary",
              "summaryfile_dn"  => "$maindir/$fn\_TEISER/DN/$fn.summary",
              "matrixfile"      => "$metadir/$fn.mimatrix",
	       "summaryfile"     => "$metadir/$fn.summary", ) ;
}

my $cmd = &get_cmd_final(\%PARAMS);
$pbs->addCmd($cmd);

my $order=0 ;
if ($exptype eq "continuous"){
    $order=1 ;
}

$pbs->addCmd("echo \"drawing the final matrix.\"");
my %PARAMS = ("pvmatrixfile_up"=> "$maindir/$fn.$ext\_TEISER/UP/$fn.$ext.matrix",
	      "pvmatrixfile_dn" => "$maindir/$fn.$ext\_TEISER/DN/$fn.$ext.matrix",
	      "consfile_up"     => "$maindir/$fn.$ext\_TEISER/UP/$fn.$ext.cons",
	      "consfile_dn"     => "$maindir/$fn.$ext\_TEISER/DN/$fn.$ext.cons",
	      "summaryfile"     => "$metadir/$fn.$ext.summary",
	      "expfile"         => "$metadir/$fn.$ext",
	      "order"           => $order,
	      "quantized"        => $quantized, ) ;
if ($ext eq ""){
    %PARAMS = ("pvmatrixfile_up"=> "$maindir/$fn\_TEISER/UP/$fn.matrix",
              "pvmatrixfile_dn" => "$maindir/$fn\_TEISER/DN/$fn.matrix",
              "consfile_up"     => "$maindir/$fn\_TEISER/UP/$fn.cons",
              "consfile_dn"     => "$maindir/$fn\_TEISER/DN/$fn.cons",
              "summaryfile"     => "$metadir/$fn.summary",
              "expfile"         => "$metadir/$fn",
              "order"           => $order,
	       "quantized"        => $quantized, ) ;
}

my $cmd = &get_cmd_drawmatrix(\%PARAMS);
$pbs->addCmd($cmd);

my $cmd = &get_cmd_drawhtmlmatrix(\%PARAMS);
$pbs->addCmd($cmd);

$pbs->addCmd("echo \"drawing motifs.\"");
my %PARAMS = ("summaryfile"     => "$metadir/$fn.$ext.summary.short") ;
%PARAMS = ("summaryfile"     => "$metadir/$fn.summary.short") if ($ext eq "") ;

my $cmd = &get_cmd_drawmotifs(\%PARAMS);
$pbs->addCmd($cmd);


$pbs->addCmd("echo \"drawing final mimatrix.\"");
my %PARAMS = ("pvmatrixfile"  => "$metadir/$fn.$ext.mimatrix",
	      "summaryfile"     => "$metadir/$fn.$ext.summary",
	      "expfile"         => "$metadir/$fn.$ext",
	      "quantized"       => $quantized,) ;
if ($ext eq ""){
    %PARAMS = ("pvmatrixfile"  => "$metadir/$fn.mimatrix",
              "summaryfile"     => "$metadir/$fn.summary",
              "expfile"         => "$metadir/$fn",
	       "quantized"       => $quantized,) ;
}

my $cmd = &get_cmd_drawmimatrix(\%PARAMS);
$pbs->addCmd($cmd);

if ($submit==0) {
  $pbs->execute ;
} elsif ($submit==1) {
  foreach my $id (@up_dn_jobs) {
    $pbs->addDepJob($id);
  }
  
  my $jobid = $pbs->submit ;
  print "Submitted job $jobid.\n";
}

sub readSpeciesData {
  my ($species, $teiserdir) = @_;
  
  my %H = ();
  open IN, "$teiserdir/TEISER_Data/species_data/$species" or die "No data file for $species: $teiserdir/TEISER_Data/species_data/$species.\n";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;    
    if ($a[1] =~ /^TEISER_Data/) {
      $a[1] = "$teiserdir/$a[1]";
    }
    $H{$a[0]} = $a[1];
  }  
  close IN;
  return \%H;
}

sub get_cmd_combine{
  my ($p) = @_;
  my $todo = "$programdir/combine_motifs -maindir $p->{maindir} -metadir $p->{metadir} -filename $p->{filename} -extname $p->{extname} -suffix $p->{suffix} -motifoutfile $p->{motifoutfile} -nseedfiles $p->{nseedfiles} -mode $p->{mode}" ;
  return $todo;
}

sub get_cmd_final{
  my ($p) = @_;
  my $todo = "$programdir/final_motif_motif_interaction -motifdir_up $p->{motifdir_up} -motifdir_dn $p->{motifdir_dn} -summaryfile_up $p->{summaryfile_up} -summaryfile_dn $p->{summaryfile_dn} -matrixfile $p->{matrixfile} -summaryfile $p->{summaryfile}" ;
  return $todo;
}

sub get_cmd_drawmatrix {
  my ($p) = @_;
  my $todo = "perl $scriptdir/teiser_draw_matrix.pl --pvmatrixfile_up=$p->{pvmatrixfile_up} --pvmatrixfile_dn=$p->{pvmatrixfile_dn} --consfile_up=$p->{consfile_up} --consfile_dn=$p->{consfile_dn} --summaryfile=$p->{summaryfile} --expfile=$p->{expfile} --quantized=$p->{quantized} --order=$p->{order}" ;
  return $todo;
}

sub get_cmd_drawhtmlmatrix {
    my ($p) = @_;
  my $todo = "perl $scriptdir/teiser_draw_matrix_html.pl --pvmatrixfile_up=$p->{pvmatrixfile_up} --pvmatrixfile_dn=$p->{pvmatrixfile_dn} --consfile_up=$p->{consfile_up} --consfile_dn=$p->{consfile_dn} --summaryfile=$p->{summaryfile} --expfile=$p->{expfile} --quantized=$p->{quantized} --order=$p->{order}" ;
    return $todo;
}

sub get_cmd_drawmimatrix {
  my ($p) = @_;
  my $todo = "perl $scriptdir/teiser_draw_mimatrix.pl --pvmatrixfile=$p->{pvmatrixfile} --summaryfile=$p->{summaryfile} --expfile=$p->{expfile} --quantized=$p->{quantized}" ;
  return $todo;
}

sub get_cmd_drawmotifs {
    my ($p) = @_;
    my $todo = "perl $scriptdir/teiser_draw_motifs.pl --summaryfile=$p->{summaryfile}" ;
    return $todo;
}
