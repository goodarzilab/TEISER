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

#use PBS ;
use threads;
use strict;
use Sets;
use Table;
use Getopt::Long;
use Data::Dumper;
use Thread::Semaphore;

if (@ARGV == 0) {
  die "Usage: perl teiser.pl --expfile=FILE --exptype=TXT --species=SP\n";
}

my @arg_copy = @ARGV;

Getopt::Long::Configure("pass_through");

my $expfile           = undef ;
my $platform          = undef ;
my $species           = undef ;
my $seedfiles         = undef ;
my $fastafile_up      = undef ;
my $fastafile_dn      = undef ;
my $fastafile_ort_up  = undef ;
my $fastafile_ort_dn  = undef ;
my $maxthreads        = 24;


my $quantized     = 1 ;
my $exptype          = "discrete" ;

my $do_up         = 1 ; 
my $do_dn         = 1 ;

my $submit        = 0 ;

GetOptions ( 'expfile=s'          => \$expfile,
	     'species=s'          => \$species,
	     'exptype=s'          => \$exptype,
	     'seedfiles=s'        => \$seedfiles,
	     'fastafile_up=s'     => \$fastafile_up,
	     'fastafile_dn=s'     => \$fastafile_dn,
	     'fastafile_ort_up=s' => \$fastafile_ort_up,
	     'fastafile_ort_dn=s' => \$fastafile_ort_dn,
	     'do_up=i'            => \$do_up,
	     'do_dn=i'            => \$do_dn,
	     'submit=s'           => \$submit,) ;

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

my $sem = Thread::Semaphore->new($maxthreads);

##################Phase 1####################
my $pwd     = `pwd`; $pwd =~ s/\n//;

my $fn = Sets::filename($expfile);
my $maindir = $expfile."_META" ;
mkdir $maindir if (! -d $maindir) ;
system("cp $expfile $maindir");

my $try ;
$try->{UP} = $do_up ;
$try->{DN} = $do_dn ;
my $hfasta ;
$hfasta->{UP}    = $fastafile_up ;
$hfasta->{DN}    = $fastafile_dn ;
$hfasta->{ORTUP} = $fastafile_ort_up ;
$hfasta->{ORTDN} = $fastafile_ort_dn ;


my @up_dn_jobs;
foreach my $do (keys %$try){
  print "Beginning loop for $do";
  my @threads = ();
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

    #$pbs->setScriptName("$t_expfile.script");
    #if ($platform eq 'tcluster') {
    #  $pbs->addCmd("setenv TEISERDIR $teiserdir");	      
    #} else {
    #  $pbs->addCmd("export TEISERDIR=$teiserdir") ;
    #}
    #$pbs->addCmd("cd $pwd") ;

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

    push @threads, async {$sem->down; system("perl $teiserdir/teiser.pl $args --seedfile=$seed --domireport=0 --dopagerun=0 --domotifint=0 --dodrawmatrix=0 --dodrawmimatrix=0 --dodrawpagematrix=0 --dodrawmotifs=0 --doconservation=0 --type=$do"); $sem->up} ;
    print "perl $teiserdir/teiser.pl $args --seedfile=$seed --domireport=0 --dopagerun=0 --domotifint=0 --dodrawmatrix=0 --dodrawmimatrix=0 --dodrawpagematrix=0 --dodrawmotifs=0 --doconservation=0 --type=$do";
  }

  foreach my $thread (@threads) {
    print "Waiting on thread $thread...";
    $thread->join;
  }

  print "\nDone with Phase 1!\n";

  ##################Phase 2##################
  system("cp $expfile $metadir");
  my @threadsp2 = ();

  #my $pbs = PBS->new ;
  #$pbs->setScriptName("$metadir/teiser_nondis_$do.script");

  #$pbs->setPlatform($platform) if (defined($platform));

  #if ($platform eq 'tcluster') {
  #  $pbs->addCmd("setenv TEISERDIR $teiserdir");	      
  #} else {
  #  $pbs->addCmd("export TEISERDIR=$teiserdir") ;
  #}


  #$pbs->addCmd("cd $pwd") ;

  my $threadcmd = "export TEISERDIR=$teiserdir";

  for (my $i=0 ; $i<@$seedfns ; $i++){
      if ($ext eq ""){
	  $threadcmd .= "; rm $maindir/$fn\_s$i.$do";
      }else{
	  $threadcmd .= "; rm $maindir/$fn\_s$i.$do.$ext";
      }
  }

  my $extname = "$do.$ext" ;
  $extname = "$do" if ($ext eq "") ;
  my $motiffileout = "$metadir/$fn.$ext.seed.bin" ;
  $motiffileout = "$metadir/$fn.seed.bin" if ($ext eq "") ;
  $threadcmd .= "; echo \"combining motifs.\"";
  my %PARAMS = ("maindir"        => $maindir,
		"metadir"         => $metadir,
		"filename"        => $fn,
		"extname"         => $extname,
		"suffix"          => "_TEISER",
		"motifoutfile"    => $motiffileout,
		"nseedfiles"      => scalar(@$seedfns),
	        "mode"            => 1, ) ;

  my $cmd = &get_cmd_combine(\%PARAMS);
  $threadcmd .= "; $cmd";
  
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
  $threadcmd .= "; perl $teiserdir/teiser.pl $args --domifind=0 --suffix=$do --type=$do" ;

  #if ($submit==0) {
  #  $pbs->execute ;
  #} elsif ($submit==1) {
  #  foreach my $id (@dep_jobs) {
  #    $pbs->addDepJob($id);
  #  }
#
#    my $jobid = $pbs->submit ;
#    push(@up_dn_jobs, $jobid) ;
#    print "Submitted job $jobid.\n";
#  }

  push @threadsp2, async {$sem->down; system($threadcmd); $sem->up};
  print $threadcmd;

  foreach my $thread (@threadsp2) {
    print "Waiting on thread $thread...";
    $thread->join;
  }
}

print "\nALL DONE WITH PART 2\n";

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

#my $pbs = PBS->new ;
my $finalcmd = "";
#$pbs->setScriptName("$metadir/teiser_nondis.script");

#$pbs->setPlatform($platform) if (defined($platform));

#if ($platform eq 'tcluster') {
#  $pbs->addCmd("setenv TEISERDIR $teiserdir");	      
#} else {
#  $pbs->addCmd("export TEISERDIR=$teiserdir") ;
$finalcmd .= "export TEISERDIR=$teiserdir; cd $pwd; echo \"combining results.\"";


#$pbs->addCmd("echo \"combining results.\""); 
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
$finalcmd .= "; $cmd";
#$pbs->addCmd($cmd);

my $order=0 ;
if ($exptype eq "continuous"){
    $order=1 ;
}

#$pbs->addCmd("echo \"drawing the final matrix.\"");
$finalcmd .= "; echo \"drawing the final matrix.\"";
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
#$pbs->addCmd($cmd);
$finalcmd .= "; $cmd";

my $cmd = &get_cmd_drawhtmlmatrix(\%PARAMS);
#$pbs->addCmd($cmd);
#$pbs->addCmd($cmd);
$finalcmd .= "; $cmd";

$finalcmd .= "; echo \"drawing motifs.\"";
my %PARAMS = ("summaryfile"     => "$metadir/$fn.$ext.summary.short") ;
%PARAMS = ("summaryfile"     => "$metadir/$fn.summary.short") if ($ext eq "") ;

my $cmd = &get_cmd_drawmotifs(\%PARAMS);
$finalcmd .= "; $cmd";


#$pbs->addCmd("echo \"drawing final mimatrix.\"");
$finalcmd .= "; echo \"drawing final mimatrix.\"";
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
$finalcmd .= "; $cmd";

#if ($submit==0) {
system($finalcmd);
print $finalcmd;
#} elsif ($submit==1) {
#  foreach my $id (@up_dn_jobs) {
#    $pbs->addDepJob($id);
#  }
  
#  my $jobid = $pbs->submit ;
#  print "Submitted job $jobid.\n";
#}

print "ALL DONE WITH EVERYTHING\n";
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
