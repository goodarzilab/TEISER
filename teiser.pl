#!/usr/bin/perl

my $teiserdir ;
BEGIN{
    if ((!$ENV{TEISERDIR}) || ($ENV{TEISERDIR} eq '')) {
	$teiserdir="./" ;
	print "The TEISERDIR environment variable is not set. It is set to default.\n";
    }
    else{
	print "The TEISERDIR environment variable is ".$ENV{TEISERDIR}."\n" ;
	$teiserdir = $ENV{TEISERDIR};
    }
}

my $cmdline = "perl teiser.pl";
foreach my $r (@ARGV) {
  $cmdline .= " $r";
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

my @argv_copy = @ARGV;

my $expfile          = undef ;
my $goindexfile      = undef ;
my $gonamesfile      = undef ;
my $seedfile         = undef ;
my $fastafile_rna    = undef ;
my $fastafile_ort    = undef ;
my $homologyfile     = undef ;
my $species          = undef ;
my $exptype          = "discrete" ;
my $ebins            = undef ; 
my $quantized        = 1 ;
my $divbins          = 50.0 ;
my $mbins            = 2 ;
my $type             = "DN" ;

my $shuffle_mifind   = 1500000 ;
my $maxp_mifind      = 0.00000001 ;
my $shuffle_mimotif  = 1000000 ;
my $maxp_mimotif     = 0.0000001 ;
my $shuffle_page     = 10000 ;
my $maxp_page        = 0.005 ;
my $maxz_mioptimize  = -1 ;
my $maxz_mifind      = -1 ;
my $jn_t             = 8 ;
my $min_r            = 2.0 ;
my $dG_t             = -1.0 ;

my $clusters         = 5 ;
my $draw_min         = -10 ;
my $draw_max         = 10 ;
my $removecols_draw  = undef ;
my $colmap_matrix    = "$scriptdir/HEATMAPS/cmap_1.txt" ;
my $colmap_mimatrix  = "$scriptdir/HEATMAPS/cmap_2.txt" ;
my $colmap_page      = "$scriptdir/HEATMAPS/cmap_3.txt" ;
my $order            = 0 ;

my $dodrawmatrix       = 1 ;
my $dodrawmimatrix     = 1 ;
my $dodrawpagematrix   = 1 ;
my $dodrawmotifs       = 1 ;
my $doskipdiscovery    = 0 ;
my $domifind           = 1 ;
my $doskipoptimization = 0 ;
my $domioptimize       = 1 ;
my $doonlypositive     = 0 ;
my $domireport         = 1 ;
my $doremovedups       = 1 ;
my $dopagerun          = 1 ;
my $domotifind         = 1 ;
my $doconservation     = 1 ;

my $suffix           = "" ;
my $submit           = 0 ;
my $platform         = undef ;
my $walltime         = "72:00:00" ;
my $queue            = undef ;
my $jobid            = undef ;

GetOptions ('expfile=s'              => \$expfile,
	    'goindexfile=s'          => \$goindexfile,
	    'gonamesfile=s'          => \$gonamesfile,
	    'seedfile=s'             => \$seedfile,
	    'fastafile_rna=s'        => \$fastafile_rna,
	    'fastafile_ort=s'        => \$fastafile_ort,
	    'species=s'              => \$species,
	    'exptype=s'              => \$exptype,
	    'ebins=i'                => \$ebins,
	    'quatized=i'             => \$quantized,
	    'divbins=i'              => \$divbins,
	    'mbins=i'                => \$mbins,
	    'type=s'                 => \$type,

	    'shuffle_mifind=i'       => \$shuffle_mifind,
	    'maxp_mifind=f'          => \$maxp_mifind,
	    'shuffle_mimotif=i'      => \$shuffle_mimotif,
	    'maxp_mimotif=f'         => \$maxp_mimotif,
	    'shuffle_page=i'         => \$shuffle_page,
	    'maxp_page=f'            => \$maxp_page,
	    'maxz_mioptimize=f'      => \$maxz_mioptimize,
	    'maxz_mifind=f'          => \$maxz_mifind,
	    'jn_t=i'                 => \$jn_t,
      'min_r=f'                => \$min_r,
      'dG_t=f'                 => \$dG_t,

	    'cluaters=i'             => \$clusters,
	    'order=i'                => \$order,
	    'draw_min=f'             => \$draw_min,
	    'draw_max=f'             => \$draw_max,
	    'removecols_draw=i'      => \$removecols_draw,

	    'dodrawmatrix=i'         => \$dodrawmatrix,
	    'dodrawmimatrix=i'       => \$dodrawmimatrix,
	    'dodrawpagematrix=i'     => \$dodrawpagematrix,
	    'dodrawmotifs=i'         => \$dodrawmotifs,
	    'doskipdiscovery=i'      => \$doskipdiscovery,
	    'doskipoptimization=i'      => \$doskipoptimization,
	    'domifind=i'             => \$domifind,
	    'domioptimize=i'         => \$domioptimize,
	    'doonlypositive=i'       => \$doonlypositive,
	    'domireport=i'           => \$domireport,
	    'doremovedups=i'         => \$doremovedups,
	    'dopagerun=i'            => \$dopagerun,
	    'domotifind=i'           => \$domotifind,
	    'doconservation=i'       => \$doconservation,

	    'suffix=s'               => \$suffix,
	    'submit=i'               => \$submit,
	    'platform=s'             => \$platform,
	    'walltime=s'             => \$walltime,) ;

if (!defined($expfile)) {
    die("Please input an expression file (--expfile=FILE).\n");
}

if ($doskipdiscovery == 1) {
    $domifind       = 0;
    $domioptimize   = 0;
}

if ($doskipoptimization == 1) {    
    $domioptimize   = 0;
}

if ($exptype eq "discrete") {
  $quantized = 1 ;
} elsif ($exptype eq "continuous"){
  $quantized = 0 ;
}

if (($quantized == 0) && (!defined($removecols_draw))) {
    $removecols_draw = 0;
}

if (defined($species)) {
  #  read species file
    my $species_data = readSpeciesData($species);
  
    $fastafile_rna    = $species_data->{"fastafile_rna"} if (!defined($fastafile_rna)) ;
    $fastafile_ort    = $species_data->{"fastafile_ort"} if (!defined($fastafile_ort)) ;
    $homologyfile     = $species_data->{"homologyfile"} if (!defined($homologyfile)) ;
    $goindexfile      = $species_data->{"goindexfile"} if  (!defined($goindexfile));
    $gonamesfile      = $species_data->{"gonamesfile"} if  (!defined($gonamesfile));
    $seedfile         = $species_data->{"seedfile"} if  (!defined($seedfile));
}

my $time = Sets::getNiceDateTime(1) ;
my $a_ref_files = Sets::getFiles($expfile);

foreach my $expfile (@$a_ref_files) {
  if (&check_input_file($expfile, $exptype) == 0) {
    die "Please correct input expression file.\n";
  }

  my $expfile_file = Sets::filename($expfile);

  my $target_dir = "$expfile"."_TEISER"."/$suffix";
  if ($suffix eq ""){
    $target_dir = "$expfile"."_TEISER";
  }
  if (! -e "$expfile"."_TEISER") {
    mkdir "$expfile"."_TEISER" ; 
  }
  if (! -e $target_dir) {
    mkdir $target_dir; 
  }

  open OUTC, ">$target_dir/cmdline.txt"
    or print "Cannot open $target_dir/cmdline.txt";
  print OUTC "$cmdline\n";
  close OUTC;

  ######## remove duplicates ########
  my $expfile_nodups     = "$target_dir/$expfile_file" ;

  my $pbs = PBS->new;
  $pbs->setPlatform($platform) if (defined($platform));
  $pbs->setQueue($queue)       if (defined($queue));
    
  $pbs->setWallTime($walltime);
  $pbs->addCmd("export TEISERDIR=$teiserdir");
  $pbs->addCmd("cd $pwd");
    
  $pbs->setScriptName("$target_dir/$expfile_file.script");
  $pbs->addCmd("date");
    
  $pbs->addCmd("echo \"Remove duplicates, create $target_dir/$expfile_file\"");
    
  if  ($doremovedups == 1) {
    my %PARAMS = ("expfile"       => $expfile,
		  "quantized"     => $quantized, 
		  "fastafile"     => $fastafile_rna,
		  "dupfile"       => undef,
		  "ebins"         => $ebins,
		  "divbins"       => $divbins,    
		  "outfile"       => $expfile_nodups) ;
	
    my $cmd = &get_cmd_removedups(\%PARAMS) ;
    $pbs->addCmd($cmd) ;
  }else{
    system("cp $expfile $expfile_nodups") ;
  }

  ######## quatize expression profile ########
  my $expfile_quant      = "$target_dir/$expfile_file.q" ;

  if ($quantized == 0) {
    $pbs->addCmd("echo \"Quantizing the input file.\"");    
    my $cmd = "perl $scriptdir/quantize_input_vector.pl --expfile=$expfile_nodups --outfile=$expfile_quant";
    if (defined($ebins)) {
      $cmd .= " -ebins $ebins ";
    }
    if (defined($divbins)) {
      $cmd .= " -divbins $divbins ";
    }
    $pbs->addCmd($cmd);
  }
  
  ######## Seed discovery ########
  my $seed_cfg_file      = "$target_dir/$expfile_file.seed.bin" ;
  my $seed_data_file     = "$target_dir/$expfile_file.seed.dat" ;

  if ($domifind == 1){
    $pbs->addCmd("echo \"step 1: seed discovery.\""); 

    my %PARAMS = ("expfile"        => $expfile_nodups,
		  "quantized"       => $quantized,
		  "rna_fastafile"   => $fastafile_rna, 
		  "shuffle"         => $shuffle_mifind,
		  "seedfile"        => $seedfile,
		  "ebins"           => $ebins,
		  "divbins"         => $divbins,
		  "seedoutfile"     => $seed_cfg_file,
		  "dataoutfile"     => $seed_data_file,
		  "max_p"           => $maxp_mifind,
		  "max_z"           => $maxz_mifind,) ;
    my $cmd = &get_cmd_mifind(\%PARAMS);
    $pbs->addCmd($cmd);
    # if no seeds, exit here
    $pbs->addCmd("CNTM=`perl $scriptdir/count_lines.pl $seed_data_file`; if [ \$CNTM = 0 ]; then echo \"No DNA motifs. Exiting now.\"; exit; fi");
  }else{
    my $cmd = "cp $seedfile $seed_cfg_file" ;
    $pbs->addCmd($cmd);
  }
  
  ######## Seed optimization ########
  my $optim_motif_file   = "$target_dir/$expfile_file.optim.bin" ;
  my $optim_data_file    = "$target_dir/$expfile_file.summary" ;

  $pbs->addCmd("echo \"step 2: seed optimization.\""); 
  my %PARAMS = ("expfile"         => $expfile_nodups,
		"quantized"       => $quantized,
		"location"        => $type,
		"rna_fastafile"   => $fastafile_rna,
		"shuffle"         => $shuffle_mifind,
		"seedfile"        => $seed_cfg_file,
		"ebins"           => $ebins,
		"dooptimize"      => $domioptimize,
		"doonlypositive"  => $doonlypositive,
		"divbins"         => $divbins,
    "minr"            => $min_r,
		"motifoutfile"    => $optim_motif_file,
		"dataoutfile"     => $optim_data_file,
		"max_p"           => $maxp_mifind,
		"jn_t"            => $jn_t,
    "dG_t"            => $dG_t,
		"max_z"           => $maxz_mioptimize) ;

  my $cmd = &get_cmd_mioptimize(\%PARAMS);
  $pbs->addCmd($cmd);

  if ($doskipdiscovery == 1) {
    $pbs->addCmd("echo \"Step 2: skip seed optimization, use -seedfile option instead.\"");
    die "Please provide -seedfile: $seedfile\n" if (! -e $seedfile);
    #$pbs->addCmd("cp $seedfile $optim_motif_file");
  }

  ######## generating reports ########
  my $motif_profile_dir  = "$target_dir/Motifs" ;
  if (! -e $motif_profile_dir) {
    mkdir $motif_profile_dir ; 
  }

  my $reportfile         = "$target_dir/$expfile_file.report" ;
  my $matrixfile         = "$target_dir/$expfile_file.matrix" ;
  if ($domireport == 1) {
    $pbs->addCmd("echo \"step 3: creating report.\"");
    my %PARAMS = ("expfile"         => $expfile_nodups,
		  "quantized"       => $quantized,
		  "rna_fastafile"   => $fastafile_rna,
		  "motiffile"       => $optim_motif_file,
		  "ebins"           => $ebins,
		  "divbins"         => $divbins,
		  "matrixfile"      => $matrixfile,
		  "reportfile"      => $reportfile,
		  "profiledir"      => $motif_profile_dir,
      "dG_t"            => $dG_t,) ;
    
    my $cmd = &get_cmd_mireport(\%PARAMS);
    $pbs->addCmd($cmd);
  }

  ######## generating motif-motif interaction matrix ########
  my $mimatrixfile      = "$target_dir/$expfile_file.mimatrix" ;
  if ($domotifind == 1) {
    $pbs->addCmd("echo \"step 4: generating motif-motif interaction matrix.\"");
    my %PARAMS = ("expfile"         => $expfile_nodups,
		  "quantized"       => $quantized,
		  "rna_fastafile"   => $fastafile_rna,
		  "motiffile"       => $optim_motif_file,
		  "ebins"           => $ebins,
		  "divbins"         => $divbins,
		  "matrixfile"      => $mimatrixfile,
		  "pvfile"          => $matrixfile,
		  "max_p"           => $maxp_mimotif,
		  "shuffle"         => $shuffle_mimotif, 
		  "profiledir"      => $motif_profile_dir,
      "dG_t"            => $dG_t,) ;
    
    my $cmd = &get_cmd_motifind(\%PARAMS);
    $pbs->addCmd($cmd);
  }

  ######## running PAGE ########
  if ($dopagerun == 1) {
    $pbs->addCmd("echo \"step 5: running PAGE.\"");
    my %PARAMS = ("dirfile"        => $motif_profile_dir,
		  "motiffile"       => $optim_motif_file,
		  "goindexfile"     => $goindexfile,
		  "gonamesfile"     => $gonamesfile,
		  "max_p"           => $maxp_page,
		  "shuffle"         => $shuffle_page, ) ;
    
    my $cmd = &get_cmd_pagerun(\%PARAMS);
    $pbs->addCmd($cmd);
  }
  ######## calculating conservation index ########
  my $consfile     = "$target_dir/$expfile_file.cons" ;
  if ($doconservation == 1){
    $pbs->addCmd("echo \"step 6: calculating conservation score.\"");
    my %PARAMS = ("motiffile"      => $optim_motif_file,
		  "fastafile1"      => $fastafile_rna,
		  "fastafile2"      => $fastafile_ort,
		  "homologyfile"    => $homologyfile,
		  "consfile"        => $consfile, 
      "dG_t"            => $dG_t,) ;
    
    my $cmd = &get_cmd_consrun(\%PARAMS);
    $pbs->addCmd($cmd);
  }

  ######## drawing matrix ########
  if ($dodrawmatrix == 1){
    $pbs->addCmd("echo \"step 7: drawing matrix.\"");
     my %PARAMS = ("pvmatrixfile"  => $matrixfile,
		  "summaryfile"     => $optim_data_file,
		  "expfile"         => $expfile_nodups,
		  "quantized"       => $quantized,
		  "colmap"          => $colmap_matrix,
		  "order"           => $order,
		  "min"             => $draw_min,
		  "max"             => $draw_max,
		  "cluster"         => $clusters,
      "suffix"         =>  "",) ;
    
    my $cmd = &get_cmd_drawmatrix(\%PARAMS);
    $pbs->addCmd($cmd);

    $pbs->addCmd("echo \"step 7: drawing matrix.\"");
     my %PARAMS = ("pvmatrixfile"  => $matrixfile.".linear",
      "summaryfile"     => $optim_data_file,
      "expfile"         => $expfile_nodups,
      "quantized"       => $quantized,
      "colmap"          => $colmap_matrix,
      "order"           => $order,
      "min"             => $draw_min,
      "max"             => $draw_max,
      "cluster"         => $clusters,
      "suffix"         =>  "linear",) ;
    
    my $cmd = &get_cmd_drawmatrix(\%PARAMS);
    $pbs->addCmd($cmd);
  }

######## drawing html matrix ########
  if ($dodrawmatrix == 1){
      $pbs->addCmd("echo \"step 7.5: drawing html matrix.\"");
     my %PARAMS = ("pvmatrixfile"  => $matrixfile,
                  "summaryfile"     => $optim_data_file,
                  "expfile"         => $expfile_nodups,
                  "quantized"       => $quantized,
                  "colmap"          => $colmap_matrix,
                  "order"           => $order,
                  "min"             => $draw_min,
                  "max"             => $draw_max,
		   "cluster"         => $clusters,) ;

      my $cmd = &get_cmd_drawhtmlmatrix(\%PARAMS);
      $pbs->addCmd($cmd);
  }

  ######## drawing mimatrix ########
  if ($dodrawmimatrix == 1){
    $pbs->addCmd("echo \"step 8: drawing mimatrix.\"");
     my %PARAMS = ("pvmatrixfile"  => $mimatrixfile,
		  "summaryfile"     => $optim_data_file,
		  "expfile"         => $expfile_nodups,
		  "quantized"       => $quantized,
		  "colmap"          => $colmap_mimatrix,
		  "order"           => $order,
		  "min"             => $draw_min,
		  "max"             => $draw_max,
		  "cluster"         => $clusters,) ;
    
    my $cmd = &get_cmd_drawmimatrix(\%PARAMS);
    $pbs->addCmd($cmd);
  }

  ######## drawing page matrix ########
  if ($dodrawpagematrix == 1){
    $pbs->addCmd("echo \"step 9: drawing page matrix.\"");
     my %PARAMS = ("summaryfile"     => $optim_data_file,
		  "colmap"          => $colmap_page,) ;
    
    my $cmd = &get_cmd_drawpagematrix(\%PARAMS);
    $pbs->addCmd($cmd);
  }

  ######## drawing motifs ########
  if ($dodrawmotifs == 1){
    $pbs->addCmd("echo \"step 10: drawing motifs.\"");
     my %PARAMS = ("summaryfile"     => $optim_data_file) ;
    
    my $cmd = &get_cmd_drawmotifs(\%PARAMS);
    $pbs->addCmd($cmd);
  }
  
  if ($submit == 0) {
    $pbs->execute;
  }else{
    $jobid = $pbs->submit; print "Submitted job $jobid.\n";
  }
}





sub readSpeciesData {
  my ($species) = @_;
  
  my %H = ();
  open IN, "$ENV{TEISERDIR}/TEISER_Data/species_data/$species" or die "No data file for $species.\n";
  while (my $l = <IN>) {
    chomp $l;
    my @a = split /\t/, $l, -1;    
    if ($a[1] =~ /^TEISER_Data/) {
      $a[1] = "$ENV{TEISERDIR}/$a[1]";
    }
    $H{$a[0]} = $a[1];
  }  
  close IN;
  
  return \%H;
}

sub check_input_file {
  my ($expfile, $exptype) = @_;
  
  my $ta = Table->new;
  $ta->loadFile($expfile);
  my $a_ref = $ta->getArray();
  
  my $r = shift @$a_ref; 
  if ($r->[1] =~ /^\d/) {
    print "WARNING: your file might not contain a header line ($r->[1]).\n";
  }

  if (scalar(@$a_ref)<2){
    print "ERROR: your expfile is not valid. Most probably it contains '\r' instead of '\n'\n" ;
    return 1 ;
  }
  
  my %H = ();
  my %V = ();
  foreach my $r (@$a_ref) {
    if (defined($H{$r->[0]})) {
      print "Your files contains multiple rows with the same gene id. Please correct that before applying TEISER.\n";
      return 0;
    }
    $H{$r->[0]} ++;
    $V{$r->[1]} = 1;
  } 
  
  
  my @v = values ( %V );
  @v = sort { $a <=> $b } @v;
  my $max = $#v;
  if (($exptype eq 'discrete') && (scalar(@v) != $max+1)) {
    my $n1 = scalar(@v);
    my $n2 = $max + 1;
    die "Problem. Your discrete vector is missing some symbols ($n1 != $n2).\n";
    return 0;
  }
  
  return 1;
}

sub get_cmd_removedups {
  my ($p) = @_;
  my $todo = "perl $scriptdir/remove_duplicates.pl --expfile=$p->{expfile} --quantized=$p->{quantized} --fastafile=$p->{fastafile} --outfile=$p->{outfile}";
  if (defined($p->{dupfile})) {
    $todo .= " --dupfile=$p->{dupfile} ";
  }
  if (defined($p->{ebins})) {
    $todo .= " --ebins=$p->{ebins} ";
  }
  if (defined($p->{divbins})) {
    $todo .= " --divbins=$p->{divbins} ";
  }
  return $todo;
}

sub get_cmd_mifind {
  my ($p) = @_;
  my $todo = "$programdir/mi_find_seed -expfile $p->{expfile} -seedfile $p->{seedfile} -quantized $p->{quantized} -rna_fastafile $p->{rna_fastafile} -shuffle $p->{shuffle} -seedoutfile $p->{seedoutfile} -dataoutfile $p->{dataoutfile} -max_p $p->{max_p}  -max_z $p->{max_z}" ;
  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} ";
  }
  if (defined($p->{divbins})) {
    $todo .= " -divbins $p->{divbins} ";
  }
  return $todo;
}

sub get_cmd_mioptimize {
  my ($p) = @_;
  my $todo = "$programdir/mi_optimize -expfile $p->{expfile} -seedfile $p->{seedfile} -quantized $p->{quantized} -rna_fastafile $p->{rna_fastafile} -shuffle $p->{shuffle} -motifoutfile $p->{motifoutfile} -dataoutfile $p->{dataoutfile} -max_p $p->{max_p} -location $p->{location} -dooptimize $p->{dooptimize} -doonlypositive $p->{doonlypositive} -max_z $p->{max_z} -minr $p->{minr} -dG_t $p->{dG_t}" ;
  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} ";
  }
  if (defined($p->{jn_t})) {
      $todo .= " -jn_t $p->{jn_t} ";
  }
  if (defined($p->{divbins})) {
    $todo .= " -divbins $p->{divbins} ";
  }
  return $todo;
}

sub get_cmd_mireport {
  my ($p) = @_;
  my $todo = "$programdir/mi_report -expfile $p->{expfile} -motiffile $p->{motiffile} -quantized $p->{quantized} -rna_fastafile $p->{rna_fastafile} -reportfile $p->{reportfile} -matrixfile $p->{matrixfile} -profiledir $p->{profiledir} -dG_t $p->{dG_t}" ;
  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} ";
  }
  if (defined($p->{divbins})) {
    $todo .= " -divbins $p->{divbins} ";
  }
  return $todo;
}

sub get_cmd_motifind {
  my ($p) = @_;
  my $todo = "$programdir/mi_motif_motif_interaction -expfile $p->{expfile} -motiffile $p->{motiffile} -quantized $p->{quantized} -rna_fastafile $p->{rna_fastafile} -matrixfile $p->{matrixfile} -pvfile $p->{pvfile} -max_p $p->{max_p} -shuffle $p->{shuffle} -profiledir $p->{profiledir} -dG_t $p->{dG_t}" ;
  if (defined($p->{ebins})) {
    $todo .= " -ebins $p->{ebins} ";
  }
  if (defined($p->{divbins})) {
    $todo .= " -divbins $p->{divbins} ";
  }
  return $todo;
}

sub get_cmd_pagerun {
  my ($p) = @_;
  my $todo = "$programdir/motif_page_run -dirfile $p->{dirfile} -motiffile $p->{motiffile} -goindexfile $p->{goindexfile} -gonamesfile $p->{gonamesfile} -max_p $p->{max_p} -shuffle $p->{shuffle}" ;
  return $todo;
}

sub get_cmd_consrun {
  my ($p) = @_;
  my $todo = "$programdir/mi_conserve -motiffile $p->{motiffile} -fastafile1 $p->{fastafile1} -fastafile2 $p->{fastafile2} -homologyfile $p->{homologyfile} -consfile $p->{consfile} -dG_t $p->{dG_t}" ;
  return $todo;
}

sub get_cmd_drawmatrix {
  my ($p) = @_;
  my $todo = "perl $scriptdir/teiser_draw_matrix.pl --pvmatrixfile=$p->{pvmatrixfile} --summaryfile=$p->{summaryfile} --expfile=$p->{expfile} --quantized=$p->{quantized} --colmap=$p->{colmap} --order=$p->{order} --min=$p->{min} --max=$p->{max} --cluster=$p->{cluster} --suffix=$p->{suffix}" ;
  return $todo;
}

sub get_cmd_drawhtmlmatrix {
    my ($p) = @_;
  my $todo = "perl $scriptdir/teiser_draw_matrix_html.pl --pvmatrixfile=$p->{pvmatrixfile} --summaryfile=$p->{summaryfile} --expfile=$p->{expfile} --quantized=$p->{quantized} --colmap=$p->{colmap} --order=$p->{order} --min=$p->{min} --max=$p->{max} --cluster=$p->{cluster}" ;
    return $todo;
}

sub get_cmd_drawmimatrix {
  my ($p) = @_;
  my $todo = "perl $scriptdir/teiser_draw_mimatrix.pl --pvmatrixfile=$p->{pvmatrixfile} --summaryfile=$p->{summaryfile} --expfile=$p->{expfile} --quantized=$p->{quantized} --colmap=$p->{colmap}" ;
  return $todo;
}

sub get_cmd_drawpagematrix {
  my ($p) = @_;
  my $todo = "perl $scriptdir/teiser_draw_pagematrix.pl --summaryfile=$p->{summaryfile} --colmap=$p->{colmap}" ;
  return $todo;
}

sub get_cmd_drawmotifs {
  my ($p) = @_;
  my $todo = "perl $scriptdir/teiser_draw_motifs.pl --summaryfile=$p->{summaryfile}" ;
  return $todo;
}
