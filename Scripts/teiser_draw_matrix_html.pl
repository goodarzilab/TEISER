my $disredir ;
BEGIN{
    if ((!$ENV{TEISERDIR}) || ($ENV{TEISERDIR} eq '')) {
	$teiserdir="./" ;
	print "The TEISERDIR environment variable is not set. It is set to default.\n";
    }
    else{
	$teiserdir = $ENV{TEISERDIR};
    }
}

use lib "$teiserdir/Scripts";
use lib "$teiserdir/Scripts/PostScript-Simple-0.07/lib";

my $programdir = $teiserdir."/Programs" ;
my $scriptdir  = $teiserdir."/Scripts" ;

use Table;
use Sets;
use Stat;
use Getopt::Long;
use PostScript::Simple;
use AggloClust;
use strict;
use Data::Dumper ;


my $pvmatrixfile      = undef;
my $pvmatrixfile_up   = undef;
my $pvmatrixfile_dn   = undef;
my $consfile_up       = undef;
my $consfile_dn       = undef;
my $summaryfile       = undef;
my $expfile           = undef;
my $datafile          = undef;
my $colmap            = "$scriptdir/HEATMAPS/cmap_1.txt";
my $cluster           = undef;
my $max_p             = undef;
my $minmax_lp         = 3;
my $quantized         = 1;
my $min               = -10 ;
my $max               = 10  ;
my $xsize             = undef;
my $xscale            = 5 ;
my $yscale            = 75 ;
my $scalefont         = 9 ;
my $h                 = undef ;
my $w                 = undef ;
my $draw_sample_heatmap = "false" ;
my $order             = 0 ;

if (@ARGV == 0) {
 die "Usage: perl teiser_draw_matrix.pl  --pvaluematrixfile=FILE --expfile=FILE --max_p=P\n";
}

GetOptions ('pvmatrixfile=s'      => \$pvmatrixfile,
	    'pvmatrixfile_up=s'   => \$pvmatrixfile_up,
	    'pvmatrixfile_dn=s'   => \$pvmatrixfile_dn,
	    'consfile_up=s'       => \$consfile_up,
	    'consfile_dn=s'       => \$consfile_dn,
	    'summaryfile=s'       => \$summaryfile,
            'cluster=s'           => \$cluster,
	    'expfile=s'           => \$expfile,
	    'quantized=s'         => \$quantized,
	    'colmap=s'            => \$colmap,
	    'minmax_lp=s'         => \$minmax_lp,
	    'min=s'               => \$min,
	    'max=s'               => \$max,
	    'xsize=s'             => \$xsize,
	    'xscale=s'            => \$xscale,
	    'yscale=s'            => \$yscale,
	    'scalefont=s'         => \$scalefont,
	    'h=s'                 => \$h,
	    'w=s'                 => \$w,
	    'max_p=s'             => \$max_p,
	    'order=s'             => \$order);

my $outfile = "$expfile".".summary.html" ;
my $dir = Sets::dirname($expfile)."/Motifs" ;
my $fn  = Sets::filename($summaryfile);
my $dname = "/Motifs";

my $ta = Table->new;
my $a_ref_M ;
my $a_ref_H ;
my $h_cons = undef ;
print "Reading MI data ... ";

my $count_rep ;
my $single_source ;
#  read in the matrix file
if (defined $pvmatrixfile){
  $ta->loadFile($pvmatrixfile);
  $a_ref_M      = $ta->getArray();
  $a_ref_H      = shift @$a_ref_M; shift @$a_ref_H;
  $count_rep    = scalar(@$a_ref_M) ;
  $single_source=1 ;
}else{
  $ta->loadFile($pvmatrixfile_up);
  $a_ref_M      = $ta->getArray() ;
  $a_ref_H      = shift @$a_ref_M; shift @$a_ref_H;
  $count_rep    = scalar(@$a_ref_M) ;
  $ta->loadFile($pvmatrixfile_dn);
  my $a_m       = $ta->getArray() ;
  shift @$a_m ;
  foreach my $r (@$a_m){
    $r->[0] += $count_rep ;
  }
  push(@$a_ref_M, @$a_m) ;
  $single_source=0 ;

  open I, "< $consfile_up" or die "$consfile_up";
  while(my $l = <I>){
    $l =~ s/\s+$// ;
    my @C = split(/\t/, $l) ;
    $h_cons->{"UP".$C[0]} = 1-$C[-2] ;
  }
  close I ;
  open I, "< $consfile_dn" or die ;
  while(my $l = <I>){
    $l =~ s/\s+$// ;
    my @C = split(/\t/, $l) ;
    $h_cons->{"DN".$C[0]} = 1-$C[-2] ;
  }
  close I ;
}

# get an 2D array
if (!defined($max_p)) {
  $max_p = 0.05 / @$a_ref_H;
}
print "Done.\n";

if (defined($cluster) && $cluster>0 && (@$a_ref_M > 5)) {

 print "Cluster rows .. ";

 my $ac = AggloClust->new;

 my @dist = ();
 my $n    = @$a_ref_M;
 for (my $i=0; $i<$n-1; $i++) {
   $dist[$i][$i] = 0;
   for (my $j=$i+1; $j<$n; $j++) {
     my @a1 = @{ $a_ref_M->[$i] }; shift @a1;
     my @a2 = @{ $a_ref_M->[$j] }; shift @a2;
     $dist[$i][$j] = 1 - Sets::pearson(\@a1, \@a2);
     $dist[$j][$i] = $dist[$i][$j];
   }
 }

 $ac->setDistanceMatrix(\@dist);
 $ac->setMaxNbClusters($cluster);
 my $a_ref_c = $ac->agglomerate_using_avg_linkage();

 my @NEWMAT = ();
 foreach my $c (@$a_ref_c) {
   print join(" ", @$c); print "\n";
   foreach my $i (@$c) {
     push @NEWMAT, $a_ref_M->[$i];
   }
 }
 $a_ref_M = \@NEWMAT;

 print "Done.";
}

#
# load color map
#

# set min max p-values
$max =  $minmax_lp if !(defined $max);
$min = -$minmax_lp if !(defined $min);

my $xsize = 600 ;
my $xbase = 120 ;
my $ybase = 180 ;
my $w = int( 0.5 + $xsize / @$a_ref_H ) ;
if ($w<6){
  $w = 6 ;
  $xsize = $w*@$a_ref_H + 1200 ;
}else{
  $xsize = $w*@$a_ref_H + 1200 ;
}

my $h = 50 ;
my $fsize= $h/4 ;

my $A_REF_COLMAP = undef;
if (defined($colmap)) {

 $ta->setDelim(" ");
 $ta->loadFile($colmap);
 $A_REF_COLMAP = $ta->getArray();
 $ta->setDelim('\t');
}

open O, "> $outfile";
print O (
	"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\"
   \"http://www.w3.org/TR/html4/strict.dtd\">\n"
);
print O ("<head>\n");
print O ("<title>TEISER results</title>\n");
print O ("<style type=\"text/css\">");
print O ("
* {
padding: 0;
margin: 0;
}
body {
font-family: Helvetica, Arial, sans-serif;
}
h1, p, p.head {
width:960px;
padding: 10px;
}
td.grid {
width: $w\px;
height: $h\px;
font-size:$fsize\px;
}
td.gridA {
font-size:$fsize\px;
font-family: Helvetica, Arial, sans-serif;
text-align: center;
padding: 0px 2px;
}
td.gridheader {
height: auto;
width: $w\px;
border: 1px solid #FFF;
text-align: center;
background-color: #000;
}
td.gridheaderA {
height: auto;
width: $w\px;
border: 1px solid #FFF;
text-align: center;
font-size: 10px;
}
td.gridheaderB {
height: auto;
width: $w\px;
border: 1px solid #FFF;
text-align: center;
background-color: #FFF;
}
td.gridlabel {
white-space: nowrap;
height: $w\px;
padding-left: 5px;
}
table.legend {
margin-left: 10px;
top: 0px;
left: 0px;
position: relative;
}
table.legend td {
text-align: center;
}
table.legend tr {
font-size:0px;
height: 2px;
}
table.main {
position: absolute;
top: $ybase\px;
left: $xbase\px;
padding: 10px;
}
");
print O ("</style>\n");
print O ("</head>\n");
print O ("<body>\n\n");

print O ("<h1><abbr title=\"Tool for Eliciting Informative Structural Elements in RNA\">TEISER</abbr> results</h1>\n");

print O ("<p class=\"head\">These results are also available as <a href=\"./$fn.pdf\"><abbr title=\"Portable Document Format\">PDF</abbr></a> and <a href=\"./$fn.eps\"><abbr title=\"Encapsulated PostScript\">EPS</abbr></a> documents.</p>\n");

print O ("<p>Depending on your display resolution, scrolling or zooming may be necessary.</p>\n");



############ Scale ############
my $scalex = $xbase/3 ;
my $hscale = $h*2/3 ;
my $W = $hscale ;

print O ("<table class=\"legend\" width=\"$hscale\px\" cellpadding=\"0\" cellspacing=\"0\">\n");
my $p = new PostScript::Simple(xsize => $W/2,
			       ysize => $h+10,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Helvetica", 10);
$p->text( { rotate => 90 }, 10, 5, "over-rep") ;
$p->output("$dir/over.eps") ;
system("convert -density 100 $dir/over.eps $dir/over.png");
print O "<tr><td style=\"border:1px #FFFFFF solid;\">" ;
print O "<img src=\"./$dname/over.png\" alt=\"over-representation\" /></td></tr>\n" ;

print O "<tr><td style=\"font-size: 12px; height: 14px; \">$max</td></tr>\n" ;

my $t = $max;
my $res = 70 ;
for (my $i=0; $i<=$res; $i++) {
  my @col = () ;
  if (!defined($colmap)) {
    if ($i>$res/2){
      @col = Sets::interp_general( $t, [0, 0, 0], [255, 0, 0], $min, 0);
    }else{
      @col = Sets::interp_general( $t, [255, 0, 0], [255, 255, 0], 0, $max);
    }
  }else{
    @col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
  }
  if ($i==$res/2){
    print O "<tr>\n" ;
    print O "<td style=\"font-size: 12px; height: 14px; \">0</td>\n" ;
    print O "</tr>\n" ;
  }else{
    my $color = &RGBDecToHex($col[0], $col[1], $col[2]);
    print O "<tr>\n" ;
    print O "<td style=\"background-color:$color\">&nbsp;</td>\n" ;
    print O "</tr>\n" ;
    $t -= ($max - $min) / $res;
  }
}

print O "<tr><td style=\"font-size: 12px; height: 14px; \">$min</td></tr>\n" ;
my $p = new PostScript::Simple(xsize => $W/2,
			       ysize => $h+10,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 10);
$p->text( { rotate => 90 }, 10, 5, "under-rep") ;
$p->output("$dir/under.eps") ;
system("convert -density 100 $dir/under.eps $dir/under.png");
print O "<tr><td>" ;
print O "<img src=\"./$dname/under.png\" alt=\"under-representation\"/></td></tr>\n" ;

print O ("\n</table>\n");


#
#  START DRAWING
#
#

$ta->loadFile($summaryfile) ;
my $midata = $ta->getArray() ;
shift(@$midata) ;
my $hsummary ;
foreach my $r (@$midata){
  if ($r->[1] eq "DN" and $single_source==0){
    $hsummary->{$count_rep+$r->[0]} = $r ;
  }else{
    $hsummary->{$r->[0]} = $r ;
  }
}


if ($order ==1){
    my @NEWMAT = ();
    my @enr ;
    my @rank = (1..@$a_ref_H) ;
    my $cnt=0 ;
    open O1, "> $summaryfile.short" or die ;
    print O1 "index\tlocation\tmotif-seq\tmotif_structure\tmi-value\tseq mi-value\tfrequency\tz-score\trobustness\tp-value\tRho\tSpearman P\n" ;
    my %spearman ;
    for (my $i=0 ; $i<@$a_ref_M ; $i++){
	my $r = $a_ref_M->[$i] ;
	my $ss = $r->[0] ;
	my @v ;
	for (my $j=1 ; $j<@$r ; $j++){
	    push(@v, $r->[$j]) ;
	}
	my ($rho,$z,$pv) =Stat::spearman(\@v, \@rank) ;

	next if ($pv > 0.01) ;

	my $min = 0 ;
	my $pos = 1 ;
	for (my $j=1 ; $j<@$r ; $j++){
	    my $lpo = -1*$r->[$j] ;
	    if ($lpo<$min){
		$min = $lpo ;
		$pos = $j ;
	    }
	}
	$enr[$cnt]->{ind} = $i ;
	$enr[$cnt]->{pos} = $pos ;
	$spearman{$ss}->{rho} = sprintf("%.4f",$rho) ;
	$spearman{$ss}->{pv} = sprintf("%.2e",$pv) ;

	$cnt++ ;
    }
    my @enr = sort {$b->{pos} <=> $a->{pos}} (@enr) ;
    
    foreach my $e (@enr){
	my $i = $e->{ind} ;
	push @NEWMAT, $a_ref_M->[$i];
    }
    $a_ref_M = \@NEWMAT ;
    print "Ordered...\n" ;
    for (my $i=0 ; $i<@$a_ref_M ; $i++){
        my $ss = $a_ref_M->[$i][0] ;
	my $dat = $hsummary->{$ss} ;
	print O1 join("\t", @$dat), "\t", $spearman{$ss}{rho}, "\t", $spearman{$ss}{pv}, "\n" ;
    }
    close(O1) ;
}

my $min_i1 =  1000000;
my $max_i2 = -1000000;
my @bins   = ();
if ($quantized == 0){
  foreach my $c (@$a_ref_H) {
    my ($i1, $i2) = $c =~ /\[(.+?)\s(.+?)\]/;
    $min_i1 = $i1 if ($i1 < $min_i1);
    $max_i2 = $i2 if ($i2 > $max_i2);
    my @a_tmp = ($i1, $i2); push @bins, \@a_tmp;
  }
  my $l = $xbase-30 ;
  my $t = $ybase+10;
  $l = sprintf("%2.2f", $l) ;
  print O "<div style=\"font-size:9px ; top:$t\px ; left:$l\px ; position:absolute ;\">$max_i2</div>" ;
  my $t = $ybase+$h*1.5-8 ;
  $l = sprintf("%2.2f", $l) ;
  print O "<div style=\"font-size:9px ; top:$t\px ; left:$l\px ; position:absolute ;\">$min_i1</div>" ;
}

print O ("<table width=\"$xsize\px\" cellpadding=\"0\" cellspacing=\"0\" class=\"main\">\n");

if ($quantized == 1){
  print O "<tr>\n" ;
  for (my $j=0; $j<@$a_ref_H; $j++) {
    print O "<td class=\"gridheaderA\" style=\"background-color: #FFFFFF;\">" ;

    my $c = $a_ref_H->[$j] ;
    my $W = 20 ;

    my $p = new PostScript::Simple(xsize => $w,
				   ysize => $W+10,
				   colour    => 1,
				   eps       => 1,
				   units     => "pt");
    $p->setcolour("black");
    $p->setfont("Arial", 6);
    $p->text( { rotate => 90 }, $w*2/3, 3, $c) ;
    $p->output("$dir/$c.eps") ;
    system("convert -density 100 $dir/$c.eps $dir/$c.png");
    print O "<img width=\"$w\" src=\"./$dname/$c.png\" alt=\"$c\" /></td>\n" ;
  }
}else{
  my $th = $max_i2 - $min_i1;
	$th = 1e-3 if ($th < 1e-3) ;

  my $H = $h*1.5 ;
  print O "<tr style=\"height:$H\px\">\n" ;
  $H-=3 ;
  my $j = 0;
  foreach my $c (@$a_ref_H) {
    my $h1 = $H * ($bins[$j]->[0] - $min_i1 ) / $th ;
    my $h2 = $H * ($bins[$j]->[1] - $min_i1 ) / $th ;

    my $s = "lower bound = ".$bins[$j]->[0]." and upper bound = ".$bins[$j]->[1] ;

    print O "<td class=\"gridheader\" onclick=\"alert('$s')\" title=\"$s\">" ;
    my $dw = $w ;
    my $dh = $h2-$h1 ;
    $dh = 1 if ($dh<1) ;
    my $y = ($h1+$h2)/2-($H-2)/2 ;
    print "$y\t$H\t$h1\n" ;
#    print O "<div style=\"width:$dw\px; height:$dh\px; background-color:#F00 ; left:0\px ; bottom:$y\px; position:relative\">" ;
    print O "<div style=\"width:100\%; height:$dh\px; background-color:#F00 ; left:0\px ; bottom:$y\px; position:relative\">" ;

    print O "</div>\n" ;
    print O "</td>\n" ;
    $j++ ;
  }
}

my $xsize = 600 ;
my $xbase = 120 ;
my $ybase = 180 ;
my $w = int( 0.5 + $xsize / @$a_ref_H ) ;
if ($w<6){
  $w = 6 ;
  $xsize = $w*@$a_ref_H + 1200 ;
}else{
  $xsize = $w*@$a_ref_H + 1200 ;
}

my $h = 50 ;


#
# draw (i,j) p-value matrix itself
#
my @col = ();
my @ss ;

$p->setfont("Arial", 9);
$p->setcolour("black");
my @head = ("id", "location", "motif", "structure", "mi", "frequency", "z-score", "robustness") ;
push(@head, "conservation") if (defined $h_cons) ;

my $W = 30 ;

my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 10);
$p->text( { rotate => 45 }, 10, 5, "id") ;
$p->output("$dir/id.eps") ;
system("convert -density 100 $dir/id.eps $dir/id.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/id.png\" style=\"top:5px ; left:25px; position:relative;\" alt=\"id\"></td>\n" ;

my $W = 50 ;

my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 10);
$p->text( { rotate => 45 }, 10, 5, "location") ;
$p->output("$dir/location.eps") ;
system("convert -density 100 $dir/location.eps $dir/location.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/location.png\" style=\"top:5px ; left:25px; position:relative;\" alt=\"location\"></td>\n" ;


my $W = 80 ;

my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 10);
$p->text( { rotate => 45 }, 10, 5, "motif") ;
$p->output("$dir/motif.eps") ;
system("convert -density 100 $dir/motif.eps $dir/motif.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/motif.png\" style=\"top:5px ; left:25px; position:relative;\" alt=\"motif\"></td>\n" ;

my $W = 80 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 10);
$p->text( { rotate => 45 }, 10, 5, "structure") ;
$p->output("$dir/structure.eps") ;
system("convert -density 100 $dir/structure.eps $dir/structure.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/structure.png\" style=\"top:5px ; left:5px; position:relative;\" alt=\"structure\"></td>\n" ;

my $W = 80 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 10);
$p->text( { rotate => 45 }, 10, 5, "MI") ;
$p->output("$dir/mi.eps") ;
system("convert -density 100 $dir/mi.eps $dir/mi.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/mi.png\" style=\"top:5px ; left:30px; position:relative;\" alt=\"mutual information\"></td>\n" ;

my $W = 80 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 10);
$p->text( { rotate => 45 }, 10, 5, "frequency") ;
$p->output("$dir/frequency.eps") ;
system("convert -density 100 $dir/frequency.eps $dir/frequency.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/frequency.png\" style=\"top:5px ; left:15px; position:relative;\" alt=\"frequency\"></td>\n" ;

my $W = 80 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 10);
$p->text( { rotate => 45 }, 10, 5, "z-score") ;
$p->output("$dir/z.eps") ;
system("convert -density 100 $dir/z.eps $dir/z.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/z.png\" style=\"top:5px ; left:15px; position:relative;\" alt=\"z-score\"></td>\n" ;

my $W = 80 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 10);
$p->text( { rotate => 45 }, 10, 0, "robustness") ;
$p->output("$dir/robustness.eps") ;
system("convert -density 100 $dir/robustness.eps $dir/robustness.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/robustness.png\" style=\"top:5px ; left:15px; position:relative;\" alt=\"robustness\"></td>\n" ;

my $W = 60 ;
my $p = new PostScript::Simple(xsize => $W,
			       ysize => $h,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setcolour("black");
$p->setfont("Arial", 9);
$p->text( { rotate => 45 }, 10, 0, "conservation") ;
$p->text( { rotate => 45 }, 38, 6, "index") ;
$p->output("$dir/cons.eps") ;
system("convert -density 100 $dir/cons.eps $dir/cons.png");
print O "<td class=\"gridheaderB\">" ;
print O "<img src=\"./$dname/cons.png\" style=\"top:5px ; left:10px; position:relative;\" alt=\"conservation index\"></td>\n" ;

print O "</tr>\n" ;


for (my $i=0; $i<@$a_ref_M; $i++) {
  print O "<tr class=\"motif\">\n" ;
  # get row
  my $r  = $a_ref_M->[$i];

 # get SS description
 my $ss = shift @$r;

 my $w1 = 50 ;

 $p->setfont("ArialBold", 12);
 $p->setcolour("black");

 print "$ss\n";
 push(@ss, $ss) ;

 # ss through entries in current row
 for (my $j=0; $j<@$r; $j++) {

   # get log pvalues for over and under rep
   my $lp = $r->[$j] ;

   # create appropriate color
   my @col = ();  #$colmap = undef;
   if (!defined($colmap)) {
       if ($lp<0) {
	 @col = Sets::interp_general( $lp, [0, 255, 0], [0, 0, 0], $min, 0);
       }else{
	 @col = Sets::interp_general( $lp, [0, 0, 0], [255, 0, 0], 0, $max);
       }
   }else{
     @col = Sets::interp_from_matlab_colormap( $lp, $A_REF_COLMAP, $min, $max);
   }
   my $color = &RGBDecToHex($col[0], $col[1], $col[2]);
   my $pv ;
   my $tv = - Sets::log10( 0.05 / scalar(@$a_ref_H) );
   if ($lp>0){
     $pv = 10**(-1*$lp) ;
   }else{
     $pv = 10**$lp ;
   }
   $pv = sprintf("p = %1.4e", $pv) ;
   
   if (abs($lp) < $tv) {
     print O "<td class=\"grid\" style=\"background-color:$color\" onclick=\"alert('$pv')\" title=\"$pv\">&nbsp;</td>\n" ;
   }else{
     if ($lp>0){
       print O "<td class=\"grid\" style=\"background-color: $color; border:1px #FF0000 solid\" onclick=\"alert('over-representation: $pv')\" title=\"over-representation: $pv\">&nbsp;</td>\n" ;
     }else{
       print O "<td class=\"grid\" style=\"background-color: $color; border:1px #0000FF solid\" onclick=\"alert('under-representation: $pv')\" title=\"under-representation: $pv\">&nbsp;</td>\n" ;
     }
   }
 }
  
 ### draw the motif
 my $dat = $hsummary->{$ss} ;
 my $myii = $dat->[0] ;
 my $mo = $dat->[2] ;
 my $myre = &seq_to_re($mo);
 my $mo = Sets::myre2wm($myre);

 open OUT, ">$dir/$myii\_motif.txt" or die "cannot open $dir/$myii\_motif.txt\n";
 print OUT $mo;
 close OUT;
 
  print O "<td class=\"gridA\">" ;
  print O $dat->[0]."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O $dat->[1]."</td>\n" ;
  my $tmp_pos = $dat->[1] ;

  system("$scriptdir/weblogo/seqlogo -f $dir/$myii\_motif.txt -F PNG  -a -c -M -n -Y -w 5 -h 3 > $dir/$myii\_motif_$tmp_pos.png");
  print O "<td class=\"gridA\">" ;
  print O "<img src=\"./$dname/$myii\_motif_$tmp_pos.png\" height=\"$h\" style=\"top:5px ; position:relative; border: 0px;\" alt=\"$myre\" title=\"$myre\"></td>\n" ;

  print O "<td class=\"gridA\">" ;
  print O $dat->[3]."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O sprintf("%4.3f", $dat->[4])."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O sprintf("%4.3f", $dat->[6])."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O sprintf("%4.3f", $dat->[7])."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O $dat->[8]."</td>\n" ;
  print O "<td class=\"gridA\">" ;
  print O (defined($h_cons)?$h_cons->{$dat->[1].$dat->[0]}:'-')."</td>\n" ;
}
print O ("\n</table>\n");
print O ("\n</body>\n");
print O ("</html>\n");

exit(0);


sub RGBDecToHex {
  my ($red, $green, $blue) = @_;
  return sprintf("#%02lx%02lx%02lx", $red, $green, $blue);
}

sub interp {
 my ($r, $min, $max) = @_;

 if ($r < $min) {
   $r = $min;
 }

 if ($r > $max) {
   $r = $max;
 }

 #  1 --> red           => "0.8  0    0",
 #  0 --> blue          => "0    0    0.8",

 #  0.8 = $max . a + b
 #  0   = $min . a + b

 my $a1 = 0.8 / ($max - $min);
 my $b1 = - $a1 * $min;

 #  0   = $max . a + b
 #  0.8 = $min . a + b

 my $a3 = -0.8 / ($max - $min);
 my $b3 = -$a3 * $max;

 my $c1 = int( 0.5 + 256 * ($a1 * $r + $b1) );
 my $c2 = 0;
 my $c3 = int( 0.5 + 256 * ($a3 * $r + $b3) );


 return ($c1, $c2, $c3);

}

sub seq_to_re {
  my $S = shift @_ ;
  $S =~ s/N/./g ;
  $S =~ s/Y/[UC]/g ;
  $S =~ s/R/[AG]/g ;
  $S =~ s/K/[UG]/g ;
  $S =~ s/M/[AC]/g ;
  $S =~ s/S/[GC]/g ;
  $S =~ s/W/[AU]/g ;
  $S =~ s/B/[GUC]/g ;
  $S =~ s/D/[GAU]/g ;
  $S =~ s/H/[ACU]/g ;
  $S =~ s/V/[GCA]/g ;

  return $S ;
}
