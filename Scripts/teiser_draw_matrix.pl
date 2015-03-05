my $teiserdir ;
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
my $suffix            = "" ;

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
	    'order=s'             => \$order,
      'suffix=s'             => \$suffix);

my $dir = Sets::dirname($expfile)."/Motifs" ;
if ( ! -d $dir ) {
    mkdir $dir ;
}

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
my $A_REF_COLMAP = undef;
if (defined($colmap)) {

 $ta->setDelim(" ");
 $ta->loadFile($colmap);
 $A_REF_COLMAP = $ta->getArray();
 $ta->setDelim('\t');
}


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
    open O, "> $summaryfile.short" or die ;
    print O "index\tlocation\tmotif-seq\tmotif_structure\tmi-value\tseq mi-value\tfrequency\tz-score\trobustness\tp-value\tRho\tSpearman P\n" ;
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
	print O join("\t", @$dat), "\t", $spearman{$ss}{rho}, "\t", $spearman{$ss}{pv}, "\n" ;
    }
    close(O) ;
}

# left and top margins
my $xbase        = 60;
my $ybase        = 100;
$ybase           = 100 if ($quantized==0) ;

# height of each entry in the matrix
$h            = 30 if ! (defined $h) ;

# size of the EPS image to be generated
my $ysize        = $ybase + $h * scalar(@$a_ref_M) + 100;
my $h_h = 10 ;

$xsize           = 1000 if ! (defined $xsize);

# width of each entry in the matrix
$w = int( 0.5 + ((1 * $xsize / 2.0) + 5)  / @$a_ref_H )-1 if (! defined($w));
$xsize += 400 ;
print "Start drawing\n";

my $p = new PostScript::Simple(xsize     => $xsize,
			                        ysize     => $ysize,
                              colour    => 1,
                              eps       => 1,
                              units     => "pt");

$p->setlinewidth(0.5);

# show header (cluster indices)
if ($quantized == 1) {
  # drawing header
  for (my $j=0; $j<@$a_ref_H; $j++) {
    $p->setcolour("black");
    $p->setfont("Arial", $w*3/4);
    if ($w*3/4 > 20){
	$p->setfont("Arial", 20);
    }
    $p->text( { rotate => 90 }, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-3), $a_ref_H->[$j]);
  }
    
}else{
  # graphical header for continuous data
  
  my $min_i1 =  1000000;
  my $max_i2 = -1000000;
  my @bins   = ();
  foreach my $c (@$a_ref_H) {
    my ($i1, $i2) = $c =~ /\[(.+?)\ (.+?)\]/;    
    $min_i1 = $i1 if ($i1 < $min_i1);
    $max_i2 = $i2 if ($i2 > $max_i2);
    my @a_tmp = ($i1, $i2); push @bins, \@a_tmp;
  }
  
  my $th = $max_i2 - $min_i1;
  my $hi = $h * 1.5;

  my $j = 0;
  foreach my $c (@$a_ref_H) {

    my $h1 = $hi * ($bins[$j]->[0] - $min_i1 ) / $th ;
    my $h2 = $hi * ($bins[$j]->[1] - $min_i1 ) / $th ;

    $p->setcolour("black");    
    $p->box({filled => 1}, 
	    $xbase + $j * $w,      $ysize - ($ybase - 5 - $hi) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ));

   
    $p->setcolour("white");
    $p->line( $xbase + $j * $w + $w,      $ysize - ($ybase - 5 - $hi) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ) );
    
    $p->setcolour("red");    
    $p->box({filled => 1}, 
	    $xbase + $j * $w,      $ysize - ($ybase - 5 - $h2) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 - $h1 ));
    
    my $cc = $c; $cc =~ s/^.+\]//; #print "$c $cc\n";
    
    $p->setcolour("black");
    $p->setfont("Arial", 6);
    $p->text( { rotate => 90 }, $xbase + $j * $w + 3*$w/4, $ysize - ($ybase-5 - $hi), $cc);  


    $j ++;
  }

  $p->setlinewidth(0.5);

  $j = 0;
   foreach my $c (@$a_ref_H) {
    
    my $h1 = $hi * ($bins[$j]->[0] - $min_i1 ) / $th ;
    my $h2 = $hi * ($bins[$j]->[1] - $min_i1 ) / $th ;
   
    $p->setcolour("white");
    $p->line( $xbase + $j * $w + $w,      $ysize - ($ybase - 5 - $hi) , 
	    $xbase + $j * $w + $w, $ysize - ($ybase - 5 + 0 ) );
    
    $j ++;
  }

  

  print "$max_i2\t$min_i1\n";

  $p->setfont("Arial", 8);

  $p->setcolour("black");    
  $p->text( {align => 'right'}, $xbase - 1 , $ysize - ($ybase - 5 - $hi + 4), $max_i2);
  $p->text( {align => 'right'}, $xbase - 1, $ysize - ($ybase - 5 - 0   + 1), $min_i1);

}



$p->setcolour("black");
$p->setfont("Arial", 8);


#
# draw (i,j) p-value matrix itself
#
# set min max p-values
$max =  $minmax_lp if !(defined $max);
$min = -$minmax_lp if !(defined $min);
my @col = ();
my @ss ;

$p->setfont("Arial", 9);
$p->setcolour("black");
my @head = ("id", "location", "motif", "structure", "MI (bits)", "MI (seq)", "frequency", "z-score", "robustness") ;
push(@head, "conservation") if (defined $h_cons) ;

$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 15 , $ysize - ($ybase-$h/4), @head[0]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 45 , $ysize - ($ybase-$h/4), @head[1]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 105 , $ysize - ($ybase-$h/4), @head[2]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 235 , $ysize - ($ybase-$h/4), @head[3]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 355 , $ysize - ($ybase-$h/4), @head[4]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 400 , $ysize - ($ybase-$h/4), @head[5]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 450 , $ysize - ($ybase-$h/4), @head[6]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 495 , $ysize - ($ybase-$h/4), @head[7]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 540 , $ysize - ($ybase-$h/4), @head[8]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 595 , $ysize - ($ybase-$h/4), @head[9]) if (defined $h_cons) ;

for (my $i=0; $i<@$a_ref_M; $i++) {
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

   # draw the matrix entry with appropriate color
   $p->setcolour(@col);
   $p->box({filled => 1},
           $xbase + $j * $w,      $ysize - ($ybase + $i*$h) ,
           $xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h)));

 }
 $p->setcolour("black");
 $p->setfont("Arial", 10);
 my $dat = $hsummary->{$ss} ;

 ### draw the motif
 my $myii = $dat->[0] ;
 my $mo = $dat->[2] ;
 my $mo = &seq_to_re($mo);
 my $mo = Sets::myre2wm($mo);

 open OUT, ">$dir/$myii\_motif.txt" or die "cannot open $dir/$myii\_motif.txt\n";
 print OUT $mo;
 close OUT;
    
 system("$scriptdir/weblogo/seqlogo -f $dir/$myii\_motif.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $dir/$myii\_motif.eps");

 my $e  = new PostScript::Simple::EPS(file => "$dir/$myii\_motif.eps");

 # get height
 my $eh = $e->height;
 my $ew = $e->width;
 # height must be $h, so scale down to $h = k * LO
 $e->scale($h / $eh);
 my $ew_new = int(0.5 + $ew * $h / $eh);
 
 $p->_add_eps($e, $xbase +(@$a_ref_H) * $w + 105 - $ew_new / 2,  $ysize - ($ybase + ($i*$h+$h))); 

 $p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 15 , $ysize - ($ybase + ($i*$h+$h/2)+2), $dat->[0]) ;
 $p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 45 , $ysize - ($ybase + ($i*$h+$h/2)+2), $dat->[1]) ;
 $p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 235 , $ysize - ($ybase + ($i*$h+$h/2)+2), $dat->[3]) ;
 $p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 355 , $ysize - ($ybase + ($i*$h+$h/2)+2), sprintf("%2.3f", $dat->[4])) ;
 $p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 400 , $ysize - ($ybase + ($i*$h+$h/2)+2), sprintf("%2.3f", $dat->[5])) ;
 $p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 450 , $ysize - ($ybase + ($i*$h+$h/2)+2), $dat->[6]) ;
 $p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 495 , $ysize - ($ybase + ($i*$h+$h/2)+2), $dat->[7]) ;
 $p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 540 , $ysize - ($ybase + ($i*$h+$h/2)+2), $dat->[8]) ;
 $p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 595 , $ysize - ($ybase + ($i*$h+$h/2)+2), $h_cons->{$dat->[1].$dat->[0]}) if (defined $h_cons);
 
}

####### draw significant boxes ##########
my $NN   = @$a_ref_H ;
for (my $i=0; $i<@$a_ref_M; $i++) {
  my $r  = $a_ref_M->[$i];
  my $w1 = 50 ;
  # ss through entries in current row
  for (my $j=0; $j<@$r; $j++) {
    
    # get log pvalues for over and under rep
    my $lp = $r->[$j] ;
    my $tp = - Sets::log10( 0.05 / $NN ) ;
    
    # create appropriate color
    my $col = undef;
    if (abs($lp) > $tp) {
      my $col = undef;
      if ($lp > 0) {
	$col = "red";
      }else{
	$col = "blue";
      }
      # draw the matrix entry with appropriate color
      $p->setcolour($col);
      $p->box({filled => 0},
	      $xbase + $j * $w,      $ysize - ($ybase + $i*$h) ,
	      $xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h)));
    }
  }
}
#
# draw scale bar
#
drawScale($xscale, $ysize / 2 - $yscale , $min, $max, 50, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, "Over-representation", "Under-representation", $scalefont, 2, 20);

my $outeps = "$expfile"."$suffix.summary.eps" ;
my $outpdf = "$expfile"."$suffix.summary.pdf" ;
my $ps2pdf = 1;


# output EPS file
print "Outputing EPS file $outeps\n";
$p->output("$outeps");

# convert to PDF
print "Convert to PDF $outpdf\n";
print("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf\n");
system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");

print "Finished.\n";
exit(0);


sub drawScale {
 my ($x, $y, $min, $max, $res, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, $upper_text, $lower_text, $scalefont, $h, $w) = @_;

 my $sep = 0;

 $p->setcolour("black");
 $p->setfont("Arial", 10);
 $p->text({align => "center"}, $x+$w/2+0, $ysize - ($y - 3), $max);

 $p->setfont("Times", $scalefont);

 $p->text({align => "left", rotate => 90}, $x+$w/2+5, $ysize - ($y - 16), $upper_text);


 my $t = $max;


 for (my $i=0; $i<=$res; $i++) {

   my @col = () ;
   if (!defined($colmap)) {
       if ($i>$res/2)
       {
	   @col = Sets::interp_general( $t, [0, 0, 0], [255, 0, 0], $min, 0);
       }
       else
       {
	   @col = Sets::interp_general( $t, [255, 0, 0], [255, 255, 0], 0, $max);
       }
   } else {
     @col = Sets::interp_from_matlab_colormap( $t, $A_REF_COLMAP, $min, $max);
   }

   $p->setcolour( @col );
   $p->box({filled => 1}, $x, $ysize - ($y + $sep + $i*$h) , $x+$w, $ysize - ($y + $sep + $i*$h + $h));
   $t -= ($max - $min) / $res;
 }

 $p->setcolour( 'black' );
 $p->setfont("Arial", 10);
 $p->text({align => "center"}, $x+$w/2+1, $ysize - ($y + $sep + $res*$h + 11), $min);
 $p->setfont("Times", $scalefont);
 $p->text({align => "right", rotate => 90}, $x+$w/2+5, $ysize - ($y + $sep + $res*$h + 11 + 10), $lower_text);


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
