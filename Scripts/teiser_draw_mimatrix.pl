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
use Getopt::Long;
use PostScript::Simple;
use AggloClust;
use strict;
use Data::Dumper ;


my $pvmatrixfile      = undef;
my $summaryfile       = undef;
my $expfile           = undef;
my $datafile          = undef;
my $colmap            = "$scriptdir/HEATMAPS/cmap_2.txt";
my $cluster           = undef;
my $max_p             = undef;
my $minmax_lp         = 3;
my $quantized         = 1;
my $min               = -8 ;
my $max               = 8  ;
my $xsize             = undef;
my $xscale            = 5 ;
my $yscale            = 75 ;
my $scalefont         = 8 ;
my $h                 = undef ;
my $w                 = undef ;
my $draw_sample_heatmap = "false" ;
my $order             = 0 ;

if (@ARGV == 0) {
 die "Usage: perl teiser_draw_matrix.pl  --pvaluematrixfile=FILE --expfile=FILE --max_p=P\n";
}

GetOptions ('pvmatrixfile=s'      => \$pvmatrixfile,
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

my $dir = Sets::dirname($expfile)."/Motifs" ;
my $ta = Table->new;

print "Reading MI data ... ";

#  read in the matrix file
$ta->loadFile($pvmatrixfile);

# get an 2D array
my $a_ref_M      = $ta->getArray();

# header
my $a_ref_H      = shift @$a_ref_M; shift @$a_ref_H;
if (!defined($max_p)) {
  $max_p = 0.05 / @$a_ref_H;
}
print "Done.\n";

# load color map
my $A_REF_COLMAP = undef;
if (defined($colmap)) {

 $ta->setDelim(" ");
 $ta->loadFile($colmap);
 $A_REF_COLMAP = $ta->getArray();
 $ta->setDelim('\t');
}

# draw (i,j) p-value matrix itself
$ta->loadFile($summaryfile) ;
my $midata = $ta->getArray() ;
shift(@$midata) ;

#
#  START DRAWING
#
#

# left and top margins
my $xbase        = 60;
my $ybase        = 100;
$ybase           = 100 if ($quantized==0) ;
$xsize           = 1000 if ! (defined $xsize);


# width of each entry in the matrix
$w = int( 0.5 + ((1 * $xsize / 2.0) + 5)  / @$a_ref_H )-1 if (! defined($w));
$h = $w ;
my $ysize        = $ybase + $h * scalar(@$a_ref_M) + 100;

$xsize += 100 ;
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
    $p->setfont("Arial", 9);
    my $dat = $midata->[$a_ref_H->[$j]] ;
    $p->text({align=>"left", rotate => 90}, $xbase + $j * $w + $w*3/4, $ysize - ($ybase-8), $dat->[0]." ".$dat->[1]);
  }
    
}
$p->setcolour("black");
$p->setfont("Arial", 9);

# set min max p-values
$max =  $minmax_lp if !(defined $max);
$min = -$minmax_lp if !(defined $min);
my @col = ();
my @ss ;

$p->setfont("Arial", 8);
$p->setcolour("black");
my @head = ("id", "location", "motif", "structure") ;

$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 15 , $ysize - ($ybase-$h/4), @head[0]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 45 , $ysize - ($ybase-$h/4), @head[1]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 105 , $ysize - ($ybase-$h/4), @head[2]) ;
$p->text ({align=>"left", rotate=>45}, $xbase + (@$a_ref_H) * $w + 185 , $ysize - ($ybase-$h/4), @head[3]) ;

for (my $i=0; $i<@$a_ref_M; $i++) {
 # get row
 my $r  = $a_ref_M->[$i];

 # get SS description
 my $ss = shift @$r;

 my $w1 = 50 ;

 $p->setfont("ArialBold", 9);
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
 $p->setfont("Arial", 8);
 my $dat = $midata->[$ss] ;

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
$p->text ({align=>"center"}, $xbase + (@$a_ref_H) * $w + 185 , $ysize - ($ybase + ($i*$h+$h/2)+2), $dat->[3]) ;

}

####### draw black boxes ##########
my $NN   = @$a_ref_H ;
for (my $i=0; $i<@$a_ref_M; $i++) {
  my $r  = $a_ref_M->[$i];
  my $w1 = 50 ;
  # ss through entries in current row
  for (my $j=0; $j<@$r; $j++) {
    # create appropriate color
    my $col = "black";
    $p->setcolour($col);
    $p->box({filled => 0},
	    $xbase + $j * $w,      $ysize - ($ybase + $i*$h) ,
	    $xbase + $j * $w + $w, $ysize - ($ybase + ($i*$h+$h)));
  }
}
#
# draw scale bar
#
drawScale($xscale, $ysize / 2 - $yscale , $min, $max, 50, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, "Over-representation", "Under-representation", $scalefont, 2, 20);

my $outeps = "$expfile".".mimatrix.eps" ;
my $outpdf = "$expfile".".mimatrix.pdf" ;
my $ps2pdf = 1;


# output EPS file
print "Outputing EPS file $outeps\n";
$p->output("$outeps");

# convert to PDF
print "Convert to PDF $outpdf\n";
system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");

print "Finished.\n";
exit(0);


sub drawScale {
 my ($x, $y, $min, $max, $res, $p, $xsize, $ysize, $colmap, $A_REF_COLMAP, $upper_text, $lower_text, $scalefont, $h, $w) = @_;

 my $sep = 0;

 $p->setcolour("black");
 $p->setfont("Arial", 8);
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
 $p->setfont("Arial", 8);
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
