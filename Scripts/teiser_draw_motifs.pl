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

use Sets ;
use strict ;
use Table ;
use Getopt::Long;
use Data::Dumper;
use PostScript::Simple;

my $summaryfile = undef ;

GetOptions ('summaryfile=s'       => \$summaryfile,) ;

my $dir = Sets::dirname($summaryfile)."/Motifs" ;

my $ta = Table->new ;
$ta->loadFile($summaryfile) ;
my $a_M = $ta->getArray() ;
my $a_H = shift @$a_M ; shift @$a_H ;

my $ylen=0;
my $maxlen=0 ;
for (my $i=0 ; $i<@$a_M ; $i++){
  $maxlen = length($a_M->[$i][1]) if (length($a_M->[$i][1])>$maxlen) ;
}

$maxlen+=5 ;
my $xbase = 40 ;
my $ybase = 20 ;
my $xsize = $xbase*3+$maxlen/3*30 ;
my $ysize = $ybase*2+$maxlen/2*30*scalar(@$a_M) ;

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setlinewidth(0.1);
$p->setfont("Times", 5);
$p->setcolour("black");
for (my $i=0 ; $i<@$a_M ; $i++){
  print $a_M->[$i][2], "\t", $a_M->[$i][3], "\n" ;

  my $myii = $a_M->[$i][0] ;
  my $e  = new PostScript::Simple::EPS(file => "$dir/$myii\_motif.eps");
  
  my $h = $maxlen*4 ;
  # get height
  my $eh = $e->height;
  my $ew = $e->width;
  # height must be $h, so scale down to $h = k * LO
  $e->scale($h / $eh);
  my $ew_new = int(0.5 + $ew * $h / $eh);
  
  $p->_add_eps($e, $xbase/5, $ysize-(($maxlen)*30/2*($i)+2*$ybase)); 

  my $c = $a_M->[$i][2] ;
  my $s = $a_M->[$i][3] ;

  my @cc = split(//, "5$c"."3") ;
  my @ss = split(//, "($s)") ;
  &drawRNA ($p, $xbase*2, $ysize-(($maxlen)*30/2*($i+1)+$ybase/2), \@ss, \@cc, 0, $#ss, 0, 0, 1, 0) ;
}


my $outeps = "$summaryfile" ;
$outeps =~ s/\.summary/.motifs.eps/ ;
my $outpdf = $outeps ;
$outpdf =~ s/\.eps/.pdf/ ;

$p->output("$outeps");
system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");
print "Finished.\n";



sub drawRNA{	
  my ($p, $xbase, $ybase, $structure, $bases, $l, $r, $lx, $ly, $rx, $ry) = @_;

  my $L = 0 ;
  my $count = 2 ;
  for (my $i=$l+1; $i<$r ; $i++){
    $L-- if ($structure->[$i] eq ")") ;
    $count++ if ($L==0) ;
    $L++ if ($structure->[$i] eq "(") ;
  }

  my $th = 2 * 3.14159 / $count;
  my $R = 1 / (2*sin($th/2));
  my $h = $R * cos($th/2);
  my ($cx, $cy) = ((($lx+$rx)/2.0)+$h*($ly-$ry),
		   (($ly+$ry)/2.0)+$h*($rx-$lx)); 
  my $deg = atan2($ly-$cy,$lx-$cx);
  
  my ($i2,$x2,$y2) = ($l,$lx,$ly);
  
  for (my $i=$l+1; $i<=$r; $i++) {
    $L-- if ($structure->[$i] eq ")");
    if ($L==0) {
      $deg -= $th;
      my ($x,$y) = ($cx+$R*cos($deg), $cy+$R*sin($deg));
      $p->setcolour(118, 118, 118) ;
      $p->setlinewidth (0.5) ;
      $p->line($x2*10+$xbase,$y2*10+$ybase,$x*10+$xbase,$y*10+$ybase);
      $p->setlinewidth (0.1) ;
      $p->setcolour("black") ;
      drawRNA($p, $xbase, $ybase, $structure, $bases, $i2,$i,$x2,$y2, $x, $y) if ($structure->[$i] eq ")");
      $p->setcolour("white") ;
      $p->box({filled=>1},$x2*10-2.5+$xbase,$y2*10+$ybase-2,$x2*10+2.5+$xbase,$y2*10+$ybase+3) ;
      $p->setcolour("black") ;
      $p->box({filled=>0},$x2*10-2.5+$xbase,$y2*10+$ybase-2,$x2*10+2.5+$xbase,$y2*10+$ybase+3) ;
      $p->text({align=>"center"}, $x2*10+$xbase,$y2*10+$ybase-1, $bases->[$i2]);
      ($x2,$y2)=($x,$y);
      $i2 = $i;
    }
    $L++ if ($structure->[$i] eq "(");
  }
  $p->setcolour(118, 118, 118) ;
  $p->setlinewidth (0.5) ;
  $p->line($x2*10+$xbase,$y2*10+$ybase,$rx*10+$xbase,$ry*10+$ybase);
  $p->setlinewidth (0.1) ;
  $p->setcolour("black") ;
  $p->setcolour("white") ;
  $p->box({filled=>1},$x2*10-2.5+$xbase,$y2*10+$ybase-2,$x2*10+2.5+$xbase,$y2*10+$ybase+3) ;
  $p->setcolour("black") ;
  $p->box({filled=>0},$x2*10-2.5+$xbase,$y2*10+$ybase-2,$x2*10+2.5+$xbase,$y2*10+$ybase+3) ;
  $p->text({align=>"center"}, $x2*10+$xbase,$y2*10+$ybase-1, $bases->[$r-1]);
  my %St = ("GU"=>1,"UG"=>1,"AU"=>2,"UA"=>2,"CG"=>3,"GC"=>3) ;
  my $w = $St{$bases->[$l].$bases->[$r]}/2 ;
  $p->setcolour(118, 118, 118) ;
  $p->line($lx*10+$xbase,$ly*10+$ybase,$rx*10+$xbase,$ry*10+$ybase) unless ($lx==0 || $ly==0);   ## 3-5 pair          
  $p->setcolour("black") ;
  $p->setcolour("white") ;
  $p->box({filled=>1},$rx*10-2.5+$xbase,$ry*10+$ybase-2,$rx*10+2.5+$xbase,$ry*10+$ybase+3) ;
  $p->setcolour("black") ;
  $p->box({filled=>0},$rx*10-2.5+$xbase,$ry*10+$ybase-2,$rx*10+2.5+$xbase,$ry*10+$ybase+3) ;

  $p->text({align=>"center"}, $rx*10+$xbase,$ry*10+$ybase-1, $bases->[$r]);
}
