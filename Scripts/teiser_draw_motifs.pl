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
print Dumper $dir ;

my $ta = Table->new ;
$ta->loadFile($summaryfile) ;
my $a_M = $ta->getArray() ;
my $a_H = shift @$a_M ; shift @$a_H ;

my $ylen=0;
my $maxlen=0 ;
for (my $i=0 ; $i<@$a_M ; $i++){
  $maxlen = length($a_M->[$i][1]) if (length($a_M->[$i][1])>$maxlen) ;
}

my $col = 5 ;
my $cnt = scalar(@$a_M) ;
$maxlen+=5 ;
my $xbase = 140 ;
my $ybase = 120 ;
my $xsize = $xbase*3+$maxlen/3*30*$col ;
my $ysize = $ybase*2+$maxlen/2*30*$cnt/$col ;

my $p = new PostScript::Simple(xsize     => $xsize,
			       ysize     => $ysize,
			       colour    => 1,
			       eps       => 1,
			       units     => "pt");
$p->setlinewidth(0.1);
$p->setfont("Times", 5);
$p->setcolour("black");

my $counti=0 ;
my $countj=-1 ;
for (my $i=0 ; $i<@$a_M ; $i++){
  $countj+=1 ;
  print $a_M->[$i][2], "\t", $a_M->[$i][3], "\t$counti\t$countj\n" ;

  my $myii = $a_M->[$i][0] ;
  my $mo = $a_M->[$i][2] ;
  my $mo = &seq_to_re($mo);
  my $mo = Sets::myre2wm($mo);

  open OUT, ">$dir/$myii\_motif.txt" or die "cannot open $dir/$myii\_motif.txt\n";
  print OUT $mo;
  close OUT;
    
  system("$scriptdir/weblogo/seqlogo -f $dir/$myii\_motif.txt -F EPS  -a -c -M -n -Y -w 5 -h 3 > $dir/$myii\_motif.eps");

  my $e  = new PostScript::Simple::EPS(file => "$dir/$myii\_motif.eps");
  
  my $h = $maxlen ;
  # get height
  my $eh = $e->height;
  my $ew = $e->width;
  # height must be $h, so scale down to $h = k * LO
  $e->scale($h / $eh);
  my $ew_new = int(0.5 + $ew * $h / $eh);
  
  $p->_add_eps($e, $xbase+$maxlen/3*30*($countj), $ysize-(($maxlen)*30/2*($counti-1)+2*$ybase)); 

  my @r = (255,113,113) ;
  my @g = (113,255,113) ;
  my @y = (255,255,113) ;
  my @b = (0,195,255) ;
  my @w = (255,255,255) ;
  my $colors = {N        => \@w,
          U        => \@r,
          A        => \@g,
          G        => \@y,
          C        => \@b} ;

  $colors->{Y} = &addColors($colors, "U", "C") ;
  $colors->{R} = &addColors($colors, "A", "G") ;
  $colors->{K} = &addColors($colors, "U", "G") ;
  $colors->{M} = &addColors($colors, "A", "C") ;
  $colors->{S} = &addColors($colors, "G", "C") ;
  $colors->{W} = &addColors($colors, "U", "A") ;
  $colors->{B} = &addColors($colors, "K", "C") ;
  $colors->{D} = &addColors($colors, "U", "R") ;
  $colors->{H} = &addColors($colors, "U", "M") ;
  $colors->{V} = &addColors($colors, "G", "M") ;

  my $c = $a_M->[$i][2] ;
  my $s = $a_M->[$i][3] ;

  my @cc = split(//, "5$c"."3") ;
  my @ss = split(//, "($s)") ;
  
  &drawRNA ($p, $xbase+$maxlen/3*30*($countj)+20, $ysize-(($maxlen)*30/2*($counti+1)+$ybase/2), \@ss, \@cc, 0, $#ss, 0, 0, 1, 0, $colors) ;
  if ($countj>$col) {
    $countj=-1 ;
    $counti += 1 ;
  }
}


my $outeps = "$summaryfile" ;
$outeps =~ s/\.summary// ;
$outeps = "$outeps.motifs.eps" ;
my $outpdf = $outeps ;
$outpdf =~ s/\.eps/.pdf/ ;

$p->output("$outeps");
system("ps2pdf -dEPSCrop -dAutoRotatePages=/None $outeps $outpdf");
print "Finished.\n";



sub drawRNA{  
  my ($p, $xbase, $ybase, $structure, $bases, $l, $r, $lx, $ly, $rx, $ry, $colors) = @_;

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
      drawRNA($p, $xbase, $ybase, $structure, $bases, $i2,$i,$x2,$y2, $x, $y, $colors) if ($structure->[$i] eq ")");
      my $col = $colors->{$bases->[$i2]} ;
      $p->setcolour(@$col) ;
      #$p->box({filled=>1},$x2*10-2.5+$xbase,$y2*10+$ybase-2,$x2*10+2.5+$xbase,$y2*10+$ybase+3) ;
      $p->circle({filled=>1},$x2*10+$xbase,$y2*10+$ybase,3) ;
      $p->setcolour("black") ;
      #$p->box({filled=>0},$x2*10-2.5+$xbase,$y2*10+$ybase-2,$x2*10+2.5+$xbase,$y2*10+$ybase+3) ;
      $p->circle({filled=>0},$x2*10+$xbase,$y2*10+$ybase,3) ;
      $p->text({align=>"center", rotate=>"0"}, $x2*10+$xbase,$y2*10+$ybase-1.5, $bases->[$i2]);
      #$p->text({align=>"center", rotate=>"90"}, $x2*10+$xbase+1.5,$y2*10+$ybase+0.5, $bases->[$i2]);
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
  my $col = $colors->{$bases->[$r-1]} ;
  $p->setcolour(@$col) ;
  #$p->box({filled=>1},$x2*10-2.5+$xbase,$y2*10+$ybase-2,$x2*10+2.5+$xbase,$y2*10+$ybase+3) ;
  $p->circle({filled=>1},$x2*10+$xbase,$y2*10+$ybase,3) ;
  $p->setcolour("black") ;
  #$p->box({filled=>0},$x2*10-2.5+$xbase,$y2*10+$ybase-2,$x2*10+2.5+$xbase,$y2*10+$ybase+3) ;
  $p->circle({filled=>0},$x2*10+$xbase,$y2*10+$ybase,3) ;
  $p->text({align=>"center", rotate=>"0"}, $x2*10+$xbase,$y2*10+$ybase-1.5, $bases->[$r-1]);
  #$p->text({align=>"center", rotate=>"90"}, $x2*10+$xbase+1.5,$y2*10+$ybase+0.5, $bases->[$r-1]);
  my %St = ("GU"=>1,"UG"=>1,"AU"=>2,"UA"=>2,"CG"=>3,"GC"=>3) ;
  my $w = $St{$bases->[$l].$bases->[$r]}/2 ;
  $p->setcolour(118, 118, 118) ;
  $p->line($lx*10+$xbase,$ly*10+$ybase,$rx*10+$xbase,$ry*10+$ybase) unless ($lx==0 || $ly==0);   ## 3-5 pair          
  $p->setcolour("black") ;
  my $col = $colors->{$bases->[$r]} ;
  $p->setcolour(@$col) ;
  #$p->box({filled=>1},$rx*10-2.5+$xbase,$ry*10+$ybase-2,$rx*10+2.5+$xbase,$ry*10+$ybase+3) ;
  $p->circle({filled=>1},$rx*10+$xbase,$ry*10+$ybase,3) ;
  $p->setcolour("black") ;
  #$p->box({filled=>0},$rx*10-2.5+$xbase,$ry*10+$ybase-2,$rx*10+2.5+$xbase,$ry*10+$ybase+3) ;
  $p->circle({filled=>0},$rx*10+$xbase,$ry*10+$ybase,3) ;    

  $p->text({align=>"center", rotate=>"0"}, $rx*10+$xbase,$ry*10+$ybase-1.5, $bases->[$r]);
  #$p->text({align=>"center", rotate=>"90"}, $rx*10+$xbase+1.5,$ry*10+$ybase+0.5, $bases->[$r]);
}

sub addColors {
  my ($c, $n1, $n2) = @_ ;
  my @a1 = @{$c->{$n1}} ;
  my @a2 = @{$c->{$n2}} ;
  my @a3 = (($a1[0]+$a2[0])/2, ($a1[1]+$a2[1])/2, ($a1[2]+$a2[2])/2) ;
  $a3[0]==255 if ($a3[0]>255) ;
  $a3[1]==255 if ($a3[1]>255) ;
  $a3[2]==255 if ($a3[2]>255) ;
  return \@a3 ;
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
