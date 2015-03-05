package Stat;

use Exporter;
@ISA = (Exporter);
@EXPORT = qq(sd);

use strict;
use warnings;
use Math::Trig;

# ----------------------- #
#   Basic stats values    #
# ----------------------- #

sub sd{
    my $ref_data = shift;
    my @data = (@$ref_data);
   
    my $sum= 0;
    my $sum_sq= 0;
    
    foreach my $datum(@data){
        $sum += $datum;
        $sum_sq += $datum**2;
    }
    
    my $n = $#data+1;
    my $avg = $sum/$n;
    my $sd = sqrt( ($sum_sq - $sum**2 / $n)/($n-1) );
    
    return($n,$sum,$sum_sq,$avg,$sd);
}


# ----------------------- #
#  Normal distribution p  #
# ----------------------- #

sub normDist{  # input $z, output $p
    
    my $z = shift;

    my @x = (
        '0.076526521133497333754640409398',
        '0.227785851141645078080496195368',
        '0.373706088715419560672548177024',
        '0.510867001950827098004364050955',
        '0.636053680726515025452836696226',
        '0.746331906460150792614305070355',
        '0.839116971822218823394529061701',
        '0.912234428251325905867752441203',
        '0.963971927277913791267666131197',
        '0.993128599185094924786122388471');

    my @w = (
        '0.152753387130725850698084331955',
        '0.149172986472603746787828737001',
        '0.142096109318382051329298325067',
        '0.131688638449176626898494499748',
        '0.118194531961518417312377377711',
        '0.101930119817240435036750135480',
        '0.083276741576704748724758143222',
        '0.062672048334109063569506535187',
        '0.040601429800386941331039952274',
        '0.017614007139152118311861962351');


    my $p1 = $z;
    if ( $z < 0 ) { $p1 = -1*$z; }

    my $xi = $p1 ; 
    my $xe = 8 + $p1;
    my $s = 0;

    for (my $i=0;$i<10;$i++){
        my $a = 0.5 * ( ($xe - $xi) * $x[$i] + ($xe + $xi ));
        my $b = 0.5 * (-($xe - $xi) * $x[$i] + ($xe + $xi ));
        $s = $s
            + $w[$i] * (1/sqrt(2*pi)*exp(1)**(-0.5*$a**2))
            + $w[$i] * (1/sqrt(2*pi)* exp(1)**(-0.5*$b**2))
        ;
     }
     
    my $p = 0.5 * $s * ($xe - $xi)*2;
    return $p;
}


# ----------------------- #
#         F tests         #
# ----------------------- #

sub ftest(){

    my ($ref_data1, $ref_data2) = @_;
    my @g1s = (@$ref_data1);
    my @g2s = (@$ref_data2);

    my ($n1,$sum1,$sum_sq1,$avg1,$sd1) = &sd(\@g1s);
    my ($n2,$sum2,$sum_sq2,$avg2,$sd2) = &sd(\@g2s);

    my $df1 = $n1-1;
    my $df2 = $n2-1;
    my ($dfl,$dfs,$f12);

    my $f1 = ($sd1**2)/($sd2**2);
    my $f2 = ($sd2**2)/($sd1**2);;

    if ($df1>$df2||$df1==$df2){$dfl=$df1;$dfs=$df2;}
    else{$dfl=$df2;$dfs=$df1;}

    if ($f1>$f2){$f12=$f1;}else{$f12=$f2;}

    #  get p value 
    my $pf = &pforf($dfl,$dfs,$f12);
    return ($f12,$pf);

}


sub pforf{   # input $dfl,$dfs,$f12, output $pf
    my ($dfl,$dfs,$f12) = @_;
    my $aa1 = 2/(9*$dfl);
    my $aa2 = 2/(9*$dfs);

    my $w = $f12**(1/3);
    my $w1 = $w+$aa1-$w*$aa2-1;
    my $w2 = $aa2*$w*$w+$aa1;
    my $z = $w1/sqrt($w2);

    $z = $z*(1+0.8*($z**4/$dfs**3)) if ($dfs<=3);
    
    my $pf = &normDist($z);
    $pf /= 2;
    return $pf;
}


# ----------------------- #
#         t tests         #
# ----------------------- #

sub unpairedt(){

    my ($ref_data1, $ref_data2) = @_;
    my @g1s = (@$ref_data1);
    my @g2s = (@$ref_data2);
    
    if ($#g1s<3||$#g2s<3){
        return (0,0);
    }

    my ($n1,$sum1,$sum_sq1,$avg1,$sd1) = &sd(\@g1s);
    my ($n2,$sum2,$sum_sq2,$avg2,$sd2) = &sd(\@g2s);

    my $utcom1 = ($n1-1)*$sd1*$sd1+($n2-1)*$sd2*$sd2;
    my $utcom2 = $n1+$n2-2;
    my $utcom3 = $utcom1/$utcom2;
    my $utcom4 = 1/$n1+1/$n2;

    my $semd = sqrt($utcom3*$utcom4);
    my $d = $avg1-$avg2;

    my $t = $d/$semd;
    my $df = $utcom2;
    my $tabs = abs($t);

    my $utp = &pfort($df,$tabs);
    return ($t,$utp);
}


sub pairedt(){

    my ($ref_data1, $ref_data2) = @_;
    my @g1s = (@$ref_data1);
    my @g2s = (@$ref_data2);
    
    if ($#g1s<3||$#g1s!=$#g2s){
        return (0,0);
    }

    my $sub=0;
    my $sum_sub=0;
    my $sum_sq_sub=0;

    for (my $i=0;$i<$#g1s+1;$i++){
        $sub = $g1s[$i]-$g2s[$i];
        $sum_sub += $sub;
        $sum_sq_sub += $sub**2;
    }

    my $n1 = $#g1s+1;

    my $avg_sub = $sum_sub/$n1;

    my $sd_sub = sqrt( ($sum_sq_sub - $sum_sub**2 / $n1)/($n1-1) );
    my $t = $avg_sub/($sd_sub/sqrt($n1));

    my $df = $n1-1;
    my $tabs = abs($t);
    
    my $ptp = &pfort($df,$tabs);
    return ($t,$ptp);
}


sub pfort{  # input $df,$tabs, output $pp
    my ($df,$tabs) = @_;
    my $w1 = (1-2/(3+8*$df));
    my $wx = $tabs*$tabs/$df+1;
    my $wy = log($wx);
    my $w2 = $df*$wy;
    my $z = $w1*sqrt($w2);
    
    my $pp = normDist($z);
    
    return $pp;
}

# ----------------------- #
#      Correlations       #
# ----------------------- #

sub pearson{
    my ($ref_data1, $ref_data2) = @_;
    my @g1s = (@$ref_data1);
    my @g2s = (@$ref_data2);
    
    if ($#g1s<3||$#g2s<3||$#g1s!=$#g2s){
        return (0,0,0);
    }

    my ($n1,$sum1,$sum_sq1,$avg1,$sd1) = &sd(\@g1s);
    my ($n2,$sum2,$sum_sq2,$avg2,$sd2) = &sd(\@g2s);

    my $subx_sum=0;
    my $sub_sq1_sum=0;
    my $sub_sq2_sum=0;

    for (my $i=0;$i<$#g1s+1;$i++){
        # upper value
        my $sub1 = $g1s[$i]-$avg1;
        my $sub2 = $g2s[$i]-$avg2;
        my $subx = $sub1*$sub2;
        $subx_sum += $subx;
        # lower value
        my $sub_sq1 = $sub1**2;
        my $sub_sq2 = $sub2**2;
        $sub_sq1_sum += $sub_sq1;
        $sub_sq2_sum += $sub_sq2;
    }

    my $r = $subx_sum / sqrt($sub_sq1_sum*$sub_sq2_sum);
    $r = 0.9999 if (1 == $r) ;
    my $z0 = (1+$r)/(1-$r);
    my $z = 0.5*log($z0) / sqrt(1/($n1-3));

    my $p = &normDist($z);
    $p=$p*2;

    return ($r,$z,$p);
}


sub spearman{
    my ($ref_data1, $ref_data2) = @_;
    my @g1s = (@$ref_data1);
    my @g2s = (@$ref_data2);
    
    if ($#g1s<3||$#g2s<3||$#g1s!=$#g2s){
        return (0,0,0);
    }
    
    my (@g1o,@g2o);
    for (0..$#g1s){
        push @g1o, [$_+1, $g1s[$_], 0];
        push @g2o, [$_+1, $g2s[$_], 0];
    }
    
    # --- Rank by value --- #
    @g1o= sort {$a->[1] <=> $b->[1]} @g1o;
    @g2o= sort {$a->[1] <=> $b->[1]} @g2o;
    
    # --- Put rank # --- #
    $g1o[$_-1]->[2] = $_ for (1..$#g1o+1);
    $g2o[$_-1]->[2] = $_ for (1..$#g2o+1);
    
    my ($ref_g1o, $mk1) = &ranker(\@g1o);
    my ($ref_g2o, $mk2) = &ranker(\@g2o);
    
    @g1o = @$ref_g1o;
    @g2o = @$ref_g2o;
    
    # --- Back to original order --- #
    @g1o= sort {$a->[0] <=> $b->[0]} @g1o;
    @g2o= sort {$a->[0] <=> $b->[0]} @g2o;
    
    my (@rank1, @rank2);
    push @rank1, $_->[2] foreach @g1o;
    push @rank2, $_->[2] foreach @g2o;
    
    my ($r,$z,$p) = &pearson(\@rank1,\@rank2);
    $p /= 2;
    return ($r,$z,$p);

}


# ----------------------- #
#      Nonparametric      #
# ----------------------- #

sub ranker{
    my $ref_array= shift;
    my @all = @$ref_array;

    my $flg = 0;
    my $same_n = 1;
    my $order_sum = 0;
    my @ks;
    for my $data_n (0..$#all-1){
        
        if ($all[$data_n]->[1] == $all[$data_n+1]->[1]){
            $order_sum += $all[$data_n]->[2];
            $same_n++;
            $flg =1;
        }else{$flg = 0;}
        
        if ($flg==0 && $order_sum>0){
            $order_sum += $all[$data_n]->[2];
            my $order_avg = $order_sum/$same_n;
            push @ks, $same_n;
            for (1..$same_n){
                $all[$data_n-$_+1]->[2] = $order_avg;;
            }
            $order_sum = 0;
            $same_n = 1;
        }
    }
    return (\@all,\@ks);
}


sub mannwhitney{
    my ($ref_data1, $ref_data2) = @_;
    my @g1s = (@$ref_data1);
    my @g2s = (@$ref_data2);
    
    if ($#g1s<3||$#g2s<3){
        return (0,0,0);
    }
    
    my $datnum = ($#g1s+$#g2s)+2;
    my $n1 = $#g1s+1;
    my $n2 = $#g2s+1;
    
    my @all;
    push @all, [1,$_,0] foreach @g1s;
    push @all, [2,$_,0] foreach @g2s;
    @all = sort {$a->[1] <=> $b->[1]} @all;
    
    # --- Put temporal order --- #
    $all[$_-1]->[2] = $_ for (1..$datnum);
    
    my ($ref_all,$ref_ks) = &ranker(\@all);
    @all = @$ref_all;
    my @ks = @$ref_ks;

    # --- Calc Ta Tb values --- #
    my ($t1, $t2);
    for (0..$datnum-1){
        $t1 += $all[$_]->[2] if $all[$_]->[0] == 1;
        $t2 += $all[$_]->[2] if $all[$_]->[0] == 2;
    }
    
    # --- Calc U value --- #
    
    my $u1 = $t1-$n1*($n1+1)/2;
    my $u2 = $t2-$n2*($n2+1)/2;
        
    my $ules = $u1<$u2 ? $u1 : $u2;
    
    # --- Calc z value --- #
    my $z;

    if ($#ks>0){
        my $k_val;
        $k_val += ($_**3-$_)/12 foreach @ks;
        $z = ( $ules-($n1*$n2)/2 ) /
        sqrt( ( $n1*$n2/($datnum*($datnum-1) ) *
               ( ($datnum**3-$datnum)/12-$k_val ) )
             );
    }else{
        $z = ( $ules-($n1*$n2)/2 ) / sqrt( $n1*$n2*($n1+$n2+1)/12 );
    }
    
    my $pmw = &normDist($z);
    return ($ules,$z,$pmw);
}


sub wilcoxson{
    my ($ref_data1, $ref_data2) = @_;
    my @g1s = (@$ref_data1);
    my @g2s = (@$ref_data2);
    
    return if ($#g1s!=$#g2s||$#g1s<3||$#g2s<3);
    
    my $m = 0;
    my $k =0;
    my @subs;
    for (0..$#g1s){
        my $sub = $g1s[$_]-$g2s[$_];
        push @subs, $sub unless $sub==0;
        $m++ unless $sub==0;
        $k++ if $sub==0;
    }
    
    @subs = sort {abs($a) <=> abs($b)} @subs;
    
    my $t_plus = 0;
    my $t_minus = 0;
    for (0..$#subs){
        $t_plus += abs($_+1) if $subs[$_] > 0;
        $t_minus += abs($_+1) if $subs[$_] < 0;
    }
    
    my $tles = 0;
    $tles = $t_plus<$t_minus ? $t_plus : $t_minus;
    
    my $z;
    my $kval = $k**3-$k;
    
    if ($m==$#g1s+1){
        $z = ( $tles - $m*($m+1)/4 ) /
            sqrt ( ($m*($m+1)*(2*$m+1)) / 24 );
    }else{
        $z = ( $tles - $m*($m+1)/4 ) /
        sqrt ( ($m*($m+1)*(2*$m+1)-$kval/2) / 24 );
    }
    
    my $pwl = &normDist($z);
    return ($tles,$z,$pwl);
}


return 1;
