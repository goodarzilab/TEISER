package AggloClust;
use Sets;
use strict;

sub new {
    my ($class) = @_;
    my ($self) = {};
    $self->{DIST}       = [];
    $self->{MAXNBCLUST} = undef;
    $self->{MIND}       = undef;
    bless $self;
    return $self;
}


sub setMaxNbClusters {
  my ($self, $m) = @_;
  $self->{MAXNBCLUST} = $m;
}


sub setDistanceMatrix {
  my ($self, $m) = @_;
  $self->{DIST} = $m;
}



sub setMinD {
  my ($self, $m) = @_;
  $self->{MIND} = $m;
}


sub agglomerate_using_avg_linkage {

  my ($self) = @_;

  my $a_ref_d = $self->{DIST};
 
  if (!defined($self->{MAXNBCLUST}) && !defined($self->{MIND})) {
    die "Please set a stopping condition.\n";
  }

  my $n       = @$a_ref_d;
  
  #print "Got $n objects\n";
  
  my @CLASSES = ();
  my @LIST    = ();
  for (my $i=0; $i<$n; $i++) {
    my @a = ($i);
    push @CLASSES, \@a;
    push @LIST   , $i;
  }
 
  my @DIST = ();
  for (my $i=0; $i<$n; $i++) {
    for (my $j=0; $j<$n; $j++) {
      $DIST[$i][$j] = $a_ref_d->[$i]->[$j];
    }
  }

  my @NEWLIST = ();

  while (1) {
    my $a_ref_ijd = $self->find_min( \@DIST, \@LIST );
    my ($i, $j, $min) = @$a_ref_ijd;

    print "best merge is $i, $j, $min\n";

    if ((defined($self->{MIND}) && ($min > $self->{MIND})) || (defined($self->{MAXNBCLUST}) && (scalar(@LIST) == $self->{MAXNBCLUST}))) {
      my @NEWCLASSES = ();
      foreach my $k (@LIST) {
	push @NEWCLASSES, $CLASSES[$k];
      }

      #
      # order the classes here ? first class is tightest ? (then use prior ??) P(X)
      #

      print "order classes\n";

      my $best_myi = undef; my $best_d = 1000000;
      my $myi      = 0;
      foreach my $r (@NEWCLASSES) {
	if (scalar(@$r) > 1) {
	  my $d = $self->get_pairwise_avg_dist($r, $r) / (scalar(@$r) / $n);
	  print scalar(@$r) . " - " . $d . "\n";
	  if ($d < $best_d) {
	    $best_myi = $myi;
	    $best_d   = $d;
	  }
	}
	$myi ++;
      }
      
      my $tmpa = $NEWCLASSES[0];
      $NEWCLASSES[0] = $NEWCLASSES[$best_myi];
      $NEWCLASSES[$best_myi] = $tmpa;

      return \@NEWCLASSES;
    }

    # merge i and j into i
    @{$CLASSES[$i]} = ( @{$CLASSES[$i]}, @{$CLASSES[$j]} );

    #print "Merge into " . join(" " , @{$CLASSES[$i]}) . "\n";

    # remove $j from list of classes
    @NEWLIST = ();
    foreach my $k (@LIST) {
      if ($k != $j) {
	push @NEWLIST, $k;
      }
    }
    @LIST = @NEWLIST;

    # update dist from new i to all other classes
    foreach my $k (@LIST) {
      next if ($k == $i); 
      $DIST[$i][$k] = $self->get_pairwise_avg_dist($CLASSES[$i], $CLASSES[$k]);
      #print "Updated dist[$i][$k] = $DIST[$i][$k]\n";
      $DIST[$k][$i] = $DIST[$i][$k];
    } 

    #<STDIN>;
  }

}

sub get_pairwise_avg_dist {
  my ($self, $a1, $a2) = @_;
  
  my $a_ref_d = $self->{DIST};

  my $n1 = @$a1;
  my $n2 = @$a2;
  
  my $sum = 0;

  for (my $i=0; $i<$n1; $i++) {
    for (my $j=0; $j<$n2; $j++) {
      $sum += $a_ref_d->[ $a1->[$i] ]->[ $a2->[$j] ];
    }
  }
  
  $sum = $sum / ($n1 * $n2);
  
  return $sum;
}

sub find_min {
  my ($self, $a_ref_d, $a_ref_l) = @_;


  my $min     = 100000;
  my $min_i   = undef;
  my $min_j   = undef;

  foreach my $i (@$a_ref_l) { 
    foreach my $j (@$a_ref_l) { 
      next if ($i == $j);
      if ($a_ref_d->[$i]->[$j] < $min) {
	$min   = $a_ref_d->[$i]->[$j];
	$min_i = $i;
	$min_j = $j;
      }
    }
  }
  
  my $a = [$min_i, $min_j, $min];
  return $a;
}


1;
