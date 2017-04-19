use lib "$ENV{TEISERDIR}/Scripts";

use Sets;
use Table;
use strict;
use Fasta;
use strict;
use Getopt::Long ;
use Data::Dumper; 

my $fastafile       = undef ;
my $outfasta       = undef ;
my $outfile         = undef ;

GetOptions ('fastafile=s'            => \$fastafile,
	    'outfasta=s'            => \$outfasta,
	    'outfile=s'              => \$outfile,) ;

my $fa = Fasta->new;
$fa->setFile($fastafile);

OF = open($outfasta, "wt") ;
OUT = open($outfile, "wt") ;
my $i=1 ;
while (my $a_ref = $fa->nextSeq()) {
  my ($n, $s) = @$a_ref;
  my $name = sprintf("seq_%d", $i) ;
  print OF ">$name\n$s\n\n";
  print OUT "$n\t$name\n";
}
close(OF) ;
close(OUT) ;
