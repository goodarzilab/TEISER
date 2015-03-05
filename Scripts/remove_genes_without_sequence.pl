use lib "$ENV{TEISERDIR}/Scripts";


use Sets;
use Table;
use strict;
use Fasta;
use strict;
use Getopt::Long ;

my $listfile       = undef ;
my $expfile         = undef ;
my $outfile         = undef ;

GetOptions ('expfile=s'              => \$expfile,
	          'listfile=s'             => \$listfile,
	          'outfile=s'              => \$outfile,
	         ) ;

my $ta = Table->new ;
$ta->loadFile($listfile) ;
my $list = $ta->getHash ;

open I, "< $expfile" ;
open O, "> $outfile" ;
my $l = <I> ;
print O "RefSeq\tvalue\n" ;

while(my $l = <I>){
  $l =~ s/\s+$// ;
  my ($ref, $v) = split(/\t/, $l) ;
  if (defined $list->{$ref}){
    print O "$ref\t$v\n" ;
  }
}
close(I) ;
close(O) ;