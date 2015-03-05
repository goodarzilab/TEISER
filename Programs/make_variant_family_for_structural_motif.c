#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>

#include "dataio.h"
#include "structures.h"
#include "nucleotides.h"
#include "create_motifs.h"
#include "read_write_motif.h"

int main(int argc, char ** argv) {
  int      i,j,k,l,m,n,o,p ;
  char     *motifoutfile ;
  s_motif  **motif_vars ;

  motifoutfile     = get_parameter(argc, argv, "-motifoutfile") ;

  int nvars=0 ;
  motif_vars = (s_motif**) malloc (10000*sizeof(s_motif*)) ;
  int nmotifs=0 ;
  for (i=4 ; i<10 ; i++){//loop lengthh
    for (j=4 ; j<8 ; j++){//stem length
      motif_vars[nvars] = (s_motif*) malloc (sizeof(s_motif)) ;
      lcl_create_structure( motif_vars[nvars], j, i );
      for (k=0 ; k<j ; k++){
        motif_vars[nvars]->phrases[k].base = _S ;
      }
      nvars++ ;
    }
  }

  printf("%i seeds written...\n", nvars) ;
  FILE *fmotif = fopen(motifoutfile, "wb") ;
  write_motifs (fmotif, motif_vars, nvars) ;
  fclose(fmotif) ;
  return (0) ;
}
