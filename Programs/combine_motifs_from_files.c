#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>
#include "sys/types.h"
#include "string.h"

#include "nucleotides.h"
#include "structures.h"
#include "dataio.h"
#include "read_write_motif.h"
#include "teiser_functions.h"

int main(int argc, char ** argv) {
  int      i,j ;

  s_motif  **motifs ;
  s_motif  **opt_motifs ;

  int      motif_count = 0 ;

  int numoffiles = atoi(argv[1]) ;
  motifs = (s_motif**) malloc (10000*sizeof(s_motif*)) ;
  for (i=0 ; i<numoffiles ; i++){
    FILE *fptr ;
    fptr = fopen ( argv[2+i], "rb") ;
    if (!fptr){
      printf("Could not open the motif file %d: %s...\n", i, argv[2+i]) ;
      continue ;
    }else{
      printf("Loaded motif file %d: %s...\n", i, argv[2+i]) ;
    }
    s_motif **temp_motifs ;
    int temp_count = read_motifs( fptr, &temp_motifs ) ;
    printf("%d motifs were loaded...\n", temp_count) ;
    for (j=0 ; j<temp_count ; j++){
      motifs[motif_count] = copy_motif(temp_motifs[j]) ;
      motif_count++ ;
    }
    fclose(fptr) ;
  }
  
  FILE *fmotif ;
  /*int opt_count=0 ;
  opt_motifs = (s_motif**) malloc (motif_count*sizeof(s_motif*)) ;
  for (i=0 ; i<motif_count ; i++){
    printf("%d=%d", i, motif_count) ; getchar() ;
    opt_motifs[opt_count++] = copy_motif(motifs[i]) ;
  }
  */
  fmotif = fopen(argv[numoffiles+2], "wb") ;
  if (!fmotif){
    printf("Cannot open motiffile: %s\n", argv[numoffiles+2]);
    return(0) ;
  }

  write_motifs (fmotif, motifs, motif_count) ;

  fclose(fmotif) ;
  return (0) ;
}

