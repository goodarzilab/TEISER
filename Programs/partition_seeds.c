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
#include "sequences.h"
#include "matchmaker.h"
#include "dataio.h"
#include "hashtable.h"
#include "readFASTA.h"
#include "read_write_motif.h"
#include "information.h"
#include "mi_library.h"
#include "teiser_functions.h"

int main(int argc, char ** argv) {
  int      i, j, k ;

  char     *seedfile ;
  char     *seedoutdir ;

  s_motif  **seeds ;

  int      seed_cnt = 0 ;
  int      nbseeds = 0 ;

  seedfile       = get_parameter(argc, argv, "-seedfile") ;
  seedoutdir     = get_parameter(argc, argv, "-seedoutdir") ;

  FILE *fptr = fopen ( seedfile, "rb") ;
  seed_cnt = read_motifs( fptr, &seeds ) ;
  printf("%d seeds were loaded...\n", seed_cnt) ;

  FILE* fseed ;
  for (i=0 ; i<seed_cnt ; i++){
    if (i%25000 == 0){
      char seedoutfile[1000] ;
      sprintf(seedoutfile, "%s/filtered.%02d.bin", seedoutdir, i/25000) ;
    
      fseed = fopen(seedoutfile, "wb") ;
      if (!fseed)
        die("Cannot open datafile\n");
      
      int cnt = seed_cnt - i ;
      if (cnt > 25000)
        cnt = 25000 ;
      
      printf ("%d seeds were written to file...\n", i) ;
      fwrite( &cnt, sizeof(int), 1, fseed );
    }
    lcl_write_motif(fseed, seeds[i]) ;
  }
  

  fclose(fseed) ;
  fclose(fptr) ;
  return (0) ;
}
