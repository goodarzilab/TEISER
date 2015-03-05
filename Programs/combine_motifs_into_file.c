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
#include "statistics.h"
#include "hashtable.h"
#include "readFASTA.h"
#include "read_write_motif.h"
#include "information.h"
#include "mi_library.h"
#include "teiser_functions.h"

int main(int argc, char ** argv) {
  int      i,j ;
  
  int      nfiles =0;
  char*    buff = (char*)malloc(100000 * sizeof(char));;
  int      motif_count =0 ;
  char     *motiffile ;
  char     **files = (char**)malloc(1000 * sizeof(char*)); ;
  char     *outfile ;

  s_motif  **motifs ;
  motifs = (s_motif**) malloc (10000000*sizeof(s_motif*)) ;

  motiffile        = get_parameter(argc, argv, "-motiffile") ;
  outfile          = get_parameter(argc, argv, "-outfile") ;

  FILE*  fp = fopen(motiffile, "rt");
  if (!fp) {
    printf("could not open set data %s\n", motiffile);
  }

  while (!feof(fp)) {
    fscanf(fp, "%s\n", buff);
    printf("%s\n", buff) ;
    FILE *fptr = fopen ( buff, "rb") ;
    
    if (!fptr){
      printf("Could not open the seed file: %s\n", files[i]) ;
      exit(0) ;
    }

    s_motif  **tmp_motifs ;
    int count = read_motifs( fptr, &tmp_motifs ) ;
    printf("%d seeds loaded\n", count) ;
    
    for (j=0 ; j<count ; j++){
      printf(".") ; fflush(stdout) ;
      motifs[motif_count++] = copy_motif(tmp_motifs[j]) ;
    }
    for (j=0 ; j<count ; j++){
      free(tmp_motifs[j]->phrases) ;
    }
    free(tmp_motifs) ;
    fclose(fptr) ;
  }

  printf("%d seeds were loaded...\n", motif_count) ;
  fflush(stdout) ;

  FILE* fmotif = fopen(outfile, "wb") ;
  if (!fmotif)
    die("Cannot open motif outfile\n");

  write_motifs (fmotif, motifs, motif_count) ;
  fclose(fmotif) ;

  return (0) ;
}
