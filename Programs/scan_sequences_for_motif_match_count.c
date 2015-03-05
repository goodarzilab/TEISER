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

  float    dG_t = -1 ;

  char     *seedfile ;
  char     *fastafile ;
  char     *dataoutfile ;
  char     *reportfile ;

  s_sequence **sequences ;
  s_motif  **seeds ;
  
  int      seq_count ;
  int      seed_cnt = 0 ;

  seqI      si ;
  char*     seq ;
  int       size ;
  char*     name ;
  int       nseqs = 0 ;
  int       seqlen = 0 ;
  
  seedfile         = get_parameter(argc, argv, "-seedfile") ;
  fastafile        = get_parameter(argc, argv, "-fastafile") ;
  if (exist_parameter(argc, argv, "-dG_t")) {
    dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  }

  printf("%s\n%s\n", seedfile, fastafile) ;

  FILE *fptr = fopen ( seedfile, "rb") ;
  if (!fptr){
    printf("Could not open the seed file...\n") ;
    exit(0) ;
  }
  
  seed_cnt = read_motifs( fptr, &seeds ) ;
  printf("%d seeds were loaded...\n", seed_cnt) ;
  fflush(stdout) ;
  fclose(fptr) ;
  
  int hits=0 ;
  int total=0 ;
  if (seqI_open(&si, fastafile) == 0) {
    printf("Error opening fastafile: %s\n", fastafile) ;
    exit(0) ;
  }
  while ((seq = seqI_nextSequence(&si, &name, &size))) {
    total++ ;
    if (total % 100000 == 0){
      printf(".") ;
      fflush(stdout) ;
    }
    seqlen = strlen(seq) ;
    s_sequence* sequence = (s_sequence*) malloc (sizeof(s_sequence)) ;
    sequence->bases = (NUCBIT*) malloc (sizeof(NUCBIT)*seqlen) ;
    sequence->name = strdup(name) ;
    sequence->length = seqlen ;
    for( i=0 ; i<seqlen ; i++ ) {
      sequence->bases[i] = lcl_get_base_id (seq[i]);
      if( !(sequence->bases[i]) ) {
        sequence->bases[i] = _N;
      }
    }
    int match = find_motif_instance ( seeds[0], sequence, dG_t ) ;
    if (match != -1){
      hits++ ;
      fflush(stdout) ;
    }
    free(sequence) ;
  }

  printf("Number of hits = %d / %d", hits, total) ;
  
  fclose(fptr) ;

  return (0) ;
}
