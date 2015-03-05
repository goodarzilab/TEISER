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
  int      i ;
  float    dG_t = -1 ;
  s_motif  **motifs ;
  int      rnd_fasta = 0 ;
  
  char     **seq_names ;
  int      seq_count ;
  int      t_seq_count ;
  int      motif_count = 0 ;

  char     *seedfile ;
  char     *inseq ;
  seedfile         = get_parameter(argc, argv, "-seedfile") ;
  inseq    = get_parameter(argc, argv, "-inseq") ;
  if (exist_parameter(argc, argv, "-dG_t")) {
    dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  }
  
  FILE *fptr = fopen ( seedfile, "rb") ;
  if (!fptr){
    printf("Could not open the seed file: %s\n", seedfile) ;
    exit(0) ;
  }
  
  s_sequence *seq = (s_sequence*) malloc (10*sizeof(s_sequence)) ;
  int seqlen = strlen(inseq) ;
  
  seq->name=strdup("test") ;
  seq->bases=(NUCBIT*) malloc (sizeof(NUCBIT)*seqlen) ;
  seq->length = seqlen ;
  for( i=0 ; i<seqlen ; i++ ) {
    seq->bases[i] = lcl_get_base_id (inseq[i]);
    if( !(seq->bases[i]) ) {
      seq->bases[i] = _N ;
    }
  }
  fflush(stdout) ;
  motif_count = read_motifs( fptr, &motifs ) ;
  s_motif *motif = copy_motif(motifs[0]) ;
  fflush(stdout) ;
  fclose(fptr) ;

  int match=0 ;
  for (i=0 ; i<motif_count ; i++){
    int tmp = find_motif_instance ( motifs[i], seq, dG_t) ;
    printf("%d\t", tmp) ;
    match+=tmp ;
  }
  printf("%d", match) ;
  
  return (0) ;
}
