#include "stdio.h"
#include "stdlib.h"

#include "structures.h"
#include "nucleotides.h"
#include "dataio.h"
#include "sequences.h"
#include "statistics.h"
#include "readFASTA.h"

NUCBIT lcl_get_base_id( char base ) {
  NUCBIT base_id = 1;

  int i=0 ;
  for( i=0 ; i<NUM_LETTERS ; i++ )
    if( base == lcl_base_templates1[i] ||
	base == lcl_base_templates2[i] ||
	base == lcl_base_templates3[i] ||
	base == lcl_base_templates4[i] )
      return base_id;
    else
      base_id <<= 1;

  if( base == lcl_N1 || base == lcl_N2 )
    return _N;
  
  return 0;
}

int read_FASTA( char* fastafile, s_sequence ***sequences, int shuffle ) {
  int       i ;
  seqI      si ;
  char*     seq ;
  int       size ;
  char*     name ;
  int       nseqs = 0 ;
  int       seqlen = 0 ;
  int       max_nseq = 1000000 ;
  	
  if (seqI_open(&si, fastafile) == 0) {
    printf("Error opening fastafile: %s\n", fastafile) ;
    exit(0) ;
  }

  *sequences = (s_sequence**) malloc (max_nseq*sizeof(s_sequence*)) ;
  while ((seq = seqI_nextSequence(&si, &name, &size))) {
    if (shuffle == 1){
      char *vnew = (char*) malloc ((strlen(seq)+1) * sizeof(char)) ;
      for (i=0; i<strlen(seq); i++) {
	int r       = (int)((double)i*default_rand());
	vnew[i] = vnew[r];
	vnew[r] = seq[i];
      }
      vnew[strlen(seq)] = '\0' ;
      strcpy(seq, vnew) ;
      free(vnew) ;
    }
    seqlen = strlen(seq) ;
    (*sequences) [nseqs] = (s_sequence*) malloc (sizeof(s_sequence)) ;
    (*sequences) [nseqs]->bases = (NUCBIT*) malloc (sizeof(NUCBIT)*seqlen) ;
    (*sequences) [nseqs]->name = strdup(name) ;
    (*sequences) [nseqs]->length = seqlen ;
    for( i=0 ; i<seqlen ; i++ ) {
      (*sequences) [nseqs]->bases[i] = lcl_get_base_id (seq[i]);
      if( !((*sequences) [nseqs]->bases[i]) ) {
	// unknown letter in sequence
	(*sequences) [nseqs]->bases[i] = _N;
      }
    }
    nseqs++ ;
    free(name) ;
    //free(seq) ;
  }

  return nseqs ;
}
