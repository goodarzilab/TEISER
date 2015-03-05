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

int* readIntSet(char* setfile, int* n) ;

int main(int argc, char ** argv) {
  int      i ;

  char     *motiffile ;
  char     *indexfile ;
  char     *outfile ;

  int      *indices ;
  s_motif  **motifs ;
  s_motif  **opt_motifs ;
  int      motif_count = 0 ;
  int      index_count = 0 ;

  motiffile        = get_parameter(argc, argv, "-motiffile") ;
  indexfile        = get_parameter(argc, argv, "-indexfile") ;
  outfile          = get_parameter(argc, argv, "-outfile") ;

  FILE *f, *fmotif ;
  FILE *fptr = fopen ( motiffile, "rb") ;
  if (!fptr){
    printf("Could not open the seed file: %s\n", motiffile) ;
    exit(0) ;
  }

  motif_count = read_motifs( fptr, &motifs ) ;
  printf("%d seeds were loaded...\n", motif_count) ;
  fflush(stdout) ;
  fclose(fptr) ;

  indices = readIntSet(indexfile, &index_count) ;

  opt_motifs = (s_motif**) malloc (index_count*sizeof(s_motif*)) ;
  
  for (i=0 ; i<index_count ; i++){
    printf("%d\t%d\t%s\n", i, indices[i], print_motif_to_char(motifs[indices[i]])) ;
    opt_motifs[i] = copy_motif(motifs[indices[i]]) ;
  }

  fmotif = fopen(outfile, "wb") ;
  if (!fmotif)
    die("Cannot open motif outfile\n");

  write_motifs (fmotif, opt_motifs, index_count) ;
  fclose(fmotif) ;

  return (0) ;
}

int* readIntSet(char* setfile, int* n) 
{

   // ok now get the set of ORFS
   
  FILE*  fp;
  int    nb = 0;
  char*  buff;
  int* oset;
  
 
    
  nb   = nbLinesInFile(setfile);
  
  *n   = nb;
 

  buff = (char*)malloc(10000 * sizeof(char));
  fp = fopen(setfile, "r");
  if (!fp) {
    printf("could not open set data %s\n", setfile);
  }
 
  // alloc memory
  oset = (int*)malloc(nb * sizeof(int));
  if (!oset) {
    printf("could not allocate memory for oset\n");
    exit(0);
  }

  nb = 0;
  while (!feof(fp)) {
    fscanf(fp, "%s\n", buff);
    oset[nb] = atoi(buff);
    nb++;
  }
  
  free(buff) ;
  fclose(fp);
  return oset;
}
