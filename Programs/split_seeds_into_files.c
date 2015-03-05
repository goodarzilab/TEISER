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
  int      i,j,k ;
  
  int      nfiles =0;
  char*    buff = (char*)malloc(100000 * sizeof(char));;
  int      motif_count =0 ;
  int      cntperfile ;
  char     *motiffile ;
  char     **files = (char**)malloc(1000 * sizeof(char*)); ;
  char     *outfile ;

  s_motif  **motifs ;
  motifs = (s_motif**) malloc (10000000*sizeof(s_motif*)) ;

  motiffile        = get_parameter(argc, argv, "-motiffile") ;
  outfile          = get_parameter(argc, argv, "-outfile") ;
  cntperfile       = atoi(get_parameter(argc, argv, "-count_per_file")) ;

  FILE *fptr = fopen ( motiffile, "rb") ;
    
  if (!fptr){
     printf("Could not open the seed file: %s\n", motiffile) ;
     exit(0) ;
  }
  int count = read_motifs( fptr, &motifs ) ;
  printf("%d seeds loaded\n", count) ;
  
  for (i=0 ; i<=count/cntperfile ; i++){
    int b = i*cntperfile ;
    int e = i*cntperfile+cntperfile ;
    if (e>=count){
      e = count ;
    }
    s_motif  **opt_motifs = (s_motif**) malloc (cntperfile*sizeof(s_motif*)) ;
    for (j=b,k=0 ; j<e ; j++,k++){
      opt_motifs[k] = copy_motif(motifs[j]) ;
    }

    char fn [10000];
    sprintf( fn, "%s.sRSM1.%3.3i.bin", outfile, i);
    FILE *fmotif ;
    fmotif = fopen(fn, "wb") ;
    if (!fmotif)
      die("Cannot open datafile\n");
    
    write_motifs (fmotif, opt_motifs, e-b) ;

    for (j=b,k=0 ; j<e ; j++,k++){
      free(opt_motifs[k]->phrases) ;
      free(opt_motifs[k]) ;
    }
    fclose(fmotif) ;
    printf("%i\t%i\t%i\t%s\n",i, i*cntperfile, e-b,fn) ;
    fflush(stdout) ;
  }


  /*FILE*  fp = fopen(motiffile, "rt");
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
  fclose(fmotif) ;*/

  return (0) ;
}
