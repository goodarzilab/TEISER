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

  char     *seeddir ;
  char     *fastafile1 ;
  char     *fastafile2 ;
  char     *dataoutfile ;
  char     *seedoutfile ;

  s_sequence **sequences1 ;
  s_sequence **sequences2 ;
  s_motif  **seeds ;
  s_motif  **passed ;

  int      seq_count1 ;
  int      seq_count2 ;
  int      seed_cnt = 0 ;
  int      nbseeds = 0 ;

  seeddir          = get_parameter(argc, argv, "-seeddir") ;
  fastafile1       = get_parameter(argc, argv, "-fastafile1") ;
  fastafile2       = get_parameter(argc, argv, "-fastafile2") ;
  dataoutfile      = get_parameter(argc, argv, "-dataoutfile") ;
  seedoutfile      = get_parameter(argc, argv, "-seedoutfile") ;
  if (exist_parameter(argc, argv, "-dG_t")) {
    dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  }

  seq_count1 = read_FASTA ( fastafile1, &sequences1, 0) ;
  seq_count2 = read_FASTA ( fastafile2, &sequences2, 0) ;
  passed = (s_motif**) malloc ( (1000000)*sizeof(s_motif*) ) ;

  FILE* f     = fopen(dataoutfile, "w") ;
  FILE* fseed = fopen(seedoutfile, "wb") ;
  if (!f || !fseed)
    die("Cannot open datafile\n");

  for (i=0 ; i<35 ; i++){
    char seedfile[1000] ;
    sprintf(seedfile, "%s/seeds.4-7.4-9.4-6.14-20.%03d.bin", seeddir, i+1) ;
    printf("seedfile: %s\n", seedfile);
    FILE *fptr = fopen ( seedfile, "rb") ;
    if (!fptr){
      printf("Could not open the seed file...\n") ;
      exit(0) ;
    } 
    seed_cnt = read_motifs( fptr, &seeds ) ;
    printf("%d seeds were loaded...\n", seed_cnt) ;
    fflush(stdout) ;
    fclose(fptr) ;

    for (j=0 ; j<seq_count1 ; j++){
      for (k=0 ; k<seed_cnt ; k++){
        int match1 = find_motif_instance ( seeds[k], sequences1[j], dG_t ) ;
        int match2 = find_motif_instance ( seeds[k], sequences2[j], dG_t ) ;
        //printf("%d\t%d\n", match1, match2) ; getchar() ;
        if ( (match1 == -1 && match2 >= 0) || (match1 >= 0 && match2 == -1) ){
          passed[nbseeds] = copy_motif (seeds[k]) ;
          nbseeds++ ;
          printf("%d\t%d\t%s\t%d\t%d\n", nbseeds, k, sequences1[j]->name, match1, match2) ;
          fprintf(f, "%d\t%d\t%s\t%d\t%d\n", nbseeds, k, sequences1[j]->name, match1, match2) ;
        }
      }
    }
    free(seeds) ;
    fclose(fptr) ;
  }
  
  printf ("%d seeds were written to file...\n", nbseeds) ;
  write_motifs (fseed, passed, nbseeds) ;

  fclose(fseed) ;
  fclose(f) ;
  return (0) ;
}
