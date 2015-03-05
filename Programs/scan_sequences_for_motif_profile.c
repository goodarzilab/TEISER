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
  
  seedfile         = get_parameter(argc, argv, "-seedfile") ;
  fastafile        = get_parameter(argc, argv, "-fastafile") ;
  dataoutfile      = get_parameter(argc, argv, "-dataoutfile") ;
  reportfile      = get_parameter(argc, argv, "-reportfile") ;
  if (exist_parameter(argc, argv, "-dG_t")) {
    dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  }
  
  seq_count = read_FASTA ( fastafile, &sequences, 0) ;
  
  FILE* f     = fopen(dataoutfile, "w") ;
  FILE* freport = fopen ( reportfile, "w") ;
  if (!f)
    die("Cannot open datafile\n");

  FILE *fptr = fopen ( seedfile, "rb") ;
  if (!fptr){
    printf("Could not open the seed file...\n") ;
    exit(0) ;
  } 
  seed_cnt = read_motifs( fptr, &seeds ) ;
  printf("%d seeds were loaded...\n", seed_cnt) ;
  fflush(stdout) ;
  fclose(fptr) ;

  for (j=0 ; j<seq_count ; j++){
    //fprintf(f, ">%s\n", sequences[j]->name) ;
    int l = sequences[j]->length ;
    int *stat = (int*) calloc (l, sizeof(int)) ;
    for (k=0 ; k<seed_cnt ; k++){
      int n_inst=0 ;
      int *inst = find_motif_instance_positions ( seeds[k], sequences[j], &n_inst, dG_t) ;
      int m=0 ;
      if (n_inst>0){
        //printf("%s\t%d\n", print_motif_to_char(seeds[k]), match) ;
        fprintf(f, "%s\t%s", sequences[j]->name, print_motif_to_char(seeds[k])) ;
        for (m=0 ; m<n_inst ; m++){
          fprintf(f, "\t%d", inst[m]) ;
          fprintf(freport, ">%s_%d_%d\n", sequences[j]->name, m, inst[m]) ;
          int x ;
          for (x=inst[m]-10 ; x<inst[m]+seeds[k]->linear_length+10 ; x++){
            switch (sequences[j]->bases[x]){
            case _U:
              //printf("T");
              fprintf(freport, "T");
              break ;
            case _A:
              //printf("A");
              fprintf(freport, "A");
              break ;
            case _G:
              //printf("G");
              fprintf(freport, "G");
              break ;
            case _C:
              //printf("C");
              fprintf(freport, "C");
              break ;
            default:
              //printf("N");
              fprintf(freport, "N");
          }
          }
          fprintf(freport, "\n");
        }
        fprintf(f, "\n") ;
      }
      free(inst) ;
    }
  }
  fclose(fptr) ;

  fclose(f) ;
  return (0) ;
}
