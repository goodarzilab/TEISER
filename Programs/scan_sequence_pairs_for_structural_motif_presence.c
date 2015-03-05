//./scan_sequence_pairs_for_structural_motif_presence -seedfile /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Results/SNPs_in_3utrs_MDA_LM2_motif_filtered_structural_seeds_non_redundant.bin -fastafile1 /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Results/SNPs_in_3utrs_MDA_LM2_cons_50bp_center_MDA.txt -fastafile2 /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Results/SNPs_in_3utrs_MDA_LM2_cons_50bp_center_LM2.txt -dataoutfile /Users/hani/Life/Projects/Ongoing_Projects/Metatstaic_SNPs/Results/SNPs_in_3utrs_MDA_LM2_motif_filtered_structural_seeds_non_redundant_report.txt

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
  char     *fastafile1 ;
  char     *fastafile2 ;
  char     *dataoutfile ;
  
  s_sequence **sequences1 ;
  s_sequence **sequences2 ;
  s_motif  **seeds ;
  
  int      seq_count1 ;
  int      seq_count2 ;
  int      seed_cnt = 0 ;
  
  seedfile         = get_parameter(argc, argv, "-seedfile") ;
  fastafile1       = get_parameter(argc, argv, "-fastafile1") ;
  fastafile2       = get_parameter(argc, argv, "-fastafile2") ;
  dataoutfile      = get_parameter(argc, argv, "-dataoutfile") ;
  if (exist_parameter(argc, argv, "-dG_t")) {
    dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  }
  
  seq_count1 = read_FASTA ( fastafile1, &sequences1, 0) ;
  seq_count2 = read_FASTA ( fastafile2, &sequences2, 0) ;
  
  FILE* f     = fopen(dataoutfile, "w") ;
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

  for (j=0 ; j<seq_count1 ; j++){
    for (k=0 ; k<seed_cnt ; k++){
      int match1 = find_motif_instance ( seeds[k], sequences1[j], dG_t ) ;
      int match2 = find_motif_instance ( seeds[k], sequences2[j], dG_t ) ;
      if ( (match1 == -1 && match2 >= 0) || (match1 >= 0 && match2 == -1) ){
        char *seq1, *seq2 ;
        seq1 = (char*) malloc ((sequences1[j]->length+1) * sizeof(char)) ;
        seq2 = (char*) malloc ((sequences2[j]->length+1) * sizeof(char)) ;
        for (i=0 ; i<sequences1[j]->length ; i++){
          switch (sequences1[j]->bases[i]){
            case _U:
              seq1[i] = 'T' ;
              break ;
            case _C:
              seq1[i] = 'C' ;
              break ;
            case _G:
              seq1[i] = 'G' ;
              break ;
            case _A:
              seq1[i] = 'A' ;
              break ;
          }
        }
        seq1[sequences1[j]->length] = '\0' ;
        for (i=0 ; i<sequences2[j]->length ; i++){
          switch (sequences2[j]->bases[i]){
            case _U:
              seq2[i] = 'T' ;
              break ;
            case _C:
              seq2[i] = 'C' ;
              break ;
            case _G:
              seq2[i] = 'G' ;
              break ;
            case _A:
              seq2[i] = 'A' ;
              break ;
          }
        }
        seq2[sequences2[j]->length] = '\0' ;
        printf("%s\t%s\t%d\t%d\t%s\t%s\n", print_motif_to_char(seeds[k]), sequences1[j]->name, match1, match2, seq1, seq2) ;
        fprintf(f, "%s\t%s\t%d\t%d\t%s\t%s\n", print_motif_to_char(seeds[k]), sequences1[j]->name, match1, match2, seq1, seq2) ;
      }
    }
  }
  free(seeds) ;
  fclose(fptr) ;

  fclose(f) ;
  return (0) ;
}
