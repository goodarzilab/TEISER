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
  float    dG_t = -1 ;

  char     *motiffile ;
  char     *fastafile1 ;
  char     *fastafile2 ;
  char     *homologyfile ;
  char     *consfile ;
  
  struct   my_hsearch_data *h_org1 ;
  struct   my_hsearch_data *h_org2 ;
  ENTRY    e ;
  ENTRY*   ep ;
  int      hashret =0 ;
  
  s_sequence **sequences_1 ;
  s_sequence **sequences_2 ;
  s_motif  **motifs ;

  char     **seq_names_1, **seq_names_2 ;
  int      seq_count_1, seq_count_2 ;
  int      t_seq_count_1, t_seq_count_2 ;
  int      motif_count ;

  int      *M_q_1 ;
  int      *M_q_2 ;
  int      shuffle = 10000 ;

  char     ***homologyTable ;
  int      num_homologs ;
  int      num_orgs ;
  
  motiffile        = get_parameter(argc, argv, "-motiffile") ;
  fastafile1       = get_parameter(argc, argv, "-fastafile1") ;
  fastafile2       = get_parameter(argc, argv, "-fastafile2") ;
  homologyfile     = get_parameter(argc, argv, "-homologyfile") ;
  consfile         = get_parameter(argc, argv, "-consfile") ;
  if (exist_parameter(argc, argv, "-dG_t")) {
    dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  }

  FILE *fptr = fopen ( motiffile, "rb") ;
  if (!fptr){
    printf("Could not open the seed file...\n") ;
    exit(0) ;
  }

  motif_count = read_motifs( fptr, &motifs ) ;
  printf("%d seeds were loaded...\n", motif_count) ;
  fflush(stdout) ;
  fclose(fptr) ;

  int rnd_fasta=0 ;
  t_seq_count_1 = read_FASTA ( fastafile1, &sequences_1, rnd_fasta) ;
  printf("%d sequences loaded from organism 1...\n", t_seq_count_1) ;
  fflush(stdout) ;

  t_seq_count_2 = read_FASTA ( fastafile2, &sequences_2, rnd_fasta) ;
  printf("%d sequences loaded from organism 2...\n", t_seq_count_2) ;
  fflush(stdout) ;

  readStringTable(homologyfile, &homologyTable, &num_homologs, &num_orgs) ;
  printf("%d homologous pairs loaded...\n", num_homologs) ;
  fflush(stdout) ;

  h_org1 = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
  hashret = my_hcreate_r(100000, h_org1);
  if (hashret == 0) {
    printf("main: couldn't make the hashtable...\n");
    exit(0);
  }
  for (i=0 ; i<num_homologs ; i++){
    e.key  = strdup(homologyTable[i][0]) ;
    e.data = (char*) i ;
    hashret = my_hsearch_r(e, ENTER, &ep, h_org1);
    if (hashret == 0){
      printf("main: couldn't add the data to hashtable...\n");
      exit(0);
    }
  }

  h_org2 = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
  hashret = my_hcreate_r(100000, h_org2);
  if (hashret == 0) {
    printf("main: couldn't make the hashtable...\n");
    exit(0);
  }
  for (i=0 ; i<num_homologs ; i++){
    e.key  = strdup(homologyTable[i][1]) ;
    e.data = (char*) i ;
    hashret = my_hsearch_r(e, ENTER, &ep, h_org2);
    if (hashret == 0){
      printf("main: couldn't add the data to hashtable...\n");
      exit(0);
    }
  }

  FILE *f = fopen (consfile, "w") ;
  if (!f){
    printf("Could not open the conservation file...\n") ;
    exit(0) ;
  }
  fprintf(f, "id\tMI\tp-value\tp-value(hyper-geom)\n") ;
  for (i=0 ; i<motif_count ; i++){
    int hits1, hits2 ;
    M_q_1 = get_motif_profile (motifs[i], sequences_1, t_seq_count_1, h_org1, &hits1, dG_t) ;
    M_q_2 = get_motif_profile (motifs[i], sequences_2, t_seq_count_2, h_org2, &hits2, dG_t) ;
    int overlap=0 ;
    for (j=0 ; j<num_homologs ; j++){
      if (M_q_1[j]==1 && M_q_2[j]==1)
	overlap++ ;
    }
    double pu = log10(cumhyper(overlap, hits1, hits2, num_homologs));
    float mi = CalculateMIbasic(M_q_1, M_q_2, num_homologs, 2, 2) ;
    double pass = evalSeed(M_q_1, num_homologs, mi, 2, M_q_2, 2, shuffle) ;
    printf("%d\t%f\t%f\t%f\n", i, mi, pass, pu) ;
    fprintf(f,"%d\t%f\t%f\t%f\n", i, mi, pass, pu) ;
  }
  
  fclose(f) ;
  fclose(fptr) ;
  return (0) ;
}
