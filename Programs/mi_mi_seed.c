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
  int      i ;
  float    dG_t = -1 ;

  char     *seedfile ;
  char     *rna_fastafile ;
  char     *expfile ;
  char     *dataoutfile ;
  
  float    *E ;
  int      *E_q ;
  int      *M_q ;

  struct   my_hsearch_data *h_rna_ind ;
  ENTRY    e ;
  ENTRY*   ep ;
  int      hashret =0 ;

  int      quantized = 1 ;
  int      ebins = 0 ;
  int      mbins = 2 ;
  int      divbins = 50 ;
  float*   E_q_bins = 0 ;
  
  int      shuffle = 1000000 ;
  
  int      rnd_fasta = 0 ;
  float    max_p = 0.0000001 ;
  float    max_z = -1 ;

  int      mincount = 5 ;

  s_sequence **sequences ;
  s_motif  **motifs ;

  char     **seq_names ;
  int      seq_count ;
  int      t_seq_count ;
  int      motif_count = 0 ;

  int      idx_eval           = -1;
  int      idx_eval_up;
  int      idx_eval_do;
  int      fastthreshold_jump = 200;

  int      last_idx_eval;
  int      nb_prev_bad;

  int*     seed_pass;

  s_seed* seed_array ;

  seedfile         = get_parameter(argc, argv, "-seedfile") ;
  rna_fastafile    = get_parameter(argc, argv, "-rna_fastafile") ;
  dataoutfile      = get_parameter(argc, argv, "-dataoutfile") ;
  
  expfile          = get_parameter(argc, argv, "-expfile") ;
  
  if (exist_parameter(argc, argv, "-max_p")) {
    max_p          = atof(get_parameter(argc, argv, "-max_p"));
  }
  if (exist_parameter(argc, argv, "-max_z")) {
    max_z          = atof(get_parameter(argc, argv, "-max_z"));
  }
  if (exist_parameter(argc, argv, "-rnd_fasta")) {
    rnd_fasta      = atoi(get_parameter(argc, argv, "-rnd_fasta"));
  }
  if (exist_parameter(argc, argv, "-ebins")) {
    ebins          = atoi(get_parameter(argc, argv, "-ebins"));
  }
  if (exist_parameter(argc, argv, "-shuffle")) {
    shuffle        = atoi(get_parameter(argc, argv, "-shuffle"));
  }
  if (exist_parameter(argc, argv, "-mincount")) {
    mincount       = atoi(get_parameter(argc, argv, "-mincount"));
  }
  if (exist_parameter(argc, argv, "-dG_t")) {
    dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  }

  FILE *f, *fseed ;
  FILE *fptr = fopen ( seedfile, "rb") ;
  if (!fptr){
    printf("Could not open the seed file...\n") ;
    exit(0) ;
  }

  motif_count = read_motifs( fptr, &motifs ) ;
  printf("%d seeds were loaded...\n", motif_count) ;
  fflush(stdout) ;
  fclose(fptr) ;

  if (rnd_fasta == 1){
    printf("Running in shuffle mode...\n") ;
  }else{
    printf("Running in discovery mode...\n") ;
  }
  t_seq_count = read_FASTA ( rna_fastafile, &sequences, rnd_fasta) ;
  printf("%d sequences loaded...\n", t_seq_count) ;
  fflush(stdout) ;
  
  E = read_expfile (expfile, sequences, t_seq_count, &seq_names, &seq_count) ;
  printf("Expfile loaded: %d values...\n", seq_count) ;
  fflush(stdout) ;
  if ((quantized == 0) && (ebins == 0)) {
    ebins = (int)(0.5 + (float)seq_count / ( divbins * mbins ));
  }
  
  if (quantized == 0) {
    printf("Adding small values...\n") ;
    add_small_values_to_identical_floats(E, seq_count);
  }

  printf("Quantizing the input vector...") ;
  E_q  = (int*)malloc((seq_count) * sizeof(int)) ;
  quantize_E(E, seq_count, quantized, &ebins, &E_q, &E_q_bins);
  printf("Done\n") ;
  fflush(stdout) ;

  h_rna_ind = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
  hashret = my_hcreate_r(100000, h_rna_ind);
  if (hashret == 0) {
    printf("main: couldn't make the hashtable...\n");
    exit(0);
  }
  for (i=0 ; i<seq_count ; i++){
    e.key  = strdup(seq_names[i]) ;
    e.data = (char*) i ;
    hashret = my_hsearch_r(e, ENTER, &ep, h_rna_ind);
    if (hashret == 0){
      printf("main: couldn't add the data to hashtable...\n");
      exit(0);
    }
  }

  seed_array = (s_seed*)malloc((motif_count) * sizeof(s_seed)) ;
  int nbseeds = 0 ;
  int hits =0;
  f     = fopen(dataoutfile, "w") ;
  if (!f){
    printf ("Couldn't: %s", dataoutfile) ;
    exit(0) ;
  }
  for (i=0 ; i<motif_count ; i++){
    if (i%100 == 0 && i>0){
      printf(".") ;
      fflush(stdout) ;
    }
    if (i%10000 == 0 && i>0){
      printf("\n") ;
    }
    int hits ;
    M_q = get_motif_profile (motifs[i], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;

    if (hits <= mincount){
      free(M_q) ;
      continue ;
    }
    float mi = CalculateMIbasic(M_q, E_q, seq_count, mbins, ebins) ;

    fprintf (f, "%f\n", mi) ;
    
    nbseeds++ ;
    free(M_q) ;
  }

  free(E) ;
  fclose(fseed) ;
  fclose(f) ;
  return (0) ;
}
