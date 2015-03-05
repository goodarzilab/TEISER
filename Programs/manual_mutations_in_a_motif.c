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

  char     *seedfile ;
  char     *rna_fastafile ;
  char     *expfile ;
  char     *motifoutfile ;

  int      dooptimize = 1 ;
  int      doonlypositive = 0 ;

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
  float    max_z = -100 ;

  float    maxfreq = 0.5;
  float    myfreq;
  float    lastmyfreq;

  s_sequence **sequences ;
  s_motif  **motifs ;

  char     **seq_names ;
  int      seq_count ;
  int      t_seq_count ;
  int      motif_count = 0 ;

  seedfile         = get_parameter(argc, argv, "-seedfile") ;
  rna_fastafile    = get_parameter(argc, argv, "-rna_fastafile") ;
  motifoutfile     = get_parameter(argc, argv, "-motifoutfile") ;

  expfile          = get_parameter(argc, argv, "-expfile") ;
  quantized        = atoi(get_parameter(argc, argv, "-quantized"));
  
  if (exist_parameter(argc, argv, "-ebins")) {
    ebins          = atoi(get_parameter(argc, argv, "-ebins"));
  }
  if (exist_parameter(argc, argv, "-shuffle")) {
    shuffle        = atoi(get_parameter(argc, argv, "-shuffle"));
  }

  FILE *fmotif ;
  FILE *fptr = fopen ( seedfile, "rb") ;
  if (!fptr){
    printf("Could not open the seed file: %s\n", seedfile) ;
    exit(0) ;
  }

  motif_count = read_motifs( fptr, &motifs ) ;
  printf("%d seeds were loaded...\n", motif_count) ;
  s_motif *motif = copy_motif(motifs[0]) ;
  fflush(stdout) ;
  fclose(fptr) ;

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

  int opt_count=0 ;
  int hits =0;
  float init_best_mymi = 0 ;
  
  M_q = get_motif_profile (motif, sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
  init_best_mymi = CalculateMIbasic(M_q, E_q, seq_count, mbins, ebins) ;
  float z = teiser_z_score_test(init_best_mymi, M_q, mbins, E_q, ebins, seq_count, 10000) ;

  print_cfg(motif) ;
  printf("\nInitial MI = %3.4f (z=%3.4f)\n\n", init_best_mymi, z) ;
  s_motif *n_motif = (s_motif*) malloc (sizeof(s_motif)) ;
  n_motif->num_phrases = motif->num_phrases-5 ;
  n_motif->linear_length = motif->linear_length-5 ;
  n_motif->phrases = (s_phrase*) malloc (n_motif->num_phrases*sizeof(s_phrase)) ;
  for (i=0 ; i< n_motif->num_phrases ; i++){
    n_motif->phrases[i].base = motif->phrases[i].base ;                                                                                                                                                                                                                                               n_motif->phrases[i].structure = motif->phrases[i].structure ;
  }
  /*int C=3 ;
  int D=5 ;
  for (i=0 ; i<C ; i++){
    n_motif->phrases[i].base = motif->phrases[i].base ;
    n_motif->phrases[i].structure = motif->phrases[i].structure ;
  }
  for ( i=motif->num_phrases ; i>D ; i--){
    n_motif->phrases[i+3].base = motif->phrases[i].base ;
    n_motif->phrases[i+3].structure = motif->phrases[i].structure ;
    }*/
  n_motif->phrases[3].base = _N ;
  /*n_motif->phrases[C].base = _N ;
  n_motif->phrases[C].structure = _leftBulge ;
  n_motif->phrases[C+1].base = _N ;
  n_motif->phrases[C+1].structure = _leftBulge ;
  n_motif->phrases[C+2].base = _N ;
  n_motif->phrases[C+2].structure = _leftBulge ;
  n_motif->phrases[C+3].base = _N ;
  n_motif->phrases[C+3].structure = _rightBulge ;
  n_motif->phrases[C+4].base = _N ;
  n_motif->phrases[C+4].structure = _rightBulge ;
  n_motif->phrases[C+5].base = _N ;
  n_motif->phrases[C+5].structure = _rightBulge ;*/

  M_q = get_motif_profile (n_motif, sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
  float mi = CalculateMIbasic(M_q, E_q, seq_count, mbins, ebins) ;
  z = teiser_z_score_test(mi, M_q, mbins, E_q, ebins, seq_count, 10000) ;
  print_cfg(n_motif) ;
  printf("\nsecondary MI = %3.4f (z = %3.4f)\n\n", mi,z) ;
  

  fmotif = fopen(motifoutfile, "wb") ;
  if (!fmotif)
    die("Cannot open datafile\n");

  fwrite( 1, sizeof(int), 1, fptr );
  lcl_write_motif(fptr, n_motif) ;

  free(E) ;
  free(E_q) ;
  fclose(fmotif) ;
  return (0) ;
}
