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
  float    dG_t=-1 ;

  char     *seedfile ;
  char     *mifile ;
  char     *rna_fastafile ;
  char     *expfile ;
  char     *dataoutfile ;
  char     *seedoutfile ;

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
  mifile           = get_parameter(argc, argv, "-mifile") ;
  rna_fastafile    = get_parameter(argc, argv, "-rna_fastafile") ;
  dataoutfile      = get_parameter(argc, argv, "-dataoutfile") ;
  seedoutfile       = get_parameter(argc, argv, "-seedoutfile") ;

  expfile          = get_parameter(argc, argv, "-expfile") ;
  quantized        = atoi(get_parameter(argc, argv, "-quantized"));
  
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


  int micnt=0 ;
  int micol=0 ;
  float **midata ;
  char **miindex ;
  char **micolnames ;
  readFloatTable(mifile, &micol, &micnt, &midata, &miindex, &micolnames, 0, 0) ;

  seed_array = (s_seed*)malloc((motif_count) * sizeof(s_seed)) ;
  int nbseeds = 0 ;
  int hits =0;

  for (i=0 ; i<motif_count ; i++){
    if (i%100 == 0 && i>0){
      printf(".") ;
      fflush(stdout) ;
    }
    if (i%10000 == 0 && i>0){
      printf("\n") ;
    }
    //printf("%d\t%f\n", i, midata[i][0]) ; getchar() ;
    seed_array[nbseeds].index = i;
    seed_array[nbseeds].score = midata[i][0] ;
    seed_array[nbseeds].hits  = hits ;
    nbseeds++ ;
  }
  qsort((void*)seed_array, nbseeds, sizeof(s_seed), CmpFunc) ;
  printf("\n") ;
  printf("number of seeds = %d/%d\n", nbseeds, motif_count) ;
  fflush(stdout) ;
  for (i=0 ; i<min(15,nbseeds) ; i++){
    printf("Seed %d, mi = %4.10f\n", seed_array[i].index, seed_array[i].score );
    fflush(stdout);
  }
  
  printf("Determining threshold ...\n");
  seed_pass = (int*)malloc(nbseeds * sizeof(int)) ;
  for (i=0; i<nbseeds; i++) {
    seed_pass[i] = -1; // don't know
  }	
  idx_eval      =  0;
  last_idx_eval = -1;
  // phase 1, find the lower limit     
  for (idx_eval=0; idx_eval<nbseeds; idx_eval+=fastthreshold_jump) {
    // evaluate idx eval
    int hits ;
    M_q = get_motif_profile (motifs[seed_array[idx_eval].index], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
    double pass = evalSeed(M_q, seq_count, seed_array[idx_eval].score, mbins, E_q, ebins, shuffle);
    float z = teiser_z_score_test( seed_array[idx_eval].score, M_q, mbins, E_q, ebins, seq_count, shuffle) ;
    free(M_q) ;
    if (pass>max_p || z<max_z) {
      seed_pass [idx_eval] = 0;
      printf("%d/%d: seed %d didn't pass (p=%3.4f, z=%3.4f)...\n", idx_eval, nbseeds, seed_array[idx_eval].index, pass, z);
      fflush(stdout) ;
      break;
    }else{
      seed_pass [idx_eval] = 1;
      last_idx_eval = idx_eval;
      printf("%d/%d: seed %d passed... (p=%3.4f, z=%3.4f)\n", idx_eval, nbseeds, seed_array[idx_eval].index, pass, z);
      fflush(stdout) ;
    }
  }

  printf("Decreasing intervals phase.\n");
  fflush(stdout) ;
  if ( last_idx_eval >= 0 ) {
    // do line-in-desert
    idx_eval_up  = max(0, idx_eval - fastthreshold_jump);
    idx_eval_do  = min(nbseeds - 1, idx_eval);
    while ((idx_eval_do - idx_eval_up) > 10) {
      idx_eval     = idx_eval_up + (idx_eval_do - idx_eval_up) / 2;
      int hits ;
      M_q = get_motif_profile (motifs[seed_array[idx_eval].index], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
      double pass = evalSeed(M_q, seq_count, seed_array[idx_eval].score, mbins, E_q, ebins, shuffle) ;
      float z = teiser_z_score_test( seed_array[idx_eval].score, M_q, mbins, E_q, ebins, seq_count, shuffle) ;
      free(M_q) ;
      if (pass<=max_p && z>max_z) {
	// go down half interval
	idx_eval_up   = idx_eval;
	last_idx_eval = idx_eval; 
	
	// record that it passed
	seed_pass[ idx_eval] = 1;
	printf("%d: up:%d down:%d, seed %d passed... (p=%3.4f, z=%3.4f)\n", idx_eval, idx_eval_up, idx_eval_do, seed_array[idx_eval].index, pass, z);
	fflush(stdout) ;
      }else{
	// go up half interval
	idx_eval_do  = idx_eval;

	// record that it did not pass
	seed_pass[ idx_eval] = 0;
	printf("%d: seed %d didn't pass (p=%3.4f, z=%3.4f)...\n", idx_eval, seed_array[idx_eval].index, pass, z);
	fflush(stdout) ;
      }
    }
  }

  //  eval next 10 to be sure
  printf("Searches for 10 consecutive 'not passed'.\n");
  fflush (stdout) ;
  
  if (last_idx_eval >= 0) {
    nb_prev_bad   = 0;
  } else {
    nb_prev_bad   = 1;
  }
  
  int check = -1 ;
  float z   = 0 ;
  i = min(nbseeds-1, last_idx_eval + 10);
  while ((i<nbseeds) && (nb_prev_bad<10)) {
    if (seed_pass[i] != -1){
      check = seed_pass[i];
    }else{
      int hits ;
      M_q = get_motif_profile (motifs[seed_array[i].index], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
      double pass = evalSeed(M_q, seq_count, seed_array[i].score, mbins, E_q, ebins, shuffle) ;
      z = teiser_z_score_test( seed_array[idx_eval].score, M_q, mbins, E_q, ebins, seq_count, shuffle) ;
      free(M_q) ;
      if (pass<max_p && z>max_z)
	check=1 ;
      else
	check=0 ;
    }
    
    if (check == 1) {
      // position of the last good one
      last_idx_eval = i;
      // reset nb bad
      nb_prev_bad = 0;
      // store
      seed_pass[i] = 1;
      printf("%d: seed %d passed (z=%f)...\n", i, seed_array[i].index, z);
      fflush(stdout) ;

      // jump 10 down (+1 because the current one is good)
      i += 10 + 1; // add 1 because 1 is going to be removed immediately below

    }else{
      nb_prev_bad ++;
      seed_pass[i] = 0;
      printf("%d: seed %d didn't pass (z=%f)...\n", i, seed_array[i].index, z);
      fflush(stdout) ;
    }
    i--;
  }
  
  f     = fopen(dataoutfile, "w") ;
  fseed = fopen(seedoutfile, "wb") ;
  if (!f || !fseed)
    die("Cannot open datafile\n");
  
  nbseeds=0 ;
  
  fprintf(f, "index\tseed-seq\tseed-struct\tmi-value\tfrequency\n") ;
  s_motif** accepted_motifs = (s_motif**) malloc ( (last_idx_eval+1)*sizeof(s_motif*) ) ;
  for (i=0 ; i<=last_idx_eval ; i++){
    accepted_motifs[nbseeds] = copy_motif (motifs[seed_array[i].index]) ;
    nbseeds++ ;
    
    fprintf(f, "%d\t%s\t%5.4f\t%3.3f\n", seed_array[i].index, print_motif_to_char(motifs[seed_array[i].index]), seed_array[i].score, ((float)seed_array[i].hits)/seq_count ) ;
    fflush(stdout) ;
  }

  write_motifs (fseed, accepted_motifs, nbseeds) ;

  free(E) ;
  fclose(fseed) ;
  fclose(f) ;
  return (0) ;
}
