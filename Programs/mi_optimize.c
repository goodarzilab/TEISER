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
  float    dG_t = -1 ;

  char     *seedfile ;
  char     *rna_fastafile ;
  char     *expfile ;
  char     *dataoutfile ;
  char     *motifoutfile ;
  char     *location ;

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

  int      mincount = 5 ;

  double   minr = 5.0 ;
  double   minratio = 0;
  int      midx;

  float    maxfreq = 0.5;
  float    myfreq;
  float    lastmyfreq;
  float    best_lastmyfreq;

  int      jn   = 10;
  int      jn_t = 6;
  int      jn_f = 3;

  s_sequence **sequences ;
  s_motif  **motifs ;
  s_motif  **opt_motifs ;

  char     **seq_names ;
  int      seq_count ;
  int      t_seq_count ;
  int      motif_count = 0 ;

  seedfile         = get_parameter(argc, argv, "-seedfile") ;
  rna_fastafile    = get_parameter(argc, argv, "-rna_fastafile") ;
  dataoutfile      = get_parameter(argc, argv, "-dataoutfile") ;
  motifoutfile     = get_parameter(argc, argv, "-motifoutfile") ;
  location         = get_parameter(argc, argv, "-location") ;

  expfile          = get_parameter(argc, argv, "-expfile") ;
  quantized        = atoi(get_parameter(argc, argv, "-quantized"));
  
  if (exist_parameter(argc, argv, "-max_p")) {
    max_p          = atof(get_parameter(argc, argv, "-max_p"));
  }
  if (exist_parameter(argc, argv, "-max_z")) {
    max_z          = atof(get_parameter(argc, argv, "-max_z"));
  }
  if (exist_parameter(argc, argv, "-minr")) {
    minr           = atof(get_parameter(argc, argv, "-minr"));
  }
  if (exist_parameter(argc, argv, "-ebins")) {
    ebins          = atoi(get_parameter(argc, argv, "-ebins"));
  }
  if (exist_parameter(argc, argv, "-jn_t")) {
    jn_t           = atoi(get_parameter(argc, argv, "-jn_t"));
  }
  if (exist_parameter(argc, argv, "-shuffle")) {
    shuffle        = atoi(get_parameter(argc, argv, "-shuffle"));
  }
  if (exist_parameter(argc, argv, "-mincount")) {
    mincount       = atoi(get_parameter(argc, argv, "-mincount"));
  }
  if (exist_parameter(argc, argv, "-dooptimize")) {    
    dooptimize     = atoi(get_parameter(argc, argv, "-dooptimize"));
  }

  if (exist_parameter(argc, argv, "-doonlypositive")) {
    doonlypositive = atoi(get_parameter(argc, argv, "-doonlypositive"));
  }
  if (exist_parameter(argc, argv, "-dG_t")) {
    dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  }

  FILE *f, *fmotif ;
  FILE *fptr = fopen ( seedfile, "rb") ;
  if (!fptr){
    printf("Could not open the seed file: %s\n", seedfile) ;
    exit(0) ;
  }

  motif_count = read_motifs( fptr, &motifs ) ;
  printf("%d seeds were loaded...\n", motif_count) ;
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
  hashret = my_hcreate_r(1000000, h_rna_ind);
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
  opt_motifs = (s_motif**) malloc (motif_count*sizeof(s_motif*)) ;
  int *jnres = (int*) malloc (motif_count*sizeof(int)) ;
  for (i=0 ; i<motif_count ; i++){
    M_q = get_motif_profile (motifs[i], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
    
    //  check how much information it adds to the previous guys
    if (opt_count > 0 && minr>0){
      minratio = minCondInfoNormalized(opt_motifs, mbins, opt_count, M_q, mbins, E_q, ebins, seq_count, minr, &midx, sequences, t_seq_count, h_rna_ind, dG_t) ;
    }else{
      minratio = minr+1 ;
    }
    if (minratio < minr){
      free(M_q) ;
      printf("seed %d killed by motif %d (ratio=%f).\n", i, midx,minratio) ;
      continue ;
    }else{
      printf("optimizing.\n") ;
    }

    if (doonlypositive==1){
      float r  = pearson_int (M_q, E_q, seq_count) ;
      if (r<0){
	free(M_q) ;
	printf("seed %d killed due to negative association (pearson=%4.3f)\n", i, r) ;
	continue ;
      }
    }

    lastmyfreq = (float)(hits) / seq_count ;
    best_lastmyfreq = lastmyfreq;

    s_motif *bestmotif = copy_motif(motifs[i]) ;
    if (dooptimize==1){
      //  initial mi value
      init_best_mymi = CalculateMIbasic(M_q, E_q, seq_count, mbins, ebins) ;
      printf("Initial MI = %3.4f\n", init_best_mymi) ;

      // create a random index
      int *k_shu ;
      int *k_inc = (int*) malloc (sizeof(int)*motifs[i]->num_phrases) ;
      int k ;
      for (k=0; k<motifs[i]->num_phrases; k++) {
	k_inc[k] = k;
      }
      k_shu = shuffleInt(k_inc, motifs[i]->num_phrases); 

      // optimize motif
      printf("Optimzing the sequence of motif %d/%d...\n", i+1, motif_count) ;
      float bestmi = init_best_mymi ;
      printf("initial motif (mi = %3.3f): %s\n", bestmi, print_motif_to_char(bestmotif)) ;
      for (k=0 ; k<motifs[i]->num_phrases ; k++){
	int pos = k_shu[k] ;
	s_motif *modified_motifs [15] ;
	modify_base( bestmotif, modified_motifs, pos ) ;
	int j=0 ;
	for (j=0 ; j<15 ; j++){
	  int h ;
	  int *M_q_t = get_motif_profile (modified_motifs[j], sequences, t_seq_count, h_rna_ind, &h, dG_t) ;
	  myfreq = (float)(h) / seq_count ;
	  float tempmi = CalculateMIbasic(M_q_t, E_q, seq_count, mbins, ebins) ;
	  if (tempmi>bestmi && h>10 && (myfreq<maxfreq || myfreq<lastmyfreq)){
	    free(bestmotif->phrases) ;
	    free(bestmotif) ;
	    bestmotif = copy_motif(modified_motifs[j]) ;
	    bestmi = tempmi ;
	    lastmyfreq = myfreq ;
	    printf("new motif (mi = %3.3f): %s\n", bestmi, print_motif_to_char(bestmotif)) ;
	  }
	  free (M_q_t) ;
	}
      }
      
      printf("Elongating the motif...\n") ;
      float premi = bestmi ;
      do{
	s_motif *modified_motifs [46] ;
	elongate_motif( bestmotif, modified_motifs) ;
	int j ;
	for (j=0 ; j<46 ; j++){
	  int h ;
	  int *M_q_t = get_motif_profile (modified_motifs[j], sequences, t_seq_count, h_rna_ind, &h, dG_t) ;
	  float tempmi = CalculateMIbasic(M_q_t, E_q, seq_count, mbins, ebins) ;
	  if (tempmi>bestmi && h>10){
	    free(bestmotif->phrases) ;
	    free(bestmotif) ;
	    bestmotif = copy_motif(modified_motifs[j]) ;
	    bestmi = tempmi ;
	    printf("new motif (mi = %3.3f): %s\n", bestmi, print_motif_to_char(bestmotif)) ;
	  }
	  free(M_q_t) ;
	}
      } while (premi >= bestmi && bestmotif->num_phrases > motifs[i]->num_phrases) ;
    }
    free(M_q) ;
    //printf("%d\t%d\n", t_seq_count, seq_count) ; getchar() ;
    M_q = get_motif_profile (bestmotif, sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
    float mi = CalculateMIbasic(M_q, E_q, seq_count, mbins, ebins) ;
    float z = teiser_z_score_test(mi, M_q, mbins, E_q, ebins, seq_count, 10000) ;
    printf("Final z-score = %4.3f (threshold=%3.3f, mi=%3.7f)\n", z, max_z,mi) ;
    if (z>max_z){
      printf("z-score passed the test\n.") ;
      int jn_test = teiser_jn_max_rank_test(M_q, mbins, E_q, ebins, seq_count, shuffle/10, max_p, jn, jn_f) ;
      if (jn_test>=jn_t){
	opt_motifs[opt_count] = copy_motif(bestmotif) ;
	jnres[opt_count] = jn_test ;
	opt_count++ ;
      }else{
	     printf("robustness=%d/%d\n.", jn_test, jn_t) ;
      }
    }

    free(bestmotif->phrases) ;
    free(bestmotif) ;
    free(M_q) ;
  }

  f      = fopen(dataoutfile, "w") ;
  fmotif = fopen(motifoutfile, "wb") ;
  if (!f || !fmotif)
    die("Cannot open datafile\n");

  fprintf(f, "index\tlocation\tmotif-seq\tmotif_structure\tmi-value\tseq mi-value\tfrequency\tz-score\trobustness\tp-value\n") ;
  for (i=0 ; i<opt_count ; i++){
    int *M_qs = get_motif_profile_seq_only (opt_motifs[i], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
    float mis = CalculateMIbasic(M_qs, E_q, seq_count, mbins, ebins) ;

    M_q = get_motif_profile (opt_motifs[i], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
    float mi = CalculateMIbasic(M_q, E_q, seq_count, mbins, ebins) ;
    double pass = evalSeed(M_q, seq_count, mi, mbins, E_q, ebins, shuffle) ;
    float z = teiser_z_score_test(mi, M_q, mbins, E_q, ebins, seq_count, 10000) ;

    char p_val[100] ;
    if (pass == 0){
      sprintf(p_val, "<%.2e", 1.0/shuffle) ;
    }else{
      sprintf(p_val, "<%.6e", pass) ;
    }

    printf("motif %d/%d\t%s\tmi=%3.6f\thits=%d\t%3.4f\n", i+1, opt_count, print_motif_to_char(opt_motifs[i]), mi, hits, z) ;
    fprintf(f, "%d\t%s\t%s\t%3.6f\t%3.6f\t%3.3f\t%3.3f\t%d/10\t%s\n", i, location, print_motif_to_char(opt_motifs[i]), mi, mis, ((float)hits)/seq_count, z, jnres[i], p_val) ;
    free(M_q) ;
    free(M_qs) ;
  }

  write_motifs (fmotif, opt_motifs, opt_count) ;

  free(E) ;
  free(E_q) ;
  fclose(fmotif) ;
  fclose(f) ;
  return (0) ;
}
