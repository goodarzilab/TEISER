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
  char     *motiffile ;
  char     *rna_fastafile ;
  char     *expfile ;
  char     *dataoutfile ;
  char     *motifoutfile ;
  
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
  s_motif  **seeds ;
  s_motif  **motifs ;
  s_motif  **opt_motifs ;

  char     **seq_names ;
  int      seq_count ;
  int      t_seq_count ;
  int      seed_count = 0 ;
  int      motif_count = 0 ;

  motiffile        = get_parameter(argc, argv, "-motiffile") ;
  seedfile         = get_parameter(argc, argv, "-seedfile") ;
  rna_fastafile    = get_parameter(argc, argv, "-rna_fastafile") ;
  dataoutfile      = get_parameter(argc, argv, "-dataoutfile") ;
  motifoutfile     = get_parameter(argc, argv, "-motifoutfile") ;
  
  expfile          = get_parameter(argc, argv, "-expfile") ;
  quantized        = atoi(get_parameter(argc, argv, "-quantized"));
  
  if (exist_parameter(argc, argv, "-minr")) {
    minr           = atof(get_parameter(argc, argv, "-minr"));
  }
  if (exist_parameter(argc, argv, "-ebins")) {
    ebins          = atoi(get_parameter(argc, argv, "-ebins"));
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
  seed_count = read_motifs( fptr, &seeds ) ;
  printf("%d seeds were loaded...\n", seed_count) ;
  fflush(stdout) ;
  fclose(fptr) ;

  fptr = fopen ( motiffile, "rb") ;
  if (!fptr){
    printf("Could not open the motiffile file: %s\n", motiffile) ;
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

  for (i=0 ; i<motif_count ; i++){
    int opt_count=0 ;
    opt_motifs = (s_motif**) malloc ((seed_count+1)*sizeof(s_motif*)) ;
    opt_motifs[opt_count++] = copy_motif(motifs[i]) ;
    
    for (j=0 ; j<seed_count ; j++){
      double ratio = CondInfoRatio(seeds[j], motifs[i], mbins, mbins, E_q, ebins, seq_count, sequences, t_seq_count, h_rna_ind, dG_t) ;
      if (ratio < minr){
        printf("seed %d killed by motif %d (ratio=%f).\n", i, j, ratio) ; fflush(stdout) ;
        opt_motifs[opt_count++] = copy_motif(seeds[j]) ;
      }
    }
    char fn [1000] ;
    sprintf(fn, "%s_%i.seeds.bin", motifoutfile, i) ;
    printf("writing motif %d family to %s...\n", i, fn) ; fflush(stdout) ;
    FILE *fmotif =fopen ( fn, "wb") ;
    if (!fmotif){
      printf("Could not open the motiffile file: %s\n", fn) ;
      exit(0) ;
    }
    write_motifs (fmotif, opt_motifs, opt_count) ;
    printf("file written\n", i, fn) ; fflush(stdout) ;
    for (j=0 ; j<opt_count ; j++){
      free(opt_motifs[j]->phrases) ;
      free(opt_motifs[j]) ;
    }

  }

  return (0) ;
}
