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
#include "dataio.h"
#include "read_write_motif.h"
#include "teiser_functions.h"

int main(int argc, char ** argv) {
  int      i,j ;

  char     *maindir ;
  char     *metadir ;
  char     *filename ;
  char     *extname ;
  char     *rna_fastafile ;
  char     *expfile ;
  char     *dataoutfile ;
  char     *motifoutfile ;
  char     *suffix ;

  int      nseedfiles ;
  int      mode = 1 ;

  s_motif  **motifs ;
  s_motif  **opt_motifs ;

  int      motif_count = 0 ;

  maindir          = get_parameter(argc, argv, "-maindir") ;
  metadir          = get_parameter(argc, argv, "-metadir") ;
  filename         = get_parameter(argc, argv, "-filename") ;
  extname          = get_parameter(argc, argv, "-extname") ;
  suffix           = get_parameter(argc, argv, "-suffix") ;
  motifoutfile     = get_parameter(argc, argv, "-motifoutfile") ;

  nseedfiles       = atoi(get_parameter(argc, argv, "-nseedfiles"));
  mode             = atoi(get_parameter(argc, argv, "-mode"));
  
  motifs = (s_motif**) malloc (10000*sizeof(s_motif*)) ;
  s_motif_list *motif_list = (s_motif_list*) malloc (10000*sizeof(s_motif_list)) ;
  for (i=0 ; i<nseedfiles ; i++){
    char S[10000] ;
    char fn[10000] ;
    FILE *fptr ;
    char summary[10000] ;

    if (mode == 1){
      sprintf (S, "%s_s%d.%s", filename, i, extname) ;
      sprintf (fn, "%s/%s%s/%s.optim.bin", maindir, S, suffix, S) ;
      fptr = fopen ( fn, "rb") ;
      if (!fptr){
	printf("Could not open the motif file %d: %s...\n", i, fn) ;
	continue ;
      }
      sprintf (summary, "%s/%s%s/%s.summary", maindir, S, suffix, S) ;
    }else if (mode == 2){
      if (i==0){
	sprintf (S, "%s/UP/%s.%s.optim.bin", metadir, filename, extname) ;
      }else{
	sprintf (S, "%s/DN/%s.%s.optim.bin", metadir, filename, extname) ;
      }
      fptr = fopen ( fn, "rb") ;
      if (!fptr){
	printf("Could not open the motif file: %s...\n", fn) ;
	continue ;
      }
      sprintf (summary, "%s/DN/%s.%s.summary", metadir, filename, extname) ;
    }
    s_motif **temp_motifs ;
    int temp_count = read_motifs( fptr, &temp_motifs ) ;
    printf("%d seeds were loaded...\n", motif_count) ;
    
    char ***summarytable ;
    int n, m ;
    readStringTable(summary, &summarytable, &n, &m) ;

    for (j=0 ; j<temp_count ; j++){
      motifs[motif_count] = copy_motif(temp_motifs[j]) ;
      motif_list[motif_count].index = motif_count ;
      printf("%s\t%s\t%s\n", summarytable[j+1][2], summarytable[j+1][3], summarytable[j+1][4]) ;
      motif_list[motif_count].score = atof(summarytable[j+1][4]) ;
      motif_count++ ;
    }
    fclose(fptr) ;
  }

  qsort((void*)motif_list, motif_count, sizeof(s_motif_list), CmpFunc) ;
  FILE *fmotif ;

  int opt_count=0 ;
  opt_motifs = (s_motif**) malloc (motif_count*sizeof(s_motif*)) ;
  for (i=0 ; i<motif_count ; i++){
    int id = motif_list[i].index ;
    opt_motifs[opt_count++] = copy_motif(motifs[id]) ;
  }

  fmotif = fopen(motifoutfile, "wb") ;
  if (!fmotif)
    die("Cannot open motiffile\n");

  write_motifs (fmotif, opt_motifs, opt_count) ;

  fclose(fmotif) ;
  return (0) ;
}
