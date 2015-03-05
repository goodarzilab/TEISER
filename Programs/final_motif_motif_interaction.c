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

int sign(float v) ;

int main(int argc, char ** argv) {
  int      i,j ;

  char     *motifdir_up ;
  char     *motifdir_dn ;
  char     *summaryfile_up ;
  char     *summaryfile_dn ;
  char     *matrixfile ;
  char     *summaryfile ;

  int      motif_cnt_up ;
  int      motif_cnt_dn ;

  float    **mimatrix ;

  int      *M_q1 ;
  int      *M_q2 ;

  struct   my_hsearch_data *h_tmp_ind ;
  struct   my_hsearch_data *h_ind ;
  struct   my_hsearch_data *h_MI ;
  ENTRY    e ;
  ENTRY    e1 ;
  ENTRY*   ep ;
  ENTRY*   ep1 ;
  int      hashret =0 ;

  int      quantized = 1 ;
  int      mbins = 2 ;
  
  int      shuffle = 10000 ;
  
  float    max_p = 0.001 ;

  FILE     *f ;

  motifdir_up      = get_parameter(argc, argv, "-motifdir_up") ;
  motifdir_dn      = get_parameter(argc, argv, "-motifdir_dn") ;
  summaryfile_up   = get_parameter(argc, argv, "-summaryfile_up") ;
  summaryfile_dn   = get_parameter(argc, argv, "-summaryfile_dn") ;
  matrixfile       = get_parameter(argc, argv, "-matrixfile") ;
  summaryfile      = get_parameter(argc, argv, "-summaryfile") ;
  
  if (exist_parameter(argc, argv, "-max_p")) {
    max_p          = atof(get_parameter(argc, argv, "-max_p"));
  }
  max_p = log10(max_p) ;
  if (exist_parameter(argc, argv, "-shuffle")) {
    shuffle        = atoi(get_parameter(argc, argv, "-shuffle"));
  }

  int msum, nsum ;
  char ***summary1, ***summary2 ;
  readStringTable(summaryfile_up, &summary1, &nsum, &msum) ;
  motif_cnt_up = nsum-1 ;
  h_MI = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
  hashret = my_hcreate_r(100000, h_MI);
  for (i=0 ; i<motif_cnt_up ; i++){
    e.key = strdup(summary1[i+1][0]);
    e.data = strdup(summary1[i+1][4]);
    hashret = my_hsearch_r(e, ENTER, &ep, h_MI) ;
    if (hashret == 0) {
      printf("Could not enter entry into hash table (table full?) ...\n");
      exit(0);
    }
  }
  readStringTable(summaryfile_dn, &summary2, &nsum, &msum) ;
  motif_cnt_dn = nsum-1 ;
  for (i=0 ; i<motif_cnt_dn ; i++){
    char ind[100] ;
    sprintf(ind, "%d", motif_cnt_up+atoi(summary2[i+1][0])) ;
    e.key = strdup(ind);
    e.data = strdup(summary2[i+1][4]);
    hashret = my_hsearch_r(e, ENTER, &ep, h_MI) ;
    if (hashret == 0) {
      printf("Could not enter entry into hash table (table full?) ...\n");
      exit(0);
    }
  }

  int n,m ;
  int **data ;
  int gene_count =0 ;
  char **gene_names ;
  char **col_names ;
  if (motif_cnt_up>0){
    char fn[1000] ;
    sprintf(fn, "%s/0.txt", motifdir_up) ;
    readIntTable(fn, &m, &n, &data, &gene_names, &col_names) ;
    gene_count=n ;
  }else{
    char fn[1000] ;
    sprintf(fn, "%s/0.txt", motifdir_dn) ;
    readIntTable(fn, &m, &n, &data, &gene_names, &col_names) ;
  }

  if (motif_cnt_up>0 && motif_cnt_dn>0){
    h_tmp_ind = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
    hashret = my_hcreate_r(100000, h_tmp_ind);
    if (hashret == 0) {
      printf("main: couldn't make the hashtable...\n");
      exit(0);
    }
    for (i=0 ; i<n ; i++){
      e.key  = strdup(gene_names[i]) ;
      e.data = (char*) i ;
      hashret = my_hsearch_r(e, ENTER, &ep, h_tmp_ind);
      if (hashret == 0){
	printf("main: couldn't add the data to hashtable...\n");
	exit(0);
      }
    }

    char fn[1000] ;
    sprintf(fn, "%s/0.txt", motifdir_dn) ;

    for (i=0 ; i<n ; i++){
      free(data[i]) ;
    }
    free(gene_names) ;
    free(col_names) ;
    free(data) ;
    
    readIntTable(fn, &m, &n, &data, &gene_names, &col_names) ;

    char **gene_names_fileterd = (char**) malloc (n * sizeof(char*)) ;
    int gene_cnt_filtered =0;
    h_ind = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
    hashret = my_hcreate_r(100000, h_ind);
    for (i=0 ; i<n ; i++){
      e.key  = strdup(gene_names[i]) ;
      my_hsearch_r(e, FIND, &ep, h_tmp_ind);
      if (ep) {
	gene_names_fileterd[gene_cnt_filtered] = strdup(gene_names[i]) ;

	e1.key  = strdup(gene_names[i]) ;
	e1.data = (char*) gene_cnt_filtered ;
	hashret = my_hsearch_r(e1, ENTER, &ep1, h_ind);
	if (hashret == 0){
	  printf("main: couldn't add the data to hashtable...\n");
	  exit(0);
	}
	gene_cnt_filtered++ ;
      }
    }
    free(gene_names) ;
    gene_names = gene_names_fileterd ;
    gene_count = gene_cnt_filtered ;
  }

  int **profile = (int**) malloc (gene_count*sizeof(int*)) ;
  for (i=0 ; i<gene_count ; i++){
    profile[i] = (int*) malloc ((motif_cnt_up+motif_cnt_dn)*sizeof(int)) ;
  }
  for (i=0 ; i<gene_count ; i++){
    for (j=0 ; j<motif_cnt_up+motif_cnt_dn ; j++){
      profile[i][j] = 0 ;
    }
  }
  int motif_count=0 ;
  for (i=0 ; i<motif_cnt_up ; i++){
    char fn[1000] ;
    sprintf(fn, "%s/%d.filtered.txt", motifdir_up, i) ;
    int m, n ;
    char **g;
    char **c ;
    int **data ;
    readIntTable(fn, &m, &n, &data, &g, &c) ;
    for (j=0 ; j<n ; j++){
      e.key  = strdup(g[j]) ;
      my_hsearch_r(e, FIND, &ep, h_ind);
      if (ep) {
	     int idx = (int)(ep->data) ;
	     profile[idx][motif_count] = data[j][0] ;
      }
    }
    motif_count++ ;
  }

  for (i=0 ; i<motif_cnt_dn ; i++){
    char fn[1000] ;
    sprintf(fn, "%s/%d.filtered.txt", motifdir_dn, i) ;
    int m, n ;
    char **g;
    char **c ;
    int **data ;
    readIntTable(fn, &m, &n, &data, &g, &c) ;
    for (j=0 ; j<n ; j++){
      e.key  = strdup(g[j]) ;
      my_hsearch_r(e, FIND, &ep, h_ind);
      if (ep) {
	int idx = (int)(ep->data) ;
	profile[idx][motif_count] = data[j][0] ;
      }
    }
    motif_count++ ;
  }
  
  

s_motif_list *sorted_motif_list ;
  sorted_motif_list = (s_motif_list*) malloc (motif_count*sizeof(s_motif_list)) ;
  for (i=0 ; i<motif_count ; i++){
    char ind[100] ;
    sprintf(ind, "%d", i) ;
    e.key = strdup(ind) ;
    my_hsearch_r(e, FIND, &ep, h_MI);
    float mi=0 ;
    if (ep) {
      mi = atof(ep->data) ;
    }
    free(e.key) ;
    sorted_motif_list[i].index = i ;
    sorted_motif_list[i].score = mi ;
    printf("%d\t%2.3f\n", i, mi) ;
  }
  printf("motifs sorted.\n") ; fflush(stdout) ;
  qsort((void*)sorted_motif_list, motif_count, sizeof(s_motif_list), CmpFunc) ;

  FILE *fsum = fopen (summaryfile, "w") ;
  fprintf(fsum, "index\tlocation\tmotif-seq\tmotif_structure\tmi-value\tseq mi-value\tfrequency\tz-score\trobustness\tp-value\n") ;
  // write a new summary file sorted by MI values
  for (i=0 ; i<motif_count ; i++){
    int id = sorted_motif_list[i].index ;
    if (id<motif_cnt_up){
      fprintf(fsum, "%s", summary1[id+1][0]) ;
      for (j=1 ; j<msum ; j++){
	fprintf(fsum, "\t%s", summary1[id+1][j]) ;
      }
      fprintf(fsum, "\n") ;
    }else{
      id = id-motif_cnt_up ;
      fprintf(fsum, "%s", summary2[id+1][0]) ;
      for (j=1 ; j<msum ; j++){
	fprintf(fsum, "\t%s", summary2[id+1][j]) ;
      }
      fprintf(fsum, "\n") ;
    }
  }
  fclose(fsum) ;

  int**modules ;
  int *module_memcount ;
  int* status ;
  int nmodules=-1 ;
  mimatrix          = (float**) malloc (motif_count*sizeof(float*)) ;
  modules           = (int**) malloc (motif_count*sizeof(int*)) ;
  status            = (int*) malloc (motif_count*sizeof(int)) ;
  module_memcount   = (int*) malloc (motif_count*sizeof(int)) ;
  for (i=0 ; i<motif_count ; i++){
    mimatrix[i] = (float*) malloc (motif_count*sizeof(float)) ;
    modules[i]  = (int*) malloc (motif_count*sizeof(int)) ;
    for (j=0 ; j<motif_count ; j++){
      modules[i][j]  = -1 ;
    }
    status[i] = -1 ;
    module_memcount[i] = 0 ;
  }
  for (i=0 ; i<motif_count ; i++){
    printf("\n") ;
    int curr = status[i] ;
    if (status[i] == -1){
      nmodules++ ;
      status[i] = nmodules ;
      curr = nmodules ;
      modules[curr][module_memcount[curr]] = i ;
      module_memcount[curr]++ ;
    }
    int id1 = sorted_motif_list[i].index ;

    M_q1 = i_matrix_column (profile, id1, gene_count) ;
    for (j=i+1 ; j<motif_count ; j++){
      printf(".") ;
      fflush(stdout) ;
      //if (status[j] != -1)
	     //continue ;

      int id2 = sorted_motif_list[j].index ;
      M_q2 = i_matrix_column (profile, id2, gene_count) ;

      //printf ("%d (%d)\t%d (%d)", i, id1, j, id2) ; getchar() ;
      
      int ov = 0 ;
      int s1 = 0 ;
      int s2 = 0 ;
      int k ;
      for (k=0 ; k<gene_count ; k++){
        if (M_q1[k]==1 && M_q2[k]==1){
          ov++ ;
        }
        if (M_q1[k]==1){
          s1++ ;
        }
        if (M_q2[k]==1){
          s2++ ;
        }
      }
      
      double po = log10(cumhyper    (ov, s1, s2, gene_count)) ;
      double pu = log10(cumhyper_u  (ov, s1, s2, gene_count)) ;

      float mi = CalculateMIbasic(M_q1, M_q2, gene_count, mbins, mbins);
      //float p  = evalSeed(M_q2, gene_count, mi, mbins, M_q2, mbins, shuffle) ;
      float r  = pearson_int (M_q1, M_q2, gene_count) ;
      //printf ("%d\t%d\t%f\n", id1, id2, mi) ; getchar() ;

      mimatrix[i][j] = mimatrix[j][i] = 0 ;
      

      if (po<max_p || pu<max_p){
        float p ;
        if (po<pu){
          p = -1 * po ;
        }else{
          p = pu ;
        }
        mimatrix[i][j] = mimatrix[j][i] = p ;
        
        if (status[j] != -1)
          continue ;
        status[j] = curr ;
        modules[curr][module_memcount[curr]] = j ;
        module_memcount[curr]++ ;
      }

      free(M_q2) ;
    }
    free(M_q1) ;
  }
  printf("number of modules = %d\n", nmodules+1) ;

  f      = fopen(matrixfile, "w") ;
  if (!f)
    die("Cannot open matrixfile\n");

  int *order ;
  int cnt=0 ;
  order = (int*) malloc (motif_count*sizeof(int)) ;
  for (i=0 ; i<=nmodules ; i++){
    for (j=0 ; j<module_memcount[i] ; j++){
      int mo = modules[i][j] ;
      //fprintf(f, "\t%d", sorted_motif_list[mo].index) ;
      fprintf(f, "\t%d", mo) ;
      order[cnt++] = mo ;
    }
  }
  fprintf(f,"\n") ;
  for (i=0 ; i<motif_count ; i++){
    int id1 = order[i] ;
    //fprintf(f,"%d", sorted_motif_list[id1].index) ;
    fprintf(f,"%d", id1) ;
    for (j=0 ; j<motif_count ; j++){
      int id2 = order[j] ;
      if (id1 == id2){
	       fprintf(f,"\tinf") ;
      }else{
	       fprintf(f,"\t%2.4f", mimatrix[id1][id2]) ;
      }
    }
    fprintf(f,"\n") ;
  }

  fclose(f) ;
  return (0) ;
}

int sign(float v) {
  return v > 0 ? 1 : (v < 0 ? -1 : 0);
}
