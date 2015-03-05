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
  char     *rna_fastafile ;
  char     *expfile ;
  char     *reportfile ;
  char     *matrixfile ;
  char     *profiledir ;

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
  
  int      mincount = 5 ;

  s_sequence **sequences ;
  s_motif  **motifs ;

  char     **seq_names ;
  int      seq_count ;
  int      t_seq_count ;
  int      motif_count = 0 ;

  int      cluster_num =0 ;
  int      *cluster_gene_count ;
  double   **pvalue, **lpvalue ;

  motiffile        = get_parameter(argc, argv, "-motiffile") ;
  rna_fastafile    = get_parameter(argc, argv, "-rna_fastafile") ;
  reportfile       = get_parameter(argc, argv, "-reportfile") ;
  matrixfile       = get_parameter(argc, argv, "-matrixfile") ;
  profiledir       = get_parameter(argc, argv, "-profiledir") ;

  expfile          = get_parameter(argc, argv, "-expfile") ;
  quantized        = atoi(get_parameter(argc, argv, "-quantized"));
  
  if (exist_parameter(argc, argv, "-ebins")) {
    ebins          = atoi(get_parameter(argc, argv, "-ebins"));
  }
  if (exist_parameter(argc, argv, "-divbins")) {
	divbins        = atoi(get_parameter(argc, argv, "-divbins"));
  }
  if (exist_parameter(argc, argv, "-dG_t")) {
    dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  }

  FILE *f, *freport ;
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

  // defining the number of clusters
  for (i=0 ; i<seq_count ; i++){
    if (E_q[i]>cluster_num)
      cluster_num = E_q[i] ;
  }
 cluster_num++ ;
 printf("Number of clusters: %d\n", cluster_num) ;

 // counting the number of genes in each cluster
 cluster_gene_count = (int*)malloc(cluster_num * sizeof(int)) ;
 for (i=0 ; i<cluster_num ; i++){
   cluster_gene_count[i] = 0 ;
 }
 for (i=0 ; i<seq_count ; i++){
   cluster_gene_count[E_q[i]]++ ;
 }

 //calculating p-values
 printf("calculating the p-value matrix.\n") ;
 printf("Allocating memory ... ");
 pvalue = (double**)malloc(cluster_num*sizeof(double*)) ;
 lpvalue = (double**)malloc(cluster_num*sizeof(double*)) ;
 for (i=0 ; i<cluster_num ; i++){
   pvalue [i] = (double*)malloc(motif_count*sizeof(double)) ;
   lpvalue [i] = (double*)malloc(motif_count*sizeof(double)) ;
 }
 printf("Done\n");

 int hits, lhits ;
 freport = fopen ( reportfile, "w") ;
 if (!freport){
   die ("Couldn't open the report file...\n") ;
 }
 
 for (i=0 ; i<motif_count ; i++){
   //get motif positions and write to file
   fprintf(freport, ">motif %d\t%s\n", i, print_motif_to_char(motifs[i])) ;
   ENTRY e;
   ENTRY *ep;
   for (j=0 ; j<t_seq_count ; j++){
     e.key = strdup (sequences[j]->name) ;
     my_hsearch_r(e, FIND, &ep, h_rna_ind);
     free (e.key) ;
     if (!ep){
       continue ;
     }
     int n_inst=0 ;
     int *inst = find_motif_instance_positions ( motifs[i], sequences[j], &n_inst, dG_t) ;
     int k=0 ;
     if (n_inst>0){
       fprintf(freport, "%s\t%d", sequences[j]->name, n_inst) ;
       for (k=0 ; k<n_inst ; k++){
	 fprintf(freport, "\t%d", inst[k]) ;
       }
       fprintf(freport, "\n") ;
     }
     free(inst) ;
   }

   M_q = get_motif_profile (motifs[i], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
   int *M_qs = get_motif_profile_seq_only (motifs[i], sequences, t_seq_count, h_rna_ind, &lhits, dG_t) ;
   FILE* fp ;
   char profile[1000] ;
   sprintf(profile, "%s/%d.txt", profiledir, i) ;
   fp = fopen(profile, "w") ;
   if (!fp)
     die("Couldn't write the profile...") ;
   fprintf(fp, "Gene\tMotif\n") ;
   for (j=0 ; j<seq_count ; j++){
     fprintf(fp, "%s\t%d\n", seq_names[j], M_q[j]) ;
   }
   fclose(fp);

   int *cluster_motif_count ; //count the targets in each cluster

   cluster_motif_count = (int*)malloc(cluster_num * sizeof(int)) ;
   for (j=0 ; j<cluster_num ; j++){
     cluster_motif_count[j]=0 ;
   }
   for (j=0 ; j<seq_count ; j++){
     if (M_q[j]==1){
       cluster_motif_count[E_q[j]]++ ;
     }
   }

   int *lcluster_motif_count ; //count the targets in each cluster

   lcluster_motif_count = (int*)malloc(cluster_num * sizeof(int)) ;
   for (j=0 ; j<cluster_num ; j++){
     lcluster_motif_count[j]=0 ;
   }
   for (j=0 ; j<seq_count ; j++){
     if (M_qs[j]==1){
       lcluster_motif_count[E_q[j]]++ ;
     }
   }

   free(M_q);
   free(M_qs);

   for (j=0 ; j<cluster_num ; j++){
    double pu = log10(cumhyper   (cluster_motif_count[j], cluster_gene_count[j], hits, seq_count));
    double pd = log10(cumhyper_u (cluster_motif_count[j], cluster_gene_count[j], hits, seq_count));
    if(pu<pd){
      pvalue[j][i] = -1*pu ;
    }else{
      pvalue[j][i] = pd ;
    }
    pu = log10(cumhyper   (lcluster_motif_count[j], cluster_gene_count[j], lhits, seq_count));
    pd = log10(cumhyper_u (lcluster_motif_count[j], cluster_gene_count[j], lhits, seq_count));
    if(pu<pd){
      lpvalue[j][i] = -1*pu ;
    }else{
      lpvalue[j][i] = pd ;
    }
   }
 }

 f = fopen(matrixfile, "w") ;
 if (!f)
   die("Cannot open datafile\n");
 
 fprintf(f, "MOTIF") ;
 // row of cluster indices
 for (i=0 ; i<cluster_num ; i++){
   if (quantized == 1) {
     fprintf(f, "\t%d", i) ;
   }else{
     fprintf(f, "\t[%3.2f %3.2f]", E_q_bins[i], E_q_bins[i+1]);
   }
 }
 fprintf(f, "\n");

 // go thru all retained categories
 for(j=0 ; j<motif_count ; j++) {
   // print desc
   fprintf(f, "%d", j) ;

   // print row of p-values
   for (i=0 ; i<cluster_num ; i++){
     fprintf(f, "\t%3.4f", pvalue[i][j]) ;
   }
   fprintf(f, "\n") ;
 }
 fclose(f) ;

 char buffer[500] ;
 sprintf (buffer, "%s.linear", matrixfile) ;
 f = fopen(buffer, "w") ;
 if (!f)
   die("Cannot open datafile\n");
 
 fprintf(f, "MOTIF") ;
 // row of cluster indices
 for (i=0 ; i<cluster_num ; i++){
   if (quantized == 1) {
     fprintf(f, "\t%d", i) ;
   }else{
     fprintf(f, "\t[%3.2f %3.2f]", E_q_bins[i], E_q_bins[i+1]);
   }
 }
 fprintf(f, "\n");

 // go thru all retained categories
 for(j=0 ; j<motif_count ; j++) {
   // print desc
   fprintf(f, "%d", j) ;

   // print row of p-values
   for (i=0 ; i<cluster_num ; i++){
     fprintf(f, "\t%3.4f", lpvalue[i][j]) ;
   }
   fprintf(f, "\n") ;
 }
 fclose(f) ;

 free(E) ;
 free(E_q) ;
 fclose(f) ;
 fclose(freport) ;
 return (0) ;
}
