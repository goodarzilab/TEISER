//$TEISERDIR/Programs/mi_report_motif_family -expfile /Users/hani/Life/Projects/Ongoing_Projects/Alternative_splicing/skipped_exons/se_events_MDA_cont.txt -motiffile /Users/hani/Life/Projects/Ongoing_Projects/Alternative_splicing/skipped_exons/se_events_MDA_cont.txt.select.bin -quantized 0 -rna_fastafile /Users/hani/Life/Projects/Ongoing_Projects/Alternative_splicing/skipped_exons/se_50bp_MDA.fa -reportfile /Users/hani/Life/Projects/Ongoing_Projects/Alternative_splicing/skipped_exons/se_events_MDA_cont.txt.report -matrixfile /Users/hani/Life/Projects/Ongoing_Projects/Alternative_splicing/skipped_exons/se_events_MDA_cont.txt.matrix -summaryfile /Users/hani/Life/Projects/Ongoing_Projects/Alternative_splicing/skipped_exons/se_events_MDA_cont.txt.summary -ebins 30  -divbins 50 
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

  char     *motiffile ;
  char     *rna_fastafile ;
  char     *expfile ;
  char     *reportfile ;
  char     *profile ;
  char     *matrixfile ;
  char     *summaryfile ;
  
  
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
  double   *pvalue ;

  motiffile        = get_parameter(argc, argv, "-motiffile") ;
  rna_fastafile    = get_parameter(argc, argv, "-rna_fastafile") ;
  reportfile       = get_parameter(argc, argv, "-reportfile") ;
  matrixfile       = get_parameter(argc, argv, "-matrixfile") ;
  summaryfile      = get_parameter(argc, argv, "-summaryfile") ;
  profile          = get_parameter(argc, argv, "-profile") ;
  
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
 pvalue = (double*)malloc(cluster_num*sizeof(double)) ;
 
 printf("Done\n");

 int hits ;
 freport = fopen ( reportfile, "w") ;
 if (!freport){
   die ("Couldn't open the report file...\n") ;
 }
 
 int *inst = (int*) calloc (10000, sizeof(int)) ;
 int ninst=0 ;
 M_q = (int*) calloc (seq_count, sizeof(int)) ;
 int nhits=0 ;
 for (i=0 ; i<motif_count ; i++){
   printf("Checking motif %d...", i) ; fflush(stdout) ;

   //get motif positions and write to file
   //fprintf(freport, ">motif %d\t%s\n", i, print_motif_to_char(motifs[i])) ;
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
     int *tmp = find_motif_instance_positions ( motifs[i], sequences[j], &n_inst, dG_t) ;

     if (n_inst==0){
      free(tmp) ;
      continue ;
     }
     //fprintf(freport, "%s\t%d", sequences[j]->name, n_inst) ;
     for (k=0 ; k<n_inst ; k++){
        fprintf(freport, ">%s_%s_%d_%d\n", sequences[j]->name, print_motif_to_char(motifs[i]), k, tmp[k]) ;
        int m ;
        for (m=tmp[k]-10 ; m<tmp[k]+motifs[i]->linear_length+10 ; m++){
          switch (sequences[j]->bases[m]){
            case _U:
              //printf("T");
              fprintf(freport, "T");
              break ;
            case _A:
              //printf("A");
              fprintf(freport, "A");
              break ;
            case _G:
              //printf("G");
              fprintf(freport, "G");
              break ;
            case _C:
              //printf("C");
              fprintf(freport, "C");
              break ;
            default:
              //printf("N");
              fprintf(freport, "N");
          }
        }
        fprintf(freport, "\n");
        //printf("\n") ;
        //printf("%s", print_motif_to_char(motifs[i])) ;
        //getchar() ;
        //fprintf(freport, "\t%d", tmp[k]) ;
     }
     fprintf(freport, "\n") ;
     fflush(freport) ;

     free(tmp) ;
   }
   printf("positions mapped.\n") ; fflush(stdout) ;

   int *M = get_motif_profile (motifs[i], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;

   for (j=0 ; j<seq_count ; j++){
      if (M[j]>0 && M_q[j]==0){
        M_q[j]=M[j] ;
        (nhits)++ ;
      }
   }
   free(M) ;
 }
 fclose(freport) ;

 k=0 ;
 /*if (ninst>0){
  fprintf(freport, "%s\t%d", sequences[j]->name, ninst) ;
  for (k=0 ; k<ninst ; k++){
    fprintf(freport, "\t%d", inst[k]) ;
  }
  fprintf(freport, "\n") ;
 }*/


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

FILE *fpro = fopen ( profile, "wt") ;
 if (!fpro){
   die ("Couldn't open the profile...\n") ;
 }
for (j=0 ; j<seq_count ; j++){
  fprintf(fpro, "%s\t%d\n", seq_names[j], M_q[j]) ;
}
fclose(fpro) ;

printf("doing stats: ") ;
float mi = CalculateMIbasic(M_q, E_q, seq_count, mbins, ebins) ;
printf("mi = %f\n", mi) ; fflush(stdout) ;
double pass = evalSeed(M_q, seq_count, mi, mbins, E_q, ebins, 10000) ;
printf("pass = %f\n", pass) ; fflush(stdout) ;
float z = teiser_z_score_test(mi, M_q, mbins, E_q, ebins, seq_count, 10000) ;
printf("z = %f\n", z) ; fflush(stdout) ;
int jn_test = teiser_jn_max_rank_test(M_q, mbins, E_q, ebins, seq_count, 1000, 0.0001, 10, 3) ;

char p_val[100] ;
if (pass == 0){
  sprintf(p_val, "<%.2e", 1.0/10000) ;
}else{
  sprintf(p_val, "<%.6e", pass) ;
}

f = fopen(summaryfile, "w") ;
fprintf(f, "index\tlocation\tmotif-seq\tmotif_structure\tmi-value\tseq mi-value\tfrequency\tz-score\trobustness\tp-value\n") ;
fprintf(f, "%d\tlocation\tfamily\tstructure\t%3.6f\t0\t%3.3f\t%3.3f\t%d/10\t%s\n", 0, mi, ((float)nhits)/seq_count, z, jn_test, p_val) ;


 free(M_q);

 for (j=0 ; j<cluster_num ; j++){
  printf("freq %d: %f\t%d\t%d\t%d\t%d\n",j,  ((float)cluster_motif_count[j])/cluster_gene_count[j], cluster_motif_count[j], cluster_gene_count[j], nhits, seq_count);
  double pu = log10(cumhyper   (cluster_motif_count[j], cluster_gene_count[j], nhits, seq_count));
  double pd = log10(cumhyper_u (cluster_motif_count[j], cluster_gene_count[j], nhits, seq_count));
  if(pu<pd){
    pvalue[j] = -1*pu ;
  }else{
    pvalue[j] = pd ;
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
 fprintf(f, "0") ;
 // print row of p-values
 for (i=0 ; i<cluster_num ; i++){
  fprintf(f, "\t%3.4f", pvalue[i]) ;
 }
 fprintf(f, "\n") ;


 free(E) ;
 free(E_q) ;
 fclose(f) ;
 fclose(freport) ;
 return (0) ;
}
