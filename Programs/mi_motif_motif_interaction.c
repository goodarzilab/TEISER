#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>
#include "sys/types.h"
#include "string.h"
#include <math.h>

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

void refine_motif_profile (int **M_q, int ind, int *E_q, int n, float **pv, float maxp) ;
int sign(float v) ;

int main(int argc, char ** argv) {
	int      i,j,k ;
	float    dG_t = -1 ;

	char     *motiffile ;
	char     *profiledir ;
	char     *expfile ;
	char     *rna_fastafile ;
	char     *matrixfile ;
	char     *pvfile ;
	
	char     **pvcolnames ;
	char     **pvrownames ;
	int      pvrow_num ;
	int      pvcol_num ;
	float    **pvmatrix ;
	
	float    **mimatrix ;
	
	float    *E ;
	int      *E_q ;
	int      *M_q1 ;
	int      *M_q2 ;
	
	struct   my_hsearch_data *h_rna_ind ;
	ENTRY    e ;
	ENTRY*   ep ;
	int      hashret =0 ;
	
	int      quantized = 1 ;
	int      ebins = 0 ;
	int      mbins = 2 ;
	int      divbins = 50 ;
	float*   E_q_bins = 0 ;
	
	int      shuffle = 10000 ;
	
	int      rnd_fasta = 0 ;
	float    max_p = 0.001 ;
	
	float    max_sig = 2.0 ;
	
	int      mincount = 5 ;
	
	s_sequence **sequences ;
	s_motif  **motifs ;
	s_motif  **opt_motifs ;
	
	char     **seq_names ;
	int      seq_count ;
	int      t_seq_count ;
	int      motif_count = 0 ;
	
	motiffile        = get_parameter(argc, argv, "-motiffile") ;
	profiledir       = get_parameter(argc, argv, "-profiledir") ;
	rna_fastafile    = get_parameter(argc, argv, "-rna_fastafile") ;
	matrixfile       = get_parameter(argc, argv, "-matrixfile") ;
	pvfile           = get_parameter(argc, argv, "-pvfile") ;
	
	expfile          = get_parameter(argc, argv, "-expfile") ;
	quantized        = atoi(get_parameter(argc, argv, "-quantized"));
	
	if (exist_parameter(argc, argv, "-max_p")) {
		max_p          = atof(get_parameter(argc, argv, "-max_p"));
	}
	max_p = log10(max_p) ;
	if (exist_parameter(argc, argv, "-ebins")) {
		ebins          = atoi(get_parameter(argc, argv, "-ebins"));
	}
	if (exist_parameter(argc, argv, "-shuffle")) {
		shuffle        = atoi(get_parameter(argc, argv, "-shuffle"));
	}
	if (exist_parameter(argc, argv, "-dG_t")) {
    	dG_t         = atof(get_parameter(argc, argv, "-dG_t"));
  	}
	

	FILE *f, *fmotif ;
	FILE *fptr = fopen ( motiffile, "rb") ;
	if (!fptr){
		printf("Could not open the seed file...\n") ;
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
	
	readFloatTable(pvfile, &pvrow_num, &pvcol_num, &pvmatrix, &pvrownames, &pvcolnames, 0, 1) ;
		
	int hits =0;
	
	s_motif_list *sorted_motif_list ;
	sorted_motif_list = (s_motif_list*) malloc (motif_count*sizeof(s_motif_list)) ;
	for (i=0 ; i<motif_count ; i++){
		M_q1 = get_motif_profile (motifs[i], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;

		float mi = CalculateMIbasic(M_q1, E_q, seq_count, mbins, ebins) ;

		sorted_motif_list[i].index = i ;
		sorted_motif_list[i].score = mi ;

		free(M_q1) ;
	}
	printf("motifs sorted.\n") ; fflush(stdout) ;
	qsort((void*)sorted_motif_list, motif_count, sizeof(s_motif_list), CmpFunc) ;

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
		M_q1 = get_motif_profile (motifs[id1], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
		refine_motif_profile (&M_q1, id1, E_q, seq_count, pvmatrix, max_sig) ;
		
		//save modified profile to file
		char fn[1000] ;
		sprintf(fn, "%s/%d.filtered.txt", profiledir, id1) ;
		FILE *ft = fopen (fn, "w") ;
		if (!ft){
			die("Couldn't write filtered profile\n") ;
		}
		fprintf(ft, "gene\tfiltered motif\n") ;
		for (j=0 ; j<seq_count ; j++){
			fprintf(ft, "%s\t%d\n", seq_names[j], M_q1[j]) ;
		}
		fclose(ft) ;
		
		for (j=i+1 ; j<motif_count ; j++){
			printf(".") ;
			fflush(stdout) ;
			
			int id2 = sorted_motif_list[j].index ;
			M_q2 = get_motif_profile (motifs[id2], sequences, t_seq_count, h_rna_ind, &hits, dG_t) ;
			refine_motif_profile (&M_q2, id2, E_q, seq_count, pvmatrix, max_sig) ;
			
			int ov = 0 ;
			int s1 = 0 ;
			int s2 = 0 ;
			for (k=0 ; k<seq_count ; k++){
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
			
			double po = log10(cumhyper    (ov, s1, s2, seq_count)) ;
			double pu = log10(cumhyper_u  (ov, s1, s2, seq_count)) ;
			
			//printf("%d-%d: (%d-%d-%d-%d) %f\t%f (%f)\n", id1, id2, ov, s1, s2, seq_count, po, pu, max_p) ;
			
			float mi = CalculateMIbasic(M_q1, M_q2, seq_count, mbins, mbins);
			//float p  = evalSeed(M_q2, seq_count, mi, mbins, M_q2, mbins, shuffle) ;
			float r  = pearson_int (M_q1, M_q2, seq_count) ;
			
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
			fprintf(f, "\t%d", sorted_motif_list[mo].index) ;
			order[cnt++] = mo ;
		}
	}
	fprintf(f,"\n") ;
	for (i=0 ; i<motif_count ; i++){
		int id1 = order[i] ;
		fprintf(f,"%d", sorted_motif_list[id1].index) ;
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
	
	free(E) ;
	free(E_q) ;
	fclose(f) ;
	return (0) ;
}

void refine_motif_profile (int **M_q, int ind, int *E_q, int n, float **pv, float maxp){
	int i ;
	for (i=0 ; i<n ; i++){
		int c = E_q[i] ;
		float p = pv[ind][c] ;
		if (p<maxp){
			if ((*M_q)[i] == 1)
				(*M_q)[i] = 0 ;
		}
	}
}

int sign(float v) {
	return v > 0 ? 1 : (v < 0 ? -1 : 0);
}
