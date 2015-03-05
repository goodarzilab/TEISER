#include "stdio.h"
#include "stdlib.h"
#include "sys/types.h"
#include "string.h"
#include "search.h"
#include "math.h"

#include "nucleotides.h"
#include "structures.h"
#include "sequences.h"
#include "matchmaker.h"
#include "dataio.h"
#include "hashtable.h"
#include "statistics.h"
#include "information.h"
#include "read_write_motif.h"
#include "teiser_functions.h"

float* read_expfile (char *expfile, s_sequence **sequences, int t_seq_count, char ***seq_names, int *seq_count){
  int i ;
  struct my_hsearch_data *h_exp ;
  ENTRY e;
  ENTRY *ep;
  int mynmax = 1000000 ;
  int   hashret =0 ;
  char* s ;
  char* buff ;
  char* genename ;
  FILE  *fp ;
  
  buff = (char*)malloc(mynmax * sizeof(char)) ;

  printf("Reading expfile...") ;
  fflush(stdout) ;
  char **rna_names = (char**) malloc ( t_seq_count * sizeof(char*)) ;
  int  rna_count = 0 ;

  fp = fopen(expfile, "r") ;
  if (!fp){
    printf("read_expfile: please enter a valid filename (%s invalid)\n", expfile) ;
    exit(0) ;
  }

  h_exp = (struct my_hsearch_data*)calloc(1,  sizeof(struct my_hsearch_data));
  hashret = my_hcreate_r(1000000, h_exp);
  if (hashret == 0) {
    printf("read_expfile: couldn't make the hashtable...\n");
    exit(0);
  }

  while (!feof(fp)){
    fgets(buff, mynmax, fp);
    if (feof(fp))
      break;

    chomp(buff);
    genename = mystrtok(buff, '\t');
    s = mystrtok(0, '\t') ;

    e.key  = strdup(genename) ;
    e.data = strdup(s) ;

    hashret = my_hsearch_r(e, ENTER, &ep, h_exp);
    if (hashret == 0){
      printf("read_expfile: couldn't add the data to hashtable...\n");
      exit(0);
    }
    free(s) ;
    free(genename) ;
  }

  float *E = (float*) malloc ( t_seq_count * sizeof(float)) ;

  for (i=0 ; i< t_seq_count ; i++){
    e.key  = strdup(sequences[i]->name) ;
    my_hsearch_r(e, FIND, &ep, h_exp);
    if (!ep){
      continue ;
    }
    s = strdup(ep->data) ;
    E[rna_count] = atof(s) ;
    rna_names[rna_count++] = strdup(sequences[i]->name) ;
    free(s) ;
  }

  (*seq_names) = (char**) malloc (rna_count * sizeof(char*)) ;
  for (i=0 ; i<rna_count ; i++){
    (*seq_names)[i] = rna_names[i] ;
  }
  *seq_count = rna_count ;

  my_hdestroy_r(h_exp) ;
  printf("Done\n") ;
  fflush(stdout) ;

  return E ;
}

int* get_motif_profile (s_motif *motif, s_sequence **sequences, int nseqs, struct my_hsearch_data *h_seq_ind, int *hits, float dG_t){
  int i ;
  ENTRY e;
  ENTRY *ep;

  *hits = 0 ;
  int *M_q = (int*) malloc (nseqs * sizeof(int)) ;
  for ( i=0 ; i<nseqs ; i++){
    M_q[i] = 0 ;
  }
  for (i=0 ; i<nseqs ; i++){
    e.key = strdup (sequences[i]->name) ;
    my_hsearch_r(e, FIND, &ep, h_seq_ind);
    free (e.key) ;
    if (!ep){
      continue ;
    }
    int index = (int) ep->data ;
    int match = find_motif_instance ( motif, sequences[i], dG_t ) ;

    if (match != -1){
      M_q[index] = 1 ;
      (*hits)++ ;
    }
  }
  return M_q ;
}

int* get_motif_profile_seq_only (s_motif *motif, s_sequence **sequences, int nseqs, struct my_hsearch_data *h_seq_ind, int *hits, float dG_t){
  int i ;
  ENTRY e;
  ENTRY *ep;
  
  s_motif *lmotif = (s_motif*) malloc (sizeof(s_motif)) ;
  lmotif->num_phrases = motif->linear_length ;
  lmotif->linear_length = motif->linear_length ;
  lmotif->phrases = (s_phrase*) malloc (lmotif->num_phrases*sizeof(s_phrase)) ;
  for (i=0 ; i< lmotif->num_phrases ; i++){
    lmotif->phrases[i].structure = _leftBulge ;
  }
  int be=0, en=lmotif->linear_length-1 ;
  for (i=0 ; i<motif->num_phrases ; i++){
    if (motif->phrases[i].structure == _leftBulge){
      lmotif->phrases[be++].base = motif->phrases[i].base ;
    }else if (motif->phrases[i].structure == _rightBulge){
      lmotif->phrases[en--].base = motif->phrases[i].base ;
    }else if (motif->phrases[i].structure == _pair){
      lmotif->phrases[be++].base = motif->phrases[i].base ;
      lmotif->phrases[en--].base = lcl_complement_base(motif->phrases[i].base) ;
    }
  }
  int cnt ;
  int *ori = get_motif_profile (motif, sequences, nseqs, h_seq_ind, &cnt, -1) ;
  int *new = get_motif_profile (lmotif, sequences, nseqs, h_seq_ind, &cnt, -1) ;
  
  *hits = 0 ;
  int *M_q = (int*) malloc (nseqs * sizeof(int)) ;
  for ( i=0 ; i<nseqs ; i++){
    M_q[i] = 0 ;
  }
  for (i=0 ; i<nseqs ; i++){
    e.key = strdup (sequences[i]->name) ;
    my_hsearch_r(e, FIND, &ep, h_seq_ind);
    free (e.key) ;
    if (!ep){
      continue ;
    }
    int index = (int) ep->data ;
    if (new[index]==1 && ori[index]==0){
      M_q[index] = 1 ;
      (*hits)++ ;
    }
  }
  return M_q ;
}

double evalSeed(int* M_q, int RNA_count, float score, int mbins, int* E_q, int ebins, int shuffle) {
  double    pass ;

  pass = teiser_max_rank_test(score, M_q, mbins, E_q, ebins, RNA_count, shuffle) ;

  if (pass > 1){
    free(M_q);
    return 1;
  }

  return pass ;
}

int CmpFunc(const void* _a, const void* _b) {
  const s_seed* a = (const s_seed*) _a;
  const s_seed* b = (const s_seed*) _b;

  if (a->score < b->score)
    return 1;
  else if(a->score == b->score)
    return  0;
  else
    return -1;
}

double minCondInfoNormalized(s_motif **motifs, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, double minr, int* midx, s_sequence **sequences, int nseqs, struct my_hsearch_data *h_seq_ind, float dG_t) {
  int    i ;
  int*   AB_q;
  int    DAE, DAB;
  double mi_ab_e, mi_a_e, mi_a_b;
  double cmi;
  double minratio = 1e6;
  double ratio;
  int*   M_q;

  DAE    = DA * DE;
  DAB    = DA * DB;
  *midx  = -1;

  for (i=0; i<nA; i++) {  // nA is the number of already kept GO cats
    int hits ;
    M_q         = get_motif_profile (motifs[i], sequences, nseqs, h_seq_ind, &hits, dG_t) ;

    // calculate conditional information I(Cnew;E|Cexist)  = I(A;E|B)  =  I(A,B;E) - I(A;E)
    AB_q        = combineQuantizedVectors(M_q, B_q, n, DA, DB);
    mi_ab_e     = CalculateMIbasic       (AB_q,   E_q, n, DAB, DE);
    mi_a_e      = CalculateMIbasic       (M_q, E_q, n, DA, DE);
    cmi         = mi_ab_e - mi_a_e;

    // calculate MI I(Cnew;Cexist)
    mi_a_b      = CalculateMIbasic       (M_q, B_q, n, DA, DB);

    if (cmi < 1e-16)
      cmi = 0.0;
    else if (mi_a_b < 1e-16)
      mi_a_b = 1e-16;

    ratio       = cmi / mi_a_b;

    if (i == 0)
      minratio = ratio;
    else {

      if (ratio < minratio) {
	minratio = ratio;
      }
    }

    free(M_q);
    free(AB_q);

    /*******  break early feature *******/
    if (minratio < minr) {
      *midx = i;  // index of category that killed the current one
      return minratio;
    }
  }

  return minratio;
}

double CondInfoRatio(s_motif *motifA, s_motif *motifB, int DA, int DB, int* E_q, int DE, int n, s_sequence **sequences, int nseqs, struct my_hsearch_data *h_seq_ind, float dG_t){
  int    i ;
  int*   AB_q;
  int    DAE, DAB;
  double mi_ab_e, mi_a_e, mi_a_b;
  double cmi;
  double ratio;
  int*   M_qA;
  int*   M_qB;

  DAE    = DA * DE;
  DAB    = DA * DB;
  
  int hitsA ;
  M_qA        = get_motif_profile (motifA, sequences, nseqs, h_seq_ind, &hitsA, dG_t) ;
  int hitsB ;
  M_qB        = get_motif_profile (motifB, sequences, nseqs, h_seq_ind, &hitsB, dG_t) ;

  AB_q        = combineQuantizedVectors(M_qA, M_qB, n, DA, DB);
  mi_ab_e     = CalculateMIbasic       (AB_q,   E_q, n, DAB, DE);
  mi_a_e      = CalculateMIbasic       (M_qA, E_q, n, DA, DE);
  cmi         = mi_ab_e - mi_a_e;

  // calculate MI I(Cnew;Cexist)
  mi_a_b      = CalculateMIbasic       (M_qA, M_qB, n, DA, DB);

  if (cmi < 1e-16)
    cmi = 0.0;
  else if (mi_a_b < 1e-16)
    mi_a_b = 1e-16;

  ratio       = cmi / mi_a_b;

  free(M_qA);
  free(M_qB);
  free(AB_q);

  return ratio;

}
double teiser_max_rank_test(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle) {
  double val = 0 ;
  int      j;
  int*     E_q_shu;
  double    mymi;
  double*   a_mi_shu;
  double    c;

  a_mi_shu = (double*)malloc(shuffle * sizeof(double));

  for (j=0; j<shuffle; j++) {
    E_q_shu = shuffleInt(E_q, nborfs);
    mymi = CalculateMIbasic(M_q, E_q_shu, nborfs, mbins, ebins);
    a_mi_shu[j] = mymi;
    free(E_q_shu);
  }

  qsort((void*)a_mi_shu, shuffle, sizeof(double), CmpDblRegular);

  c  = a_mi_shu[ shuffle - 1 ];

  if ((c < 1e-10) && (score < 1e-10))  // shortcut
    val = shuffle;
  else if (score <= a_mi_shu[0]) {  // other shortcut
    val = shuffle;
  } else {
    j = shuffle - 1;
    while ((j >= 0) && (score <= a_mi_shu[j])) {   // go while the shuffled score is higher than the real one
      j --;
    }
    val = shuffle - j - 1;

  }
  free(a_mi_shu);

  val = val/shuffle ;
  return val;
}

float teiser_z_score_test(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle) {
  double val = 0 ;
  int      j,k;
  int*     E_q_shu;
  float    mymi;
  double    c;
  float sum=0, sum_2=0 ;

  for (j=0; j<shuffle; j++) {
    E_q_shu = shuffleInt(E_q, nborfs);
    mymi = CalculateMIbasic(M_q, E_q_shu, nborfs, mbins, ebins);
    sum+=mymi ;
    sum_2 += mymi*mymi ;
    free(E_q_shu);
  }
  float ave = sum/shuffle ;
  float std = sqrt((sum_2-sum*sum/shuffle) / (shuffle-1));

  float z = (score-ave)/std ;
  return z;
}

int teiser_jn_max_rank_test(int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, float max_p, int jn, int jn_f) {
  int*  M_q_cross;
  int*  E_q_cross;
  int   pass = 0;
  int   l;
  int   nborfs_retained;
  float r;
  float p ;
  int   j;
  double mymi;

  M_q_cross = (int*)malloc(nborfs * sizeof(int));
  E_q_cross = (int*)malloc(nborfs * sizeof(int));

  for (l=0; l<jn; l++) {
    // keep random 90% of the data
    nborfs_retained = 0;
    for (j=0; j<nborfs; j++) {
      r       = default_rand();

      if (r > 1.0/jn_f) {
        M_q_cross[ nborfs_retained ] = M_q[j];
        E_q_cross[ nborfs_retained ] = E_q[j];
        nborfs_retained ++;
      }
    }

    // eval mi
    mymi      = CalculateMIbasic(M_q_cross, E_q_cross, nborfs_retained, mbins, ebins);
    p         = teiser_max_rank_test(mymi, M_q_cross, mbins, E_q_cross, ebins, nborfs_retained, shuffle);
    if (p<max_p)
      pass++ ;
  }

  free(M_q_cross);
  free(E_q_cross);

  return pass ;
}

void modify_base ( s_motif *source_motif, s_motif **modified_motifs, int position ) {
  // creates an array of *s_motif having all possible combinations of bases at the specified position
  // modified_motifs will have 15 members
  // source_motif itself will also be one of the members of modified_motifs
  // position must be >= 0 and < source_motif ->num_phrases
  /*
    example of calling modify_base:
    s_motif *modified_motifs[ 15 ];
    modify_base( seeds[ this_seed ], modified_motifs, 0 );
  */
  memset( modified_motifs, 0, sizeof(s_motif*) * 15 );
  if( position < 0 || position >= source_motif ->num_phrases )
    return;

  int num_motifs = 0;
  NUCBIT base = 0 ;
  for ( base = 1; base <= _N; base ++ ) {
    // copy the contents of the source motif
    modified_motifs [num_motifs] = copy_motif (source_motif) ;

    // change the specified position to the new base composition
    modified_motifs [num_motifs]->phrases[position].base = base;
    num_motifs ++;
  }
}

void elongate_motif (s_motif *source_motif, s_motif **modified_motifs) {
  // creates an array of *s_motif each having one additional phrase compared to source_motif
  // modified_motifs will have 46 members
  // source_motif itself will also be one of the members of modified_motifs
  /*
    example of calling modify_base:
    
    s_motif *modified_motifs[ 46 ];
    modify_base( seeds[ this_seed ], modified_motifs );
  */
  memset( modified_motifs, 0, sizeof(s_motif*) * 46 );
  
  int num_motifs = 0;
  // the first motif is the same as the source_motif
  modified_motifs [num_motifs] = copy_motif (source_motif) ;

  num_motifs++;
  BYTE structure = 0 ;
  NUCBIT base = 1 ;
  for( structure = 1; structure <= 3; structure ++ )
    for( base = 1; base <= _N; base ++ ) {
      modified_motifs [num_motifs] = (s_motif*) malloc (sizeof(s_motif)) ;
      
      // the new motif has one more phrase compared to source_motif
      modified_motifs [num_motifs]->num_phrases = source_motif->num_phrases + 1;
      // determine the linear length of the new motif
      modified_motifs [num_motifs]->linear_length = source_motif->linear_length + ( (structure==_pair)? 2:1 );

      // copy the structure and sequence of the source_motif
      modified_motifs [num_motifs]->phrases = (s_phrase*) malloc (modified_motifs[num_motifs]->num_phrases*sizeof(s_phrase)) ;
      memcpy (modified_motifs[num_motifs]->phrases+1, source_motif->phrases, sizeof(s_phrase)*source_motif->num_phrases );
      
      // set the new phrase
      modified_motifs [num_motifs]->phrases[0].structure = structure;
      modified_motifs [num_motifs]->phrases[0].base = base;

      num_motifs ++;
    }
}
