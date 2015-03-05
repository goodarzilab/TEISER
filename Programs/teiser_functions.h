typedef struct _s_seed {
  int    index;      // where is it found
  double score;
  int    hits;
} s_seed;

typedef struct _s_motif_list {
  int    index;      // where is it found
  double score;
} s_motif_list;

float* read_expfile (char *expfile, s_sequence **sequences, int t_seq_count, char ***seq_names, int *seq_count) ;
int* get_motif_profile (s_motif *motif, s_sequence **sequences, int nseqs, struct my_hsearch_data *h_seq_ind, int *hits, float dG_t) ;
int* get_motif_profile_seq_only (s_motif *motif, s_sequence **sequences, int nseqs, struct my_hsearch_data *h_seq_ind, int *hits, float dG_t) ;
int CmpFunc(const void* _a, const void* _b) ;
double minCondInfoNormalized(s_motif **motifs, int DA, int nA, int* B_q, int DB, int* E_q, int DE, int n, double minr, int* midx, s_sequence **sequences, int nseqs, struct my_hsearch_data *h_seq_ind, float dG_t) ;
double CondInfoRatio(s_motif *motifA, s_motif *motifB, int DA, int DB, int* E_q, int DE, int n, s_sequence **sequences, int nseqs, struct my_hsearch_data *h_seq_ind, float dG_t) ;
double evalSeed(int* M_q, int RNA_count, float score, int mbins, int* E_q, int ebins, int shuffle) ;
double teiser_max_rank_test(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle) ;
int teiser_jn_max_rank_test(int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle, float max_p, int jn, int jn_f) ;
void modify_base ( s_motif *source_motif, s_motif **modified_motifs, int position ) ;
void elongate_motif (s_motif *source_motif, s_motif **modified_motifs) ;
float teiser_z_score_test(double score, int* M_q, int mbins, int* E_q, int ebins, int nborfs, int shuffle) ;
