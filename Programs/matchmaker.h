void get_linear_length( s_motif *motif ) ;
int lcl_is_a( NUCBIT child, NUCBIT parent ) ;
int lcl_pair_bases( NUCBIT base1, NUCBIT base2 ) ;
int lcl_match_motif_seq( s_motif *motif, NUCBIT *sequence, int b, int l, float dG_t ) ;
int find_motif_instance( s_motif *motif, s_sequence *sequence, float dG_t ) ;
int* find_motif_instance_positions ( s_motif *motif, s_sequence *sequence, int *hits, float dG_t ) ;
int lcl_match_motif_seq_only( s_motif *motif, NUCBIT *sequence ) ;
int find_motif_instance_seq_only( s_motif *motif, s_sequence *sequence ) ;
NUCBIT lcl_complement_base( NUCBIT base ) ;

