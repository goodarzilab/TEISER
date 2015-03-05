double lcl_get_pairing_energy( NUCBIT base1, NUCBIT base2 ) ;
double lcl_get_free_folding_energy( s_motif *motif, NUCBIT *sequence ) ;
void lcl_get_surrounding_sequence(NUCBIT *source_seq, int seq_length, int exclude_start_pos, int exclude_length, NUCBIT *surrounding_seq, int *surrounding_length ) ;
double lcl_get_minimum_energy( NUCBIT *sequence, int seq_length ) ;
double get_dG_ratio_to_optimum( s_motif *motif, NUCBIT *sequence, int seq_length, int motif_start_pos ) ;