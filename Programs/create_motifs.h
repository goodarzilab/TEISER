void lcl_create_structure ( s_motif *motif, int stem_length, int loop_length ) ;
int lcl_next_sequence ( NUCBIT *sequence, int length ) ;
int lcl_count_informative_bases ( NUCBIT *sequence, int length ) ;
void lcl_initialize_sequence ( NUCBIT *sequence, int length ) ;
int lcl_copy_sequence ( s_motif *motif, NUCBIT *sequence, int length ) ;
int lcl_create_motifs( char *output_file, int *file_num, int num_motifs_per_file, s_motif **motifs, int num_motifs, int *num_total_motifs, int stem_length, int loop_length, int min_inf_bases, int max_inf_bases, double minI, double maxI, FILE *fptr_distribution ) ;
int create_and_write_motifs( char *output_file, int num_motifs_per_file, s_motif **motifs, int min_stem_length, int max_stem_length, int min_loop_length, int max_loop_length, int min_inf_bases, int max_inf_bases, double minI, double maxI ) ;
double lcl_calculate_probability( s_motif *motif ) ;
double calcluate_information (s_motif *motif) ;
