CFLAGS = -O3 -g
CC = gcc

all : seed_creator mi_find_seed mi_optimize mi_report motif_page_run mi_motif_motif_interaction combine_motifs final_motif_motif_interaction mi_conserve

mi_cond_inf_ratio : mi_cond_inf_ratio.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_cond_inf_ratio mi_cond_inf_ratio.c  teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

scan_sequences_for_motif_match_count : scan_sequences_for_motif_match_count.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o scan_sequences_for_motif_match_count scan_sequences_for_motif_match_count.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o -lm $(LFLAGS)

make_variant_family_for_structural_motif : make_variant_family_for_structural_motif.c read_write_motif.o create_motifs.o dataio.o
	$(CC) $(CFLAGS) -Wall -o make_variant_family_for_structural_motif make_variant_family_for_structural_motif.c read_write_motif.o create_motifs.o dataio.o -lm $(LFLAGS)

mi_report_gene_list : mi_report_gene_list.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_report_gene_list mi_report_gene_list.c  teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

mi_report_motif_family : mi_report_motif_family.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_report_motif_family mi_report_motif_family.c  teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

mi_report_motif_family_count_based : mi_report_motif_family_count_based.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_report_motif_family_count_based mi_report_motif_family_count_based.c  teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

print_motifs : print_motifs.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o print_motifs print_motifs.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

combine_motifs_from_files : combine_motifs_from_files.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o combine_motifs_from_files combine_motifs_from_files.c  teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

scan_sequences_for_motif_profile : scan_sequences_for_motif_profile.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o scan_sequences_for_motif_profile scan_sequences_for_motif_profile.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o -lm $(LFLAGS)

scan_sequence_pairs_for_structural_motif_presence : scan_sequence_pairs_for_structural_motif_presence.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o scan_sequence_pairs_for_structural_motif_presence scan_sequence_pairs_for_structural_motif_presence.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o -lm $(LFLAGS)

find_threshold_for_mi_values : find_threshold_for_mi_values.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o find_threshold_for_mi_values find_threshold_for_mi_values.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o -lm $(LFLAGS)

partition_seeds : partition_seeds.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o partition_seeds partition_seeds.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o -lm $(LFLAGS)

mi_mi_seed : mi_mi_seed.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_mi_seed mi_mi_seed.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o -lm $(LFLAGS)

scan_sequence_pairs_for_seed_presence : scan_sequence_pairs_for_seed_presence.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o scan_sequence_pairs_for_seed_presence scan_sequence_pairs_for_seed_presence.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o -lm $(LFLAGS)

get_variant_family_for_structural_motif : get_variant_family_for_structural_motif.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o get_variant_family_for_structural_motif get_variant_family_for_structural_motif.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

manual_mutations_in_a_motif : manual_mutations_in_a_motif.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o manual_mutations_in_a_motif manual_mutations_in_a_motif.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

match_motif_to_seq : match_motif_to_seq.c match_motif_to_seq.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o match_motif_to_seq match_motif_to_seq.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

split_seeds_into_files : split_seeds_into_files.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o split_seeds_into_files split_seeds_into_files.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

combine_motifs_into_file : combine_motifs_into_file.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o combine_motifs_into_file combine_motifs_into_file.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

get_motif_at_index : get_motif_at_index.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o get_motif_at_index get_motif_at_index.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

mi_conserve : mi_conserve.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_conserve mi_conserve.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o -lm $(LFLAGS)

final_motif_motif_interaction : final_motif_motif_interaction.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o final_motif_motif_interaction final_motif_motif_interaction.c  teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

combine_motifs : combine_motifs.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o combine_motifs combine_motifs.c  teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

mi_motif_motif_interaction : mi_motif_motif_interaction.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_motif_motif_interaction mi_motif_motif_interaction.c  teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

motif_page_run: motif_page_run.c dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o  motif_page_run motif_page_run.c dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

mi_report : mi_report.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_report mi_report.c  teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

mi_optimize : mi_optimize.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_optimize mi_optimize.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o folding_energy.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o -lm $(LFLAGS)

mi_find_seed : mi_find_seed.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o
	$(CC) $(CFLAGS) -Wall -o mi_find_seed mi_find_seed.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o folding_energy.o read_write_motif.o readFASTA.o -lm $(LFLAGS)

seed_creator : seed_creator.c read_write_motif.o create_motifs.o dataio.o
	$(CC) $(CFLAGS) -Wall -o seed_creator seed_creator.c read_write_motif.o create_motifs.o dataio.o -lm $(LFLAGS)

teiser_functions.o : teiser_functions.c teiser_functions.h dataio.o sequences.o hashtable.o matchmaker.o read_write_motif.o information.o statistics.o
	$(CC) $(CFLAGS) -Wall -c teiser_functions.c

readFASTA.o : readFASTA.c readFASTA.h dataio.o sequences.o statistics.o
	$(CC) $(CFLAGS) -Wall -c readFASTA.c

matchmaker.o : matchmaker.c matchmaker.h folding_energy.o
	$(CC) $(CFLAGS) -Wall -c matchmaker.c

folding_energy.o : folding_energy.c folding_energy.h
	$(CC) $(CFLAGS) -Wall -c folding_energy.c

read_write_motif.o : read_write_motif.c read_write_motif.h
	$(CC) $(CFLAGS) -Wall -c read_write_motif.c

create_motifs.o : create_motifs.c create_motifs.h
	$(CC) $(CFLAGS) -Wall -c create_motifs.c

mi_library.o : mi_library.c mi_library.h
	$(CC) $(CFLAGS) -Wall -c mi_library.c

statistics.o : statistics.c statistics.h
	$(CC) $(CFLAGS) -Wall -c statistics.c

sequences.o : sequences.c sequences.h
	$(CC) $(CFLAGS) -Wall -c sequences.c

dataio.o : dataio.c dataio.h
	$(CC) $(CFLAGS) -Wall -c dataio.c

information.o : information.c information.h
	$(CC) $(CFLAGS) -Wall -c information.c

hashtable.o : hashtable.c hashtable.h
	$(CC) $(CFLAGS) -Wall -c hashtable.c

clean: 
	rm *.o
