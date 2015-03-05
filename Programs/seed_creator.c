#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>

#include "dataio.h"
#include "structures.h"
#include "nucleotides.h"
#include "create_motifs.h"
#include "read_write_motif.h"

int main(int argc, char ** argv) {
  // GENERAL
  int i ;
  
  int num_motifs_per_file = 2000000 ;

  int min_stem_length ; //Minimum stem length
  int max_stem_length ; //Maximum stem length
  int min_loop_length ; //Minimum loop length
  int max_loop_length ; //Maximum loop length
  int min_inf_seq ;     //Minimum number of informative bases
  int max_inf_seq ;     //Maximum number of informative bases
  float min_inf ;       //Minimum information
  float max_inf ;       //Maximum information

  char *outfile ;       //Output file

  if (argc == 1) {
    printf("Usage : seed_creator -min_stem_length INT -max_stem_length INT -min_loop_length INT -max_loop_length INT -min_inf_seq INT -max_inf_seq INT -max_inf FLOAT -min_inf FLOAT -outfile FILE\n");
    exit(0);
  }

  min_stem_length = atoi(get_parameter(argc, argv, "-min_stem_length")) ;
  max_stem_length = atoi(get_parameter(argc, argv, "-max_stem_length")) ;
  min_loop_length = atoi(get_parameter(argc, argv, "-min_loop_length")) ;
  max_loop_length = atoi(get_parameter(argc, argv, "-max_loop_length")) ;
  min_inf_seq     = atoi(get_parameter(argc, argv, "-min_inf_seq")) ;
  max_inf_seq     = atoi(get_parameter(argc, argv, "-max_inf_seq")) ;
  min_inf         = atof(get_parameter(argc, argv, "-min_inf")) ;
  max_inf         = atof(get_parameter(argc, argv, "-max_inf")) ;

  outfile         = get_parameter(argc, argv, "-outfile");

  s_motif **motifs;
  //int num_motifs = create_motifs( NULL, min_stem_length, max_stem_length, min_loop_length, max_loop_length, min_inf_seq, max_inf_seq ) ;
  //printf("%d seeds would be created.\n", num_motifs) ;

  motifs = (s_motif**) malloc (num_motifs_per_file * sizeof(s_motif*)) ;

  printf("Creating motifs... ") ;
  fflush(stdout) ;
  int num_motifs = create_and_write_motifs( outfile, num_motifs_per_file, motifs, min_stem_length, max_stem_length, min_loop_length, max_loop_length, min_inf_seq, max_inf_seq, min_inf, max_inf );
  printf("%d seeds were created.\n", num_motifs) ;

  free(motifs) ;

  return 0 ;
}
