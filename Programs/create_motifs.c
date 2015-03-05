#include <stdio.h>
#include <stdlib.h>
#include "structures.h"
#include "nucleotides.h"
#include "create_motifs.h"
#include "teiser_functions.h"

void lcl_create_structure( s_motif *motif, int stem_length, int loop_length ) {
// creates a stem-loop structure with the specified stem length and loop length
// all bases will be N
  int i=0 ;
  motif->num_phrases = stem_length + loop_length;
  motif->linear_length = stem_length * 2 + loop_length;

  motif->phrases = (s_phrase*) malloc (motif->num_phrases*sizeof(s_phrase)) ;

  for ( i=0; i<stem_length; i++ ) {
    motif->phrases[i].structure = _pair;
    motif->phrases[i].base = _N;
  }

  for ( i=0; i<loop_length; i++ ){
    motif->phrases[stem_length+i].structure = _leftBulge;
    motif->phrases[stem_length+i].base = _N;
  }
}


int lcl_next_sequence( NUCBIT *sequence, int length ) {
// returns true if there is a next sequence, otherwise false
/* A sequence is treated as a base-4 number, where the 5' nucleotide is the right digit and the 3' nucleotide is the left digit.
The function next_sequence generates sequences by increasing this number by 1. */
  int i=0 ;
  for ( i=0 ; i<length ; i++ ) {
    if ( sequence[i] < LAST_LETTER ) {
      sequence[i] <<= 1;
      return 1;
    }else if( sequence[i] != _N ){
      sequence[i] = _N;
      return 1;
    }else
      sequence[i] = 1;
  }
  
  return 0;
}

int lcl_count_informative_bases( NUCBIT *sequence, int length ) {
// counts the number of bases that are not N
  int count=0; int i=0 ;
  for ( i=0; i<length; i++ )
    if ( sequence[i] != _N )
      count++;

  return count;
}

void lcl_initialize_sequence( NUCBIT *sequence, int length ) {
// sets all bases of the sequence to U
  int i=0 ;
  for( i=0 ; i<length; i++ )
    sequence[i] = 1;
}

int lcl_copy_sequence( s_motif *motif, NUCBIT *sequence, int length ) {
// copies the specified sequence to the motif
// returns true if successful, false otherwise

  if( motif->num_phrases != length ){
    printf("ERROR: Invalid copy request.") ;
    return 0;
  }
  int i=0 ;
  for( i=0; i<length; i++ )
    motif->phrases[i].base = sequence[i];

  return 1;
}

int lcl_create_motifs( char *output_file, int *file_num, int num_motifs_per_file, s_motif **motifs, int num_motifs, int *num_total_motifs, int stem_length, int loop_length, int min_inf_bases, int max_inf_bases, double minI, double maxI, FILE *fptr_distribution ) {
  // every time that num_motifs_per_file motifs are created, they will be written in the output and the array is flushed
  // returns the number of motifs in the array that are not written yet
  // if motifs==NULL, it only counts the number of motifs and does not store or write them
  
  int num_phrases = stem_length + loop_length;
  NUCBIT *sequence = (NUCBIT*) malloc (num_phrases*sizeof(NUCBIT)) ;
  lcl_initialize_sequence( sequence, num_phrases );
  
  while (-1) {
    int num_inf_bases = lcl_count_informative_bases( sequence, num_phrases );
    if( num_inf_bases >= min_inf_bases && num_inf_bases <= max_inf_bases ) {
      // the number of informative bases in this sequence is in the specified range
      motifs [num_motifs] = (s_motif*) malloc (sizeof(s_motif));
      lcl_create_structure( motifs[ num_motifs ], stem_length, loop_length );
      lcl_copy_sequence( motifs[ num_motifs ], sequence, num_phrases );
      
      // calculate the total information of the new motif and compare it with the specified range
      double I = calcluate_information( motifs[ num_motifs ] );
      if( fptr_distribution ){
	fprintf( fptr_distribution, "%i\t%i\t%i\t%f\n", stem_length, loop_length, num_inf_bases, I );
      }

      if( I < minI || I > maxI ) {
	// out of range. the new motif should be deleted
	free(motifs[num_motifs]->phrases) ;
	free(motifs[num_motifs]) ;
      }else{
	num_motifs ++;
	(*num_total_motifs) ++;
	if (num_motifs % 10000 == 0){
	  printf(".") ; 
	  fflush(stdout) ;
	}

	if( num_motifs >= num_motifs_per_file ) {
	  // the maximum number of motifs per file is reached
	  // write the motifs
	  write_and_release_motifs( output_file, *file_num, motifs, num_motifs );

	  // adjust the file number accordingly
	  (*file_num) ++;
	  num_motifs = 0;
	}
      }
    }

    if( !lcl_next_sequence( sequence, num_phrases ) )
      break;
  }

  free(sequence) ;
  return num_motifs;
}

int create_and_write_motifs( char *output_file, int num_motifs_per_file, s_motif **motifs, int min_stem_length, int max_stem_length, int min_loop_length, int max_loop_length, int min_inf_bases, int max_inf_bases, double minI, double maxI ) {
  // if motifs==NULL, it only counts the number of motifs and does not store or write them
  // generate the file name
  
  char distribution_output [10000];
  sprintf( distribution_output, "%s.dist.txt", output_file );

  int num_total_motifs = 0;
  int num_motifs = 0;
  int file_num = 1;

  int stem_length, loop_length ;
  for( stem_length = min_stem_length; stem_length <= max_stem_length; stem_length ++ ){
    printf("\nSTEM LENGTH = %d\n", stem_length) ;
    for( loop_length = min_loop_length; loop_length <= max_loop_length; loop_length ++ )
      num_motifs = lcl_create_motifs( output_file, &file_num, num_motifs_per_file, motifs, num_motifs, &num_total_motifs, stem_length, loop_length, min_inf_bases, max_inf_bases, minI, maxI, NULL );
  }

  if( num_motifs )
    write_and_release_motifs( output_file, file_num, motifs, num_motifs );
	
  return num_total_motifs;
}

double lcl_calculate_probability( s_motif *motif ) {
  double motif_p = 1.0;
  int i ;
  for (i=0; i<motif->num_phrases ; i++ ) {
    double phrase_p = 0;
    if (motif->phrases[i].base & _U) {
      if( motif->phrases[i].structure == _pair )
	phrase_p += ( 0.25 * 0.5 ); // U can create a base pair with either A or G
      else
	phrase_p += 0.25;
    }
    if( motif->phrases[i].base & _C ) {
      if( motif->phrases[i].structure == _pair )
	phrase_p += ( 0.25 * 0.25 ); // C can create a base pair with G
      else
	phrase_p += 0.25;
    }
    if( motif->phrases[i].base & _A ) {
      if( motif->phrases[i].structure == _pair )
	phrase_p += ( 0.25 * 0.25 ); // A can create a base pair with U
      else
	phrase_p += 0.25;
    }
    if( motif->phrases[i].base & _G ) {
      if( motif->phrases[i].structure == _pair )
	phrase_p += ( 0.25 * 0.5 ); // G can create a base pair with either U or C
      else
	phrase_p += 0.25;
    }
    motif_p *= phrase_p;
  }

  return motif_p;
}

double calcluate_information (s_motif *motif) {
  double motif_p = lcl_calculate_probability( motif );
  if( motif_p > 0 )
    return - ( log(motif_p)/log(2.0) );
  return -1;
}
