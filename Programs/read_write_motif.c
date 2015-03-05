#include <stdio.h>
#include <stdlib.h>

#include "structures.h"
#include "nucleotides.h"
#include "read_write_motif.h"

void lcl_write_motif( FILE *fptr, s_motif *motif ) {
  // writes the motif in fptr
  fwrite( &(motif->num_phrases), sizeof(int), 1, fptr );
  fwrite( motif->phrases, sizeof(s_phrase), motif->num_phrases, fptr );
  fwrite( &(motif->linear_length), sizeof(int), 1, fptr );
}

void lcl_write_motifs( FILE *fptr, s_motif **motifs, int num_motifs ) {
  // writes all motifs of the specified array in fptr. the file starts with the number of motifs, followed by motif data
  fwrite( &num_motifs, sizeof(int), 1, fptr );

  int i ;
  for (i=0 ; i<num_motifs ; i++ )
    lcl_write_motif( fptr, motifs[i] );
}


void write_motifs( FILE *fptr, s_motif **motifs, int num_motifs ) {
// writes all motifs of the specified array in fptr. the file starts with the number of motifs, followed by motif data
  fwrite( &num_motifs, sizeof(int), 1, fptr );
  int i=0 ;
  printf("writing motifs to file...") ;
  for( i=0; i<num_motifs; i++ ){
    lcl_write_motif( fptr, motifs[i] );
  }
  printf("\n") ;
}

void write_and_release_motifs( char *output_file, int file_num, s_motif **motifs, int num_motifs ) {
  int i=0 ;
  // generate the file name
  char this_output[1000];
  sprintf( this_output, "%s.%3.3i.bin", output_file, file_num );
  
  // opening the output file
  FILE *fptr_out = fopen( this_output, "wb" );
  
  if (!fptr_out) {
    printf("Warning: cannot open the seed output file: %s\n", this_output) ;
    return;
  }

  // write the output
  printf("\nWriting motif set no. %d...\n", file_num) ;
  lcl_write_motifs( fptr_out, motifs, num_motifs );

  // delete the written motifs
  for( i = 0; i < num_motifs; i ++ ) {
    free (motifs[i]->phrases) ;
    free (motifs[i]);
  }
  fclose( fptr_out );
}


void lcl_read_motif( FILE *fptr, s_motif *motif ) {
// reads a single motif from fptr
  fread( &(motif->num_phrases), sizeof(int), 1, fptr );

  motif->phrases = (s_phrase*) malloc (motif->num_phrases*sizeof(s_phrase)) ;
  fread ( motif->phrases, sizeof(s_phrase), motif->num_phrases, fptr ) ;

  fread( &(motif->linear_length), sizeof(int), 1, fptr );
}

int read_motifs( FILE *fptr, s_motif ***motifs ) {
// reads all the motifs in fptr and returns the number of read motifs
// allocates the necessary memory to "motifs"

  int num_motifs = 0;
  fread( &num_motifs, sizeof(int), 1, fptr );

  *motifs = (s_motif**) malloc (num_motifs*sizeof(s_motif*)) ;
  
  int i=0 ;
  for( i=0 ; i<num_motifs; i++ ) {
    (*motifs)[i] = (s_motif*) malloc (sizeof(s_motif)) ;
    lcl_read_motif( fptr, (*motifs)[i] );
  }

  return num_motifs;
}

s_motif* copy_motif (s_motif *m1){
	int i=0 ;
	s_motif *m2 = (s_motif*) malloc (sizeof(s_motif)) ;
	m2->num_phrases = m1->num_phrases ;
	m2->linear_length = m1->linear_length ;
	m2->phrases = (s_phrase*) malloc (m1->num_phrases*sizeof(s_phrase)) ;
	for (i=0 ; i<m1->num_phrases ; i++){
		m2->phrases[i].base = m1->phrases[i].base ;
		m2->phrases[i].structure = m1->phrases[i].structure ;
	}
	return m2 ;
}

void print_cfg(s_motif *motif) {
  int length = motif->linear_length ;
  int i=0, b=0, e=length-1 ;
  for( i=0; i<motif->num_phrases; i++ ){
    printf("phrase %d: ", i) ;
    if (motif->phrases[i].structure == _leftBulge){
      printf("leftBulge\t") ;
    }else if (motif->phrases[i].structure == _rightBulge){
      printf("rightBulge\t") ;
    }else if (motif->phrases[i].structure == _pair){
      printf("basePair\t") ;
    }
    printf("_%c\n",nucbit_map[(motif->phrases[i].base)-1]) ;
  }
}

void print_motif(s_motif *motif) {
  int length = motif->linear_length ;
  int i=0, b=0, e=length-1 ;
  char *temp = (char*) malloc ((length+1)*sizeof(char)) ;
  char *str = (char*) malloc ((length+1)*sizeof(char)) ;
  temp[length] = '\0' ;
  str[length] = '\0' ;
  for( i=0; i<motif->num_phrases; i++ ){
    if (motif->phrases[i].structure == _leftBulge){
      str[b] = '.' ;
    }else if (motif->phrases[i].structure == _rightBulge){
      str[e] = '.' ;
    }else if (motif->phrases[i].structure == _pair){
      str[b] = '(' ;
      str[e] = ')' ;
    }
    if (motif->phrases[i].structure == _leftBulge){
      temp[b++] = nucbit_map[(motif->phrases[i].base)-1] ;
    }else if (motif->phrases[i].structure == _rightBulge){
      temp[e--] = nucbit_map[(motif->phrases[i].base)-1] ;
    }else if (motif->phrases[i].structure == _pair){
      temp[b++] = nucbit_map[(motif->phrases[i].base)-1] ;
      temp[e--] = nucbit_map_rc[(motif->phrases[i].base)-1] ;
    }

  }
  
  printf("%s\n%s", temp, str) ;
  free(temp) ;
  free(str) ;
}

char* print_motif_to_char(s_motif *motif) {
  
  int length = motif->linear_length ;
  int i=0, b=0, e=length-1 ;
  char *temp = (char*) malloc ((length+1)*sizeof(char)) ;
  char *str = (char*) malloc ((length+1)*sizeof(char)) ;
  temp[length] = '\0' ;
  str[length] = '\0' ;
	for( i=0; i<motif->num_phrases; i++ ){
		if (motif->phrases[i].structure == _leftBulge){
			str[b] = '.' ;
		}else if (motif->phrases[i].structure == _rightBulge){
			str[e] = '.' ;
		}else if (motif->phrases[i].structure == _pair){
			str[b] = '(' ;
			str[e] = ')' ;
		}
        if (motif->phrases[i].structure == _leftBulge){
            temp[b++] = nucbit_map[(motif->phrases[i].base)-1] ;
        }else if (motif->phrases[i].structure == _rightBulge){
            temp[e--] = nucbit_map[(motif->phrases[i].base)-1] ;
        }else if (motif->phrases[i].structure == _pair){
            temp[b++] = nucbit_map[(motif->phrases[i].base)-1] ;
            temp[e--] = nucbit_map_rc[(motif->phrases[i].base)-1] ;
        }
	}
	
  char *m = (char*) malloc (length*3*sizeof(char)) ;
  sprintf(m, "%s\t%s", temp, str) ;
  free(temp) ;
  free(str) ;
  return m ;
}

char* print_motif_to_seq(s_motif *motif) {
  int length = motif->linear_length ;
  int i=0, b=0, e=length-1 ;
  char *temp = (char*) malloc ((length+1)*sizeof(char)) ;
  temp[length] = '\0' ;
	for( i=0; i<motif->num_phrases; i++ ){
        if (motif->phrases[i].structure == _leftBulge){
            temp[b++] = nucbit_map[(motif->phrases[i].base)-1] ;
        }else if (motif->phrases[i].structure == _rightBulge){
            temp[e--] = nucbit_map[(motif->phrases[i].base)-1] ;
        }else if (motif->phrases[i].structure == _pair){
            temp[b++] = nucbit_map[(motif->phrases[i].base)-1] ;
            temp[e--] = nucbit_map_rc[(motif->phrases[i].base)-1] ;
        }
	}
	
  return temp ;
}
