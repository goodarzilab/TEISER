#include "stdlib.h"
#include "stdio.h"

#include "structures.h"
#include "nucleotides.h"
#include "matchmaker.h"
#include "folding_energy.h"

void get_linear_length( s_motif *motif ) {
  // sets the value for the linear length of a motif
  motif->linear_length = 0;
  int i=0 ;
  for( i=0 ; i<motif->num_phrases ; i++ )
    if ( motif->phrases[i].structure == _pair )
      motif->linear_length += 2; // each base pair is composed of two bases
    else
      motif->linear_length ++; // other structures are composed of one base each
}

int lcl_is_a ( NUCBIT child, NUCBIT parent ) {
// returns true if child is an instance of parent, false otherwise
  if( (child & parent) != child ) // each 'on' bit in child should also be 'on' in parent
    return 0;
  return 1;
}

int lcl_pair_bases( NUCBIT base1, NUCBIT base2 ) {
// returns true if the two bases can form a Watson-Crick (including G-U) base pair, false otherwise
  if( ( base1 == _U && base2 == _A ) ||
      ( base1 == _C && base2 == _G ) ||
      ( base1 == _G && base2 == _C ) ||
      ( base1 == _A && base2 == _U ) ||
      ( base1 == _U && base2 == _G ) ||
      ( base1 == _G && base2 == _U ) )
    return 1;

  return 0;
}

NUCBIT lcl_complement_base( NUCBIT base ) {
  NUCBIT complement = 0;
  
  if( lcl_is_a( _U, base ) )
    {
      complement |= _A;
      complement |= _G;
    }

  if( lcl_is_a( _A, base ) )
    complement |= _U;

  if( lcl_is_a( _C, base ) )
    complement |= _G;

  if( lcl_is_a( _G, base ) )
    {
      complement |= _C;
      complement |= _U;
    }

  return complement;
}

int lcl_match_motif_seq( s_motif *motif, NUCBIT *sequence, int b, int l, float dG_t ) {
// returns true if the motif can be matched with the sequence
  // set the initial left and right bases
  NUCBIT *left = sequence;
  NUCBIT *right = sequence + motif->linear_length-1;

  int i ;
  for( i=0; i<motif->num_phrases ; i++ )
    switch( motif->phrases[i].structure ) {
    case _pair:
      //if( !lcl_is_a ( *left, motif->phrases[i].base ) || !lcl_is_a( *right, lcl_complement_base( motif->phrases[i].base ) ) )
      if( !lcl_is_a ( *left, motif->phrases[i].base ) || !lcl_pair_bases( *left, *right ) )
	// either the sequence does not match or the left and right bases cannot form a Watson-Crick base pair
	     return 0;
      left ++;
      right --;
      break;
    case _leftBulge:
      if( !lcl_is_a( *left, motif->phrases[i].base ) ) // the sequence does not match
	     return 0;
      left ++;
      break;
    case _rightBulge:
      if( !lcl_is_a( *right, motif->phrases[i].base ) ) // the sequence does not match
	return 0;
      right --;
      break;
    default:
      return 0;
  }
	
  if (dG_t > 0){
    double r = get_dG_ratio_to_optimum( motif, sequence-b, l, 0 ) ;
    if (r<dG_t)
      return 0 ;
  }
  return -1;
}

int find_motif_instance ( s_motif *motif, s_sequence *sequence, float dG_t ) {
  // returns true if an instance of the motif is found in the sequence, false otherwise
  int i=0 ;
  for( i=0 ; i<sequence->length-motif->linear_length+1 ; i++ ){
    int b = 0 ;
    if (i<b){
      b=i ;
    }
    int l = motif->linear_length+0 ;
    if (i+l>=sequence->length){
      l -= (i+l-sequence->length+1) ;
    }
    
    if( lcl_match_motif_seq( motif, sequence->bases + i, b,l, dG_t ) )
      return i;
  }
  return -1;
}

int* find_motif_instance_positions ( s_motif *motif, s_sequence *sequence, int *hits, float dG_t) {
  // returns true if an instance of the motif is found in the sequence, false otherwise
  int i=0 ;
  *hits =0 ;

  int *temp = (int*) malloc (1000*sizeof(int)) ;
  for( i=0 ; i<sequence->length-motif->linear_length+1 ; i++ ){
    int b = 0 ;
    if (i<b){
      b=i ;
    }
    int l = motif->linear_length+0 ;
    if (i+l>=sequence->length){
      l -= (i+l-sequence->length+1) ;
    }
    if( lcl_match_motif_seq( motif, sequence->bases + i, b,l,dG_t) ){
      temp[*hits] = i ;
      (*hits)++ ;
    }
  }

  
  if (*hits > 0){
    int *m = (int*) malloc (sizeof(int)*(*hits)) ;
    for (i=0 ; i<*hits ; i++){
      m[i] = temp[i] ;
    }
    free(temp) ;
    return m ;
  }
  free(temp) ;
  return NULL ;
}

int lcl_match_motif_only_seq( s_motif *motif, NUCBIT *sequence ) {
  // returns true if the motif, regardless of structure, can be matched with the sequence
  // set the initial left and right bases
  NUCBIT *left = sequence;
  NUCBIT *right = sequence + motif ->linear_length - 1;
  int i ;
  for( i=0; i<motif->num_phrases ; i++ )
    switch( motif->phrases[i].structure ) {
    case _pair:
      if( !lcl_is_a ( *left, motif->phrases[i].base ))
	// either the sequence does not match or the left and right bases cannot form a Watson-Crick base pair
	return 0;
      left ++;
      right --;
      break;
    case _leftBulge:
      if( !lcl_is_a( *left, motif->phrases[i].base ) ) // the sequence does not match
	return 0;
      left ++;
      break;
    case _rightBulge:
      if( !lcl_is_a( *right, motif->phrases[i].base ) ) // the sequence does not match
	return 0;
      right --;
      break;
    default:
      return 0;
    }
	
  return -1;
}


int find_motif_instance_seq_only( s_motif *motif, s_sequence *sequence ) {
  // returns true if an instance of the motif is found in the sequence, false otherwise
  int i=0 ;
  for( i=0 ; i<sequence->length-motif->linear_length+1 ; i++ )
    if( lcl_match_motif_only_seq( motif, sequence->bases + i ) )
      return i;
  return -1;
}
