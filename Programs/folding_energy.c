#include "stdlib.h"
#include "stdio.h"

#include "structures.h"
#include "nucleotides.h"
#include "folding_energy.h"

/////////////////////////////////////////////////////////////////////////////////////////
double lcl_get_pairing_energy( NUCBIT base1, NUCBIT base2 )
// returns the energy associated with pairing of base1 with base2
// if the two bases cannot pair, returns +10000
{
	if( ( base1 == _C && base2 == _G ) ||
		( base1 == _G && base2 == _C ) )
		return -3;
	
	if( ( base1 == _U && base2 == _A ) ||
		( base1 == _A && base2 == _U ) )
		return -2;
		
	if( ( base1 == _U && base2 == _G ) ||
		( base1 == _G && base2 == _U ) )
		return -1;

	return 10000;
}

double lcl_get_free_folding_energy( s_motif *motif, NUCBIT *sequence )
// returns the free energy for the "sequence" taking the fold that is determined by the "motif"
// returns +10000 if an unknown CFG phrase is encountered
{
	int i ;
	// set the initial left and right bases
	NUCBIT *left = sequence;
	NUCBIT *right = sequence + motif ->linear_length - 1;
	
	double energy = 0;

	for( i = 0; i < motif ->num_phrases; i ++ )
		switch( motif ->phrases[ i ].structure ){
		case _pair:
			energy += lcl_get_pairing_energy( *left, *right );
			left ++;
			right --;
			break;
		case _leftBulge:
			left ++;
			break;
		case _rightBulge:
			right --;
			break;
		default:
			return 10000;
	}
	return energy;
}

void lcl_get_surrounding_sequence(NUCBIT *source_seq, int seq_length, int exclude_start_pos, int exclude_length, NUCBIT *surrounding_seq, int *surrounding_length )
// excludes the specified range from the source sequence, replaces it with four N's,
//  and copies the result into surrounding_seq
{
	int replaced = 0;
	*surrounding_length = 0;
	
	int i;
	for( i = 0; i < seq_length; i ++ )
		if( i < exclude_start_pos || i >= exclude_start_pos + exclude_length ){
			surrounding_seq[ *surrounding_length ] = source_seq[ i ];
			(*surrounding_length) ++;
		}else if( !replaced ){ // the exclude region is encountered, but it hasn't been replaced with N's yet
			int j;
			for( j = 0; j < 4; j ++ ){
				surrounding_seq[ *surrounding_length ] = _N;
				(*surrounding_length) ++;
			}
		}				
}

double lcl_get_minimum_energy( NUCBIT *sequence, int seq_length )
{
	static int max_length = 0; // the maximum length that has been passed to this function so far
	static double *matrix[ MAX_SEQ_LENGTH ]; // the matrix for dynamic programming

	if( seq_length <= 0 )
		return 0; // the sequence length is zero, and therefore there is nothing to do
		
	if( seq_length > max_length ){ // this sequence length is larger than any length before
	// therefore, the matrix memory needs to be adjusted
		// release the previously allocated memory
		int i;
		for( i = 0; i < max_length; i ++ )
			free(matrix[i]) ;
			
		// allocate new memory
		max_length = seq_length;
		for( i = 0; i < max_length; i ++ )
			matrix[i] = (double*) malloc (max_length*sizeof(double)) ;
	}
	
	// now calculate the matrix values
	// start from bottom-right
	int start_row;
	for( start_row = seq_length - 1; start_row >= 0; start_row -- ){
		int row, col;
		// start from start_row and start_col, and go diagonally to left-top of the matrix
		for( row = start_row, col = seq_length - 1; row >= 0 && col >= 0; row --, col -- ){
			if( col - row < 4 ) // at least a three-base spacing is needed for the bases to pair
				matrix[ row ][ col ] = 0;
			else{
				double min = matrix[ row + 1 ][ col ];

				if( matrix[ row ][ col - 1 ] < min )
					min = matrix[ row ][ col - 1 ];

				if( (
						matrix[ row + 1 ][ col - 1 ] +
						lcl_get_pairing_energy( sequence[ row ], sequence[ col ] )
					) < min )
					min =
						matrix[ row + 1 ][ col - 1 ] +
						lcl_get_pairing_energy( sequence[ row ], sequence[ col ] );
				
				int k;
				for( k = row + 1; k < col; k ++ )
					if( (
							matrix[ row ][ k ] + matrix[ k + 1 ][ col ]
						) < min )
						min = matrix[ row ][ k ] + matrix[ k + 1 ][ col ];
						
				matrix[ row ][ col ] = min;
			}
		}
	}
	
	return matrix[0][seq_length - 1];
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////  This is the function that should be called by TEISER /////////////////
/////////////////////////////////////////////////////////////////////////////////////////
double get_dG_ratio_to_optimum( s_motif *motif, NUCBIT *sequence, int seq_length, int motif_start_pos )
// returns the ratio of the dG for a fold that takes the shape of the "motif" at the specified position
//  relative to the global minimum dG for folding
// returns a negative value if unsuccessful
//
// motif: contains the structure of the motif
// sequence: the sequence that contains an instance of the motif + some of the surrounding sequence
// seq_length: the length of "sequence"
// motif_start_pos: the start position of the motif instance (0-based index, relative to the start of "sequence")
{
	// get the sequence that surrounds the motif instance
	NUCBIT surrounding_seq[ MAX_SEQ_LENGTH ];
	int surrounding_length;
	lcl_get_surrounding_sequence(
		sequence, seq_length,
		motif_start_pos, motif ->linear_length,
		surrounding_seq,
		&surrounding_length );
	
	// the folding dG of the whole sequence without any constraint
	double min_dG =
		lcl_get_minimum_energy( sequence, seq_length );
		
	if( min_dG >= 0 )
		return -1;
		
	// the folding dG with the constraint of having the motif instance in the specified position
	// This is equal to the folding dG of the motif instance + the folding dG of the surrounding sequence
	double this_dG =
		lcl_get_free_folding_energy( motif, sequence + motif_start_pos ) +
		lcl_get_minimum_energy( surrounding_seq, surrounding_length );
		
	return this_dG / min_dG;
}
