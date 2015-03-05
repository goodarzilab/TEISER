#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <search.h>
#include <limits.h>
#include <time.h>
#include "sys/types.h"
#include "string.h"

#include "nucleotides.h"
#include "structures.h"
#include "sequences.h"
#include "matchmaker.h"
#include "dataio.h"
#include "statistics.h"
#include "hashtable.h"
#include "readFASTA.h"
#include "read_write_motif.h"
#include "information.h"
#include "mi_library.h"
#include "teiser_functions.h"

int* readIntSet(char* setfile, int* n) ;

int main(int argc, char ** argv) {
    int      i ;

    char     *motiffile ;

    s_motif  **motifs ;
    int      motif_count = 0 ;

    motiffile        = get_parameter(argc, argv, "-motiffile") ;
    FILE *f, *fmotif ;
    FILE *fptr = fopen ( motiffile, "rb") ;
    if (!fptr){
    printf("Could not open the seed file: %s\n", motiffile) ;
        exit(0) ;
    }

    motif_count = read_motifs( fptr, &motifs ) ;
    printf("%d seeds were loaded...\n", motif_count) ;
    fflush(stdout) ;
    fclose(fptr) ;

    for (i=0 ; i<motif_count ; i++){
      printf("%d\t%s\n", i, print_motif_to_char(motifs[i])) ;
      // print_motif(motifs[i]) ;
      //printf("\n") ;
    }
    return (0) ;
}
