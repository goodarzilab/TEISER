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

int main(int argc, char ** argv) {
  int      i,j,k,l,m,n,o,p ;

  char     *seedfile ;
  char     *rna_fastafile ;
  char     *expfile ;
  char     *motifoutfile ;

  int      dooptimize = 1 ;
  int      doonlypositive = 0 ;

  float    *E ;
  int      *E_q ;
  int      *M_q ;

  struct   my_hsearch_data *h_rna_ind ;
  ENTRY    e ;
  ENTRY*   ep ;
  int      hashret =0 ;

  int      quantized = 1 ;
  int      ebins = 0 ;
  int      mbins = 2 ;
  int      divbins = 50 ;
  float*   E_q_bins = 0 ;

  int      shuffle = 1000000 ;
  
  int      rnd_fasta = 0 ;
  float    max_p = 0.0000001 ;
  float    max_z = -100 ;

  float    maxfreq = 0.5;
  float    myfreq;
  float    lastmyfreq;

  s_motif  **motifs ;
  s_motif  **motif_vars ;

  char     **seq_names ;
  int      seq_count ;
  int      t_seq_count ;
  int      motif_count = 0 ;

  seedfile         = get_parameter(argc, argv, "-seedfile") ;
  motifoutfile     = get_parameter(argc, argv, "-motifoutfile") ;

  FILE *fmotif ;
  FILE *fptr = fopen ( seedfile, "rb") ;
  if (!fptr){
    printf("Could not open the seed file: %s\n", seedfile) ;
    exit(0) ;
  }

  motif_count = read_motifs( fptr, &motifs ) ;
  printf("%d seeds were loaded...\n", motif_count) ;
  s_motif *motif = copy_motif(motifs[0]) ;
  fflush(stdout) ;
  fclose(fptr) ;

  print_cfg(motif) ;
  getchar() ;

  int nvars=0 ;
  motif_vars = (s_motif**) malloc (100*sizeof(s_motif*)) ;
  
  //motif_vars[nvars]=copy_motif(motif) ;
  motif_vars[nvars] = (s_motif*) malloc (sizeof(s_motif)) ;    
  motif_vars[nvars]->num_phrases = 13 ;
  motif_vars[nvars]->linear_length = 17 ;
  motif_vars[nvars]->phrases = (s_phrase*) malloc (motif_vars[nvars]->num_phrases*sizeof(s_phrase)) ;

  motif_vars[nvars]->phrases[0].structure = _pair ;
  motif_vars[nvars]->phrases[1].structure = _pair ;
  motif_vars[nvars]->phrases[2].structure = _pair ;
  motif_vars[nvars]->phrases[3].structure = _pair ;
  motif_vars[nvars]->phrases[4].structure = _leftBulge ;
  motif_vars[nvars]->phrases[5].structure = _leftBulge ;
  motif_vars[nvars]->phrases[6].structure = _leftBulge ;
  motif_vars[nvars]->phrases[7].structure = _leftBulge ;
  motif_vars[nvars]->phrases[8].structure = _leftBulge ;
  motif_vars[nvars]->phrases[9].structure = _leftBulge ;
  motif_vars[nvars]->phrases[10].structure = _leftBulge ;
  motif_vars[nvars]->phrases[11].structure = _leftBulge ;
  motif_vars[nvars]->phrases[12].structure = _leftBulge ;
  

  motif_vars[nvars]->phrases[0].base = _N ;
  motif_vars[nvars]->phrases[1].base = _N ;
  motif_vars[nvars]->phrases[2].base = _S ;
  motif_vars[nvars]->phrases[3].base = _N ;
  motif_vars[nvars]->phrases[4].base = _C ;
  motif_vars[nvars]->phrases[5].base = _S ;
  motif_vars[nvars]->phrases[6].base = _N ;
  motif_vars[nvars]->phrases[7].base = _N ;
  motif_vars[nvars]->phrases[8].base = _N ;
  motif_vars[nvars]->phrases[9].base = _S ;
  motif_vars[nvars]->phrases[10].base = _S ;
  motif_vars[nvars]->phrases[11].base = _S ;
  motif_vars[nvars]->phrases[12].base = _N ;
  //motif_vars[nvars]->phrases[9].base = _S ;
  print_cfg(motif_vars[nvars]) ;
  nvars++ ;
/*
  motif_vars[nvars]=copy_motif(motif) ;
  motif_vars[nvars]->phrases[0].base = _C ;
  nvars++ ;

  motif_vars[nvars]=copy_motif(motif) ;
  motif_vars[nvars]->phrases[0].base = _U ;
  nvars++ ;

  motif_vars[nvars]=copy_motif(motif) ;
  motif_vars[nvars]->phrases[1].base = _A ;
  nvars++ ;

  motif_vars[nvars]=copy_motif(motif) ;
  motif_vars[nvars]->phrases[1].base = _C ;
  nvars++ ;

  motif_vars[nvars]=copy_motif(motif) ;
  motif_vars[nvars]->phrases[1].base = _U ;
  nvars++ ;

  motif_vars[nvars]=copy_motif(motif) ;
  motif_vars[nvars]->phrases[2].base = _A ;
  nvars++ ;

  motif_vars[nvars]=copy_motif(motif) ;
  motif_vars[nvars]->phrases[2].base = _C ;
  nvars++ ;

  motif_vars[nvars]=copy_motif(motif) ;
  motif_vars[nvars]->phrases[2].base = _U ;
  nvars++ ;

  motif_vars[nvars++]=copy_motif(motif) ;

  */
  //Deletions in loop from phrase 6-9
  //Single loop deletions
  /*for (i=4;i<11;i++){
    s_motif *n_motif = (s_motif*) malloc (sizeof(s_motif)) ;    
    n_motif->num_phrases = motif->num_phrases-1 ;
    n_motif->linear_length = motif->linear_length-1 ;
    n_motif->phrases = (s_phrase*) malloc (n_motif->num_phrases*sizeof(s_phrase)) ;

    int cnt=0 ;
    for (j=0 ; j<motif->num_phrases ; j++){
      if (j==i) continue ;
      n_motif->phrases[cnt].base = motif->phrases[j].base ;
      n_motif->phrases[cnt].structure = motif->phrases[j].structure ;
      cnt++;
        
    }
    motif_vars[nvars++]=copy_motif(n_motif) ;
    free(n_motif->phrases) ;
    free(n_motif) ;
  }

  //Double loop deletions
  for (i=4;i<9;i++){
    s_motif *n_motif = (s_motif*) malloc (sizeof(s_motif)) ;    
    n_motif->num_phrases = motif->num_phrases-2 ;
    n_motif->linear_length = motif->linear_length-2 ;
    n_motif->phrases = (s_phrase*) malloc (n_motif->num_phrases*sizeof(s_phrase)) ;
    int cnt=0 ;
    for (j=0 ; j<motif->num_phrases ; j++){
      if (j==i) continue ;
      if (j==i+1) continue ;
     n_motif->phrases[cnt].base = motif->phrases[j].base ;
     n_motif->phrases[cnt].structure = motif->phrases[j].structure ;
     cnt++;
    }
    motif_vars[nvars++]=copy_motif(n_motif) ;
    free(n_motif->phrases) ;
    free(n_motif) ;
  }

  //triple loop deletions
  for (i=4;i<8;i++){
    s_motif *n_motif = (s_motif*) malloc (sizeof(s_motif)) ;    
    n_motif->num_phrases = motif->num_phrases-3 ;
    n_motif->linear_length = motif->linear_length-3 ;
    n_motif->phrases = (s_phrase*) malloc (n_motif->num_phrases*sizeof(s_phrase)) ;
    int cnt=0 ;
    for (j=0 ; j<motif->num_phrases ; j++){
      if (j==i) continue ;
      if (j==i+1) continue ;
      if (j==i+2) continue ;
     n_motif->phrases[cnt].base = motif->phrases[j].base ;
     n_motif->phrases[cnt].structure = motif->phrases[j].structure ;
     cnt++;
    }
    motif_vars[nvars++]=copy_motif(n_motif) ;
    free(n_motif->phrases) ;
    free(n_motif) ;
  }

  //quad loop deletions
  for (i=4;i<7;i++){
    s_motif *n_motif = (s_motif*) malloc (sizeof(s_motif)) ;    
    n_motif->num_phrases = motif->num_phrases-4 ;
    n_motif->linear_length = motif->linear_length-4 ;
    n_motif->phrases = (s_phrase*) malloc (n_motif->num_phrases*sizeof(s_phrase)) ;
    int cnt=0 ;
    for (j=0 ; j<motif->num_phrases ; j++){
      if (j==i) continue ;
      if (j==i+1) continue ;
      if (j==i+2) continue ;
      if (j==i+3) continue ;
     n_motif->phrases[cnt].base = motif->phrases[j].base ;
     n_motif->phrases[cnt].structure = motif->phrases[j].structure ;
     cnt++;
    }
    motif_vars[nvars++]=copy_motif(n_motif) ;
    free(n_motif->phrases) ;
    free(n_motif) ;
  }*/

  //Insertions in the loop (1-3) at three positions
  /*int ncnt=nvars ;
  for(i=0 ;i<ncnt;i++){
    s_motif *v_motif=copy_motif(motif_vars[i]) ;
    
    //print_cfg(v_motif) ;
    //printf("HERE") ;getchar() ;
    for (j=0;j<v_motif->num_phrases;j++){
      if (v_motif->phrases[j].structure==_leftBulge && v_motif->phrases[j].base!=_N){
        int anchor=j ;
        //printf("anchor=%i",anchor) ; getchar() ;
        for(k=1 ;k<3 ;k++){
          s_motif *n_motif = (s_motif*) malloc (sizeof(s_motif)) ;    
          n_motif->num_phrases = v_motif->num_phrases+k ;
          n_motif->linear_length = v_motif->linear_length+k ;
          n_motif->phrases = (s_phrase*) malloc (n_motif->num_phrases*sizeof(s_phrase)) ;
          int cnt=0 ;
          for (l=0 ; l<v_motif->num_phrases ; l++){
            if (l==anchor){
              for(m=0 ;m<k ;m++){
                n_motif->phrases[cnt].base = _N ;
                n_motif->phrases[cnt].structure = _leftBulge ;
                cnt++ ;
              }
            }
            n_motif->phrases[cnt].base = v_motif->phrases[l].base ;
            n_motif->phrases[cnt].structure = v_motif->phrases[l].structure ;
            cnt++;
          }
          //printf("HERE: %i: %i-%i\n", i, anchor, k) ;
          //fflush(stdout) ;
          //print_cfg(v_motif) ;getchar() ;
          //print_cfg(n_motif) ;getchar() ;
          motif_vars[nvars++]=copy_motif(n_motif) ;
          free(n_motif->phrases) ;
          free(n_motif) ; 
        }
      }
    }
    free(v_motif->phrases) ;
    free(v_motif) ; 
  }*/

  //Bulges (left & right) in the ste, (1-3) 
  /*ncnt=nvars ;
  for(i=0 ;i<ncnt;i++){
    s_motif *v_motif=copy_motif(motif_vars[i]) ;
    
    //print_cfg(v_motif) ;
    //printf("HERE") ;getchar() ;
    for (j=0;j<v_motif->num_phrases;j++){
      if (v_motif->phrases[j].structure==_pair){
        int anchor=j ;
        //printf("anchor=%i",anchor) ; getchar() ;
        for(k=1 ;k<2 ;k++){
          s_motif *n_motif = (s_motif*) malloc (sizeof(s_motif)) ;    
          n_motif->num_phrases = v_motif->num_phrases+k ;
          n_motif->linear_length = v_motif->linear_length+k ;
          n_motif->phrases = (s_phrase*) malloc (n_motif->num_phrases*sizeof(s_phrase)) ;
          int cnt=0 ;
          for (l=0 ; l<v_motif->num_phrases ; l++){
            if (l==anchor){
              for(m=0 ;m<k ;m++){
                n_motif->phrases[cnt].base = _N ;
                n_motif->phrases[cnt].structure = _leftBulge ;
                cnt++ ;
              }
            }
            n_motif->phrases[cnt].base = v_motif->phrases[l].base ;
            n_motif->phrases[cnt].structure = v_motif->phrases[l].structure ;
            cnt++;
          }
          //printf("HERE: %i: %i-%i\n", i, anchor, k) ;
          //fflush(stdout) ;
          //print_cfg(v_motif) ;getchar() ;
          //print_cfg(n_motif) ;getchar() ;
          motif_vars[nvars++]=copy_motif(n_motif) ;
          free(n_motif->phrases) ;
          free(n_motif) ; 
        }
        for(k=1 ;k<2 ;k++){
          s_motif *n_motif = (s_motif*) malloc (sizeof(s_motif)) ;    
          n_motif->num_phrases = v_motif->num_phrases+k ;
          n_motif->linear_length = v_motif->linear_length+k ;
          n_motif->phrases = (s_phrase*) malloc (n_motif->num_phrases*sizeof(s_phrase)) ;
          int cnt=0 ;
          for (l=0 ; l<v_motif->num_phrases ; l++){
            if (l==anchor){
              for(m=0 ;m<k ;m++){
                n_motif->phrases[cnt].base = _N ;
                n_motif->phrases[cnt].structure = _rightBulge ;
                cnt++ ;
              }
            }
            n_motif->phrases[cnt].base = v_motif->phrases[l].base ;
            n_motif->phrases[cnt].structure = v_motif->phrases[l].structure ;
            cnt++;
          }
          //printf("HERE: %i: %i-%i\n", i, anchor, k) ;
          //fflush(stdout) ;
          //print_cfg(v_motif) ;getchar() ;
          //print_cfg(n_motif) ;getchar() ;
          motif_vars[nvars++]=copy_motif(n_motif) ;
          free(n_motif->phrases) ;
          free(n_motif) ; 
        }
      }
    }
    free(v_motif->phrases) ;
    free(v_motif) ; 
  }

  //extend stem from 1-4
  ncnt=nvars ;
  for(i=0 ;i<ncnt;i++){
    s_motif *v_motif=copy_motif(motif_vars[i]) ;
    
    //print_cfg(v_motif) ;
    //printf("HERE") ;getchar() ;
    for (j=0;j<v_motif->num_phrases;j++){
        for(k=1 ;k<3 ;k++){
          s_motif *n_motif = (s_motif*) malloc (sizeof(s_motif)) ;    
          n_motif->num_phrases = v_motif->num_phrases+k ;
          n_motif->linear_length = v_motif->linear_length+k*2 ;
          n_motif->phrases = (s_phrase*) malloc (n_motif->num_phrases*sizeof(s_phrase)) ;
          int cnt=0 ;
          for(m=0 ;m<k ;m++){
                n_motif->phrases[cnt].base = _N ;
                n_motif->phrases[cnt].structure = _pair ;
                cnt++ ;
          }
          for (l=0 ; l<v_motif->num_phrases ; l++){
            n_motif->phrases[cnt].base = v_motif->phrases[l].base ;
            n_motif->phrases[cnt].structure = v_motif->phrases[l].structure ;
            cnt++;
          }
          //printf("HERE: %i: %i-%i\n", i, anchor, k) ;
          //fflush(stdout) ;
          //print_cfg(v_motif) ;getchar() ;
          //print_cfg(n_motif) ;getchar() ;
          motif_vars[nvars++]=copy_motif(n_motif) ;
          free(n_motif->phrases) ;
          free(n_motif) ; 
        }
      }
    free(v_motif->phrases) ;
    free(v_motif) ; 
  }*/

  /*int C=3 ;
  int D=5 ;
  for (i=0 ; i<C ; i++){
    n_motif->phrases[i].base = motif->phrases[i].base ;
    n_motif->phrases[i].structure = motif->phrases[i].structure ;
  }
  for ( i=motif->num_phrases ; i>D ; i--){
    n_motif->phrases[i+3].base = motif->phrases[i].base ;
    n_motif->phrases[i+3].structure = motif->phrases[i].structure ;
    }*/
  //n_motif->phrases[3].base = _N ;
  /*n_motif->phrases[C].base = _N ;
  n_motif->phrases[C].structure = _leftBulge ;
  n_motif->phrases[C+1].base = _N ;
  n_motif->phrases[C+1].structure = _leftBulge ;
  n_motif->phrases[C+2].base = _N ;
  n_motif->phrases[C+2].structure = _leftBulge ;
  n_motif->phrases[C+3].base = _N ;
  n_motif->phrases[C+3].structure = _rightBulge ;
  n_motif->phrases[C+4].base = _N ;
  n_motif->phrases[C+4].structure = _rightBulge ;
  n_motif->phrases[C+5].base = _N ;
  n_motif->phrases[C+5].structure = _rightBulge ;*/

  //fmotif = fopen(motifoutfile, "wb") ;
  //if (!fmotif)
    //die("Cannot open datafile\n");

  //fwrite( 1, sizeof(int), 1, fptr );
  //lcl_write_motif(fptr, n_motif) ;
  /*printf("\n%i motifs", nvars) ; getchar() ;
  for (i=0;i<nvars;i++){
    printf("%i:%i\n", i,motif_vars[i]->num_phrases);fflush(stdout) ;
    print_cfg(motif_vars[i]) ;
  }*/
  printf("%i seeds written...\n", nvars) ;
  fmotif = fopen(motifoutfile, "wb") ;
  write_motifs (fmotif, motif_vars, nvars) ;
  fclose(fmotif) ;
  return (0) ;
}
