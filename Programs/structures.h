#include <string.h>

#define MAX_MOTIFS		10000000
#define MAX_SEQUENCES		30000
#define MAX_SEQ_LENGTH		1000

typedef unsigned char BYTE;
typedef BYTE NUCBIT;

typedef struct _s_phrase	{
  // an RNA molecule is shown by a string of s_phrase
  NUCBIT base;      // _U | _C | _G | _A; if structure = _pair, base shows the left nucleotide in the paired bases
  BYTE structure;  // _pair ^ _leftBulge ^ _rightBulge
} s_phrase;

typedef struct _s_sequence {
  char *name;

  NUCBIT *bases;
  int length;
} s_sequence ;

typedef struct _s_motif {
  s_phrase *phrases;
  int num_phrases;

  int linear_length;
} s_motif ;
