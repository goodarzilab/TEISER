#define _U				0x01
#define _C				0x02
#define _G				0x04
#define _A				0x08
#define _N				0x0F

#define _Y				0x03 //UC
#define _R				0x0C //AG
#define _K				0x05 //UG
#define _M				0x0A //AC
#define _S				0x06 //GC
#define _W				0x09 //AU

#define _B				0x07 //GUC
#define _D				0x0D //GAU
#define _H				0x0B //ACU
#define _V				0x0E //GCA

#define NUM_LETTERS		4
#define LAST_LETTER		_A

static char lcl_base_templates1[] = "TCGA";
static char lcl_base_templates2[] = "tcga";
static char lcl_base_templates3[] = "UCGA";
static char lcl_base_templates4[] = "ucga";

static char lcl_N1 = 'N';
static char lcl_N2 = 'n';

#define _pair			1
#define _leftBulge		2
#define _rightBulge		3

/*
NUCBIT	IUPAC	IUPAC_RC
1	_U	_R
2	_C	_G
3	_Y	_R
4	_G	_Y
5	_K	_N
6	_S	_B
7	_B	_N
8	_A	_U
9	_W	_D
10	_M	_K
11	_H	_D
12	_R	_Y
13	_D	_N
14	_V	_B
15	_N	_N
*/
static char nucbit_map[] = {'U', 'C', 'Y', 'G', 'K', 'S', 'B', 'A', 'W', 'M', 'H', 'R', 'D', 'V', 'N'};
static char nucbit_map_rc[] = {'R', 'G', 'R', 'Y', 'N', 'B', 'N', 'U', 'D', 'K', 'D', 'Y', 'N', 'B', 'N'};
