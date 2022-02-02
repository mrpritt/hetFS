/*
 Program to generate the 120 flow shop instances proposed by
 Taillard, E.D.: "Benchmarks for basic scheduling problems",
 EJOR vol. 64, pp. 278-285, 1993
 
 Originaly written by Taillard in PASCAL, re-written in C by
 Dirk C. Mattfeld, University of Bremen <dirk@uni-bremen.de> and
 Rob J.M. Vaessens, Eindhoven University of Technology <robv@win.tue.nl>.
 Last update : 2/7/1996.

 modified by Alexander Benavides for heterogeneous flowshop for disabled workers.

 compile with$ gcc taHFSPgen.c /path/to/libm.a -o taHFSPgen


     Verification file 'ta001i2' (if VERIFY == 1 and FIRMACIND == 1)
     
     20 5 5
      0  83  99  62  93 inf 1  96  84 133  99 145 2  31 inf inf  21  32 3  99 129 127 inf  74 4  88  74 102 116  62
      0 113 127 114 105 inf 1   4   4   3   5   3 2  93 inf inf 163 173 3  84  71  99 inf  68 4 104  93  83 107  70
      0  27  29  17  19 inf 1  17  11  22  15  21 2  84 inf inf  59  86 3  31  45  41 inf  54 4  36  34  26  24  38
      0 126  94 134 119 inf 1 171 144 195 106 186 2  22 inf inf  27  18 3 136  91  99 inf 134 4 154  88 103 114  95
      0 138  78  84 154 inf 1  95  75  58  68 107 2 135 inf inf 113 160 3 142 140 137 inf 112 4  63  89  66 104  58
      0  62  72  58  37 inf 1 102 103  88  91 100 2  90 inf inf  49  89 3 115 166 127 inf 111 4  58  53  68  52  59
      0  54  85  92  93 inf 1 196 132 134 160 141 2 102 inf inf  83 100 3  23  21  19 inf  20 4  88  89  68  97  62
      0  43  45  71  39 inf 1  82  86  64 111  63 2  45 inf inf  23  35 3  75  92  76 inf 118 4  65  82  79  70  48
      0  31  52  40  31 inf 1   5   5   7   8   8 2  72 inf inf  57  96 3  83  83  95 inf  69 4 108  95  70  95  71
      0 126 139 107 137 inf 1  72 109 103  72  58 2 124 inf inf  76  99 3 143  95 113 inf 146 4  19  15  22  18  21
      0  81  76 124 107 inf 1   3   5   5   4   3 2  10 inf inf   8   7 3 147 107 159 inf 142 4 106 128 120  95 122
      0 134 114 160 161 inf 1 102  87  75  99  99 2   2 inf inf   1   2 3  11  10  13 inf  14 4 132 132  92 108  83
      0  23  22  18  24 inf 1  93 129 122 133 146 2  69 inf inf 105  67 3  66  60  58 inf  56 4  11  12  11  10  11
      0  56  32  46  43 inf 1 101  86  88 133 108 2  48 inf inf  77  42 3  80  53  68 inf  57 4  60  57  83  97  55
      0  23  13  16  17 inf 1  59  83  88  75  77 2  80 inf inf 108  89 3  63 112  85 inf  99 4  64  85  66  73  61
      0 107 118 137 116 inf 1  23  26  27  22  23 2  91 inf inf  91  63 3  69  59  59 inf  69 4 106 153 135 136 167
      0  51  50  38  60 inf 1  29  21  21  30  35 2  36 inf inf  26  49 3  65  99 108 inf  85 4 101 106 115  59  66
      0 117 133 122 147 inf 1 118 128 139  90 115 2 130 inf inf 110 122 3  89 154 127 inf 143 4  18  18  24  27  23
      0 127 117  73 101 inf 1   8  10   9   7   6 2  86 inf inf 132  78 3  94  83  64 inf  55 4  97 106  70  89 109
      0 177 149 110 159 inf 1 101 135  78 151  82 2  49 inf inf  63  65 3  45  47  38 inf  52 4  51  38  40  29  33

*/ 

#define ANSI_C 0     /* 0:   K&R function style convention */
#define VERIFY 0     /* 1:   produce the verification file */ 
#define FIRMACIND 0  /* 0,1: first machine index           */ 

#include <stdio.h>
#include <math.h>

struct problem {
  long rand_time;      /* random seed for jobs */ 
  short num_jobs;      /* number of jobs */ 
  short num_mach;      /* number of machines */ 
};

#if VERIFY == 1

struct problem S[] = {
  {         0,  0, 0},
  { 873654221, 20, 5},
  {         0,  0, 0}};

#else /* VERIFY */ 
    
struct problem S[] = {
{         0,     0,  0},
                         /* 20 jobs  5 machines */ 
{ 873654221,    20,  5},  
{ 379008056,    20,  5}, 
{ 1866992158,   20,  5}, 
{ 216771124,    20,  5}, 
{ 495070989,    20,  5}, 
{ 402959317,    20,  5}, 
{ 1369363414,   20,  5}, 
{ 2021925980,   20,  5},
{ 573109518,    20,  5}, 
{ 88325120,     20,  5}, 
                          /* 20 jobs  10 machines */ 
{ 587595453,    20, 10},
{ 1401007982,   20, 10},
{ 873136276,    20, 10}, 
{ 268827376,    20, 10}, 
{ 1634173168,   20, 10},
{ 691823909,    20, 10}, 
{ 73807235,     20, 10}, 
{ 1273398721,   20, 10}, 
{ 2065119309,   20, 10}, 
{ 1672900551,   20, 10},
                          /* 20 jobs 20 machines */
{ 479340445,    20, 20},  
{ 268827376,    20, 20},
{ 1958948863,   20, 20},
{ 918272953,    20, 20},
{ 555010963,    20, 20},
{ 2010851491,   20, 20},
{ 1519833303,   20, 20},
{ 1748670931,   20, 20},
{ 1923497586,   20, 20},
{ 1829909967,   20, 20},
                          /* 50 jobs  5 machines */  
{ 1328042058,   50,  5}, 
{ 200382020,    50,  5},
{ 496319842,    50,  5},
{ 1203030903,   50,  5},
{ 1730708564,   50,  5},
{ 450926852,    50,  5},
{ 1303135678,   50,  5},
{ 1273398721,   50,  5},
{ 587288402,    50,  5},
{ 248421594,    50,  5},
                          /* 50 Jobs 10 machines */ 
{ 1958948863,   50, 10},
{ 575633267,    50, 10},
{ 655816003,    50, 10}, 
{ 1977864101,   50, 10},
{ 93805469,     50, 10},
{ 1803345551,   50, 10},  
{ 49612559,     50, 10},
{ 1899802599,   50, 10},
{ 2013025619,   50, 10},
{ 578962478,    50, 10},
                          /* 50 jobs 20 machines */ 
{ 1539989115,   50, 20},
{ 691823909,    50, 20},
{ 655816003,    50, 20}, 
{ 1315102446,   50, 20}, 
{ 1949668355,   50, 20},
{ 1923497586,   50, 20},
{ 1805594913,   50, 20},
{ 1861070898,   50, 20}, 
{ 715643788,    50, 20}, 
{ 464843328,    50, 20}, 
                          /* 100 jobs  5 machines */ 
{ 896678084,   100,  5},
{ 1179439976,  100,  5}, 
{ 1122278347,  100,  5}, 
{ 416756875,   100,  5},
{ 267829958,   100,  5}, 
{ 1835213917,  100,  5}, 
{ 1328833962,  100,  5}, 
{ 1418570761,  100,  5}, 
{ 161033112,   100,  5},
{ 304212574,   100,  5}, 
                          /* 100 jobs 10 machines */ 
{ 1539989115,  100, 10},
{ 655816003,   100, 10}, 
{ 960914243,   100, 10}, 
{ 1915696806,  100, 10},
{ 2013025619,  100, 10}, 
{ 1168140026,  100, 10}, 
{ 1923497586,  100, 10}, 
{ 167698528,   100, 10}, 
{ 1528387973,  100, 10}, 
{ 993794175,   100, 10}, 
                          /* 100 jobs 20 machines */
{ 450926852,   100, 20},
{ 1462772409,  100, 20}, 
{ 1021685265,  100, 20}, 
{ 83696007,    100, 20}, 
{ 508154254,   100, 20}, 
{ 1861070898,  100, 20}, 
{ 26482542,    100, 20}, 
{ 444956424,   100, 20}, 
{ 2115448041,  100, 20}, 
{ 118254244,   100, 20}, 
                          /* 200 jobs 10 machines */ 
{ 471503978,   200, 10},
{ 1215892992,  200, 10}, 
{ 135346136,   200, 10}, 
{ 1602504050,  200, 10}, 
{ 160037322,   200, 10}, 
{ 551454346,   200, 10}, 
{ 519485142,   200, 10}, 
{ 383947510,   200, 10}, 
{ 1968171878,  200, 10}, 
{ 540872513,   200, 10}, 
                          /* 200 jobs 20 machines */
{ 2013025619,  200, 20},
{ 475051709,   200, 20}, 
{ 914834335,   200, 20}, 
{ 810642687,   200, 20},  
{ 1019331795,  200, 20}, 
{ 2056065863,  200, 20}, 
{ 1342855162,  200, 20}, 
{ 1325809384,  200, 20}, 
{ 1988803007,  200, 20}, 
{ 765656702,   200, 20}, 
                          /* 500 jobs 20 machines */
{ 1368624604,  500, 20},
{ 450181436,   500, 20}, 
{ 1927888393,  500, 20}, 
{ 1759567256,  500, 20}, 
{ 606425239,   500, 20}, 
{ 19268348,    500, 20}, 
{ 1298201670,  500, 20}, 
{ 2041736264,  500, 20},
{ 379756761,   500, 20},
{ 28837162,    500, 20},
{          0,    0,  0}};
#endif /* VERIFY */

/* generate a random number uniformly between low and high */

#if ANSI_C == 1
int unif (long *seed, short low, short high)
#else
int unif (seed, low, high)
long *seed; short low, high;
#endif
{
  static long m = 2147483647, a = 16807, b = 127773, c = 2836;
  double  value_0_1;              

  long k = *seed / b;
  *seed = a * (*seed % b) - k * c;
  if(*seed < 0) *seed = *seed + m;
  value_0_1 =  *seed / (double) m;

  return (int) (low + floor(value_0_1 * (high - low + 1)));
}

/* Maximal 500 jobs and 20 machines are provided. */
/* For larger problems extend array sizes.        */ 

int d[21][501];                       /* duration */ 
int dw[21][501][21];                    /* perturbed duration for workers */ 
int wm[2][20*20+1];                   /* workers machine assignation */
int ff[21][21];

#if ANSI_C == 1
void generate_flow_shop(short p)          /* Fill d and M according to S[p] */ 
#else
void generate_flow_shop(p)
short p;
#endif
{
  short i, j;
  long time_seed = S[p].rand_time;

  for(i = 0; i < S[p].num_mach; ++i)      /* determine a random duration */ 
    for (j = 0; j < S[p].num_jobs; ++j)   /* for all operations */ 
      d[i][j] = unif(&time_seed, 1, 99);  /* 99 = max. duration of op. */
}

#if ANSI_C == 1
void perturbate_flow_shop(short p)          /* Perturbate d for workers according to S[p] */ 
#else
void perturbate_flow_shop(p)
short p;
#endif
{
  short i, j, w, m, t;
  long time_seed = S[p].rand_time;

  for(i = 0; i < S[p].num_mach; ++i)      /* perturbate duration randomly */
    for (j = 0; j < S[p].num_jobs; ++j)   /* for all operations */
      for (w = 0; w < S[p].num_mach; ++w)   /* for all workers */
        dw[i][j][w] = unif(&time_seed, d[i][j], 5*d[i][j]);  /* perturbates between d and n*d */

  /* ensure one worker and one machine */
  for(i = 0; i < S[p].num_mach; ++i) ff[0][i]=i;
  for(i = 0; i < S[p].num_mach; ++i)
  {
    j = unif(&time_seed, 0, S[p].num_mach-1);
    w = ff[0][j];
    ff[0][j] = ff[0][i];
    ff[0][i] = w;
  }

  /* create incompatibilities */
  j=S[p].num_mach*(S[p].num_mach-1);
  for (i = 0; i < j; ++i)
  {
    wm[0][i] = (ff[0][i/S[p].num_mach] == i%S[p].num_mach)?
      S[p].num_mach-1: i/S[p].num_mach;
    wm[1][i] = i%S[p].num_mach;
  }

  for (i = 0; i < j; ++i)
  {
    t = unif(&time_seed, 0, j-1);
    w = wm[0][t];
    m = wm[1][t];
    wm[0][t] = wm[0][i];
    wm[1][t] = wm[1][i];
    wm[0][i] = w;
    wm[1][i] = m;
  }
}

#if ANSI_C == 1
void write_problem(short p, short incompat)  /* write out problem */ 
#else
void write_problem(p, incompat)
short p, incompat;
#endif
{
  short i, j, k, w;
  int z;
  FILE *f = NULL;
  char name[11];

  sprintf(name,"ta%03dI%d.txt", p, incompat); /* file name construction */ 
  if(!(f = fopen(name,"w"))) {               /* open file for writing  */ 
    fprintf(stderr,"file %s error\n", name);
    return;
  }
  fprintf(f,"%d %d %d\n", S[p].num_jobs, S[p].num_mach, S[p].num_mach); /* write header line */ 

  z = S[p].num_mach;
  z = z * z / 10;
  z = z * incompat;
  for(i = 0; i < S[p].num_mach; ++i)
    for(w = 0; w < S[p].num_mach; ++w)
      ff[i][w] = 0;
  for(i = 0; i < z; ++i)
    ff[wm[0][i]][wm[1][i]] = 1;

  for(j = 0; j < S[p].num_jobs; ++j) {
    for(i = 0; i < S[p].num_mach; ++i) {
      fprintf(f," %d", i+FIRMACIND);   /* write machine */ 
      for(w = 0; w < S[p].num_mach; ++w) {
        if(ff[i][w]/* != 0 */){
          fprintf(f," inf");}
        else{
          fprintf(f," %3d", dw[i][j][w]);}
      }
    }
    fprintf(f,"\n");                         /* newline == End of job */ 
  }
  fclose(f);                                 /* close file */ 
}


int main()                                    
{
  short i = 1;
  while(S[i].rand_time) {                    /* for i == 1 up to NULL entry */
    generate_flow_shop(i);                   /* generate problem i  */ 
    perturbate_flow_shop(i);                 /* generate problem i  */ 
    write_problem(i,0);                      /* write out problem i */ 
    write_problem(i,1);                      /* write out problem i */ 
    write_problem(i,2);                      /* write out problem i */ 
    ++i;                                     /* increment i */ 
  }
  return 0;
}

