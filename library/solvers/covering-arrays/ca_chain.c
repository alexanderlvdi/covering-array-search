/* ca_chain.c — guided ejection-chain / min-conflicts solver for strength-3 covering arrays
 * over V symbols: the move-selection recipe that cracked CA(89;3,42,3) (attempt 7 of
 * problems/ca-3-3-k42-n89, chain.py), rewritten in C so a multi-seed sweep is cheap.
 *
 * State: N x K array over {0..V-1}.  CNT[triple,code] = number of rows carrying that tuple.
 * Move: flip one cell (r,c) to another symbol.  Exact dU = #(lost entries with CNT==1)
 *       - #(gained entries with CNT==0)   (disjoint, because old and new symbols differ).
 * Search: attack ONE uncovered tuple; enumerate only its Hamming-1 closers (rows matching the
 *       wanted tuple in 2 of its 3 columns); pick min dU, ties by max `protection` (singletons
 *       turned into doubles: a pay-forward ejection chain).  Tabu on (row,col).  Hamming-2 rows
 *       are the fallback when a tuple has no one-cell closer.  Two random flips (a kick) when
 *       the walk is flat at its best value.
 *
 * usage: ca_chain N K V in out seed secs [tenure] [kick]
 *        prints: uncovered_tuples=<best>
 * search-kind: aimed
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int N, K, V, VS, NT, NCOLE;
static unsigned char (*A)[64];
static unsigned char *KEY;          /* [N][NT] tuple code of row r on triple ti */
static short *CNT;                  /* [NT][VS] */
static int *colstart, *colti, *colmul;
static int (*tri)[3];
static int *unc; static int U, UNCCAP;
static unsigned char *isunc;
static int tabu[256][64];
static unsigned long rs = 88172645463325252UL;
static inline unsigned long rnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static void die(const char *m){ fprintf(stderr, "%s\n", m); exit(1); }
static int nC3(int k){ return k*(k-1)*(k-2)/6; }

static void unc_add(int idx){
  if(isunc[idx]) return;
  isunc[idx]=1;
  if(U>=UNCCAP){ UNCCAP*=2; unc=realloc(unc, (size_t)UNCCAP*sizeof(int)); }
  unc[U++]=idx;
}
static void unc_del(int idx){
  if(!isunc[idx]) return;
  isunc[idx]=0;
  for(int i=0;i<U;i++) if(unc[i]==idx){ unc[i]=unc[--U]; return; }
}

/* exact delta of flipping (r,c) to v; prot = #new tuples whose count is 1 (become 2) */
static inline void evalc(int r,int c,int v,int *dUo,int *proto){
  int d = v - A[r][c], lost=0, gain=0, prot=0;
  for(int j=colstart[c]; j<colstart[c+1]; j++){
    int ti=colti[j], base=ti*VS, ok=KEY[r*NT+ti], nk=ok+d*colmul[j];
    if(CNT[base+ok]==1) lost++;
    if(CNT[base+nk]==0) gain++;
    if(CNT[base+nk]==1) prot++;
  }
  *dUo = lost-gain; *proto = prot;
}

static void apply(int r,int c,int v){
  int d = v - A[r][c];
  if(!d) return;
  for(int j=colstart[c]; j<colstart[c+1]; j++){
    int ti=colti[j], base=ti*VS, mul=colmul[j];
    int ok=KEY[r*NT+ti], nk=ok+d*mul;
    if(CNT[base+ok]==1) unc_add(base+ok);
    CNT[base+ok]--;
    if(CNT[base+nk]==0) unc_del(base+nk);
    CNT[base+nk]++;
    KEY[r*NT+ti]=(unsigned char)nk;
  }
  A[r][c]=(unsigned char)v;
}

static void write_arr(const char *p, unsigned char (*M)[64]){
  FILE *f=fopen(p,"w"); if(!f) die("cannot write out");
  for(int r=0;r<N;r++){ for(int j=0;j<K;j++) fprintf(f, j==K-1?"%d\n":"%d ", M[r][j]); }
  fclose(f);
}

/* closers of uncovered (ti,code): Hamming-1 rows, else up to 12 Hamming-2 rows */
static int closing(int ti,int code,int *cr,int *cc,int *cv,int *dist){
  int cols[3]={tri[ti][0],tri[ti][1],tri[ti][2]};
  int want[3]={code/(V*V),(code/V)%V,code%V};
  int n=0, nh2=0;
  *dist=1;
  for(int r=0;r<N;r++){
    int d=(A[r][cols[0]]!=want[0])+(A[r][cols[1]]!=want[1])+(A[r][cols[2]]!=want[2]);
    if(d==1){
      for(int p=0;p<3;p++) if(A[r][cols[p]]!=want[p]){ cr[n]=r; cc[n]=cols[p]; cv[n]=want[p]; n++; break; }
    }
  }
  if(n) return n;
  *dist=2;
  for(int r=0;r<N && nh2<12;r++){
    int d=(A[r][cols[0]]!=want[0])+(A[r][cols[1]]!=want[1])+(A[r][cols[2]]!=want[2]);
    if(d==2){
      nh2++;
      for(int p=0;p<3;p++) if(A[r][cols[p]]!=want[p]){ cr[n]=r; cc[n]=cols[p]; cv[n]=want[p]; n++; }
    }
  }
  return n;
}

int main(int argc,char**argv){
  if(argc<8){ fprintf(stderr,"usage: ca_chain N K V in out seed secs [tenure] [kick]\n"); return 1; }
  N=atoi(argv[1]); K=atoi(argv[2]); V=atoi(argv[3]);
  const char *inp=argv[4], *outp=argv[5];
  rs ^= (unsigned long)strtoul(argv[6],0,10)*2654435761UL + 12345UL;
  int secs=atoi(argv[7]);
  int tenure = argc>8?atoi(argv[8]):8;
  int kick   = argc>9?atoi(argv[9]):500;
  if(K>64||N>256) die("N/K too big");
  VS=V*V*V; NT=nC3(K); NCOLE=nC3(K-1);
  A=calloc(256,64); if(!A) die("oom A");
  KEY=malloc((size_t)N*NT); if(!KEY) die("oom key");
  CNT=calloc((size_t)NT*VS,sizeof(short)); if(!CNT) die("oom cnt");
  isunc=calloc((size_t)NT*VS,1); unc=malloc(64*sizeof(int)); UNCCAP=64;
  tri=malloc((size_t)NT*3*sizeof(int)); if(!tri) die("oom tri");
  int t=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++){ tri[t][0]=a;tri[t][1]=b;tri[t][2]=c;t++; }
  colstart=malloc((size_t)(K+1)*sizeof(int));
  colti=malloc((size_t)K*NCOLE*sizeof(int)); colmul=malloc((size_t)K*NCOLE*sizeof(int));
  int *nl=calloc(K,sizeof(int));
  int *cti=malloc((size_t)K*NCOLE*sizeof(int)); int *cpos=malloc((size_t)K*NCOLE*sizeof(int));
  for(int ti2=0;ti2<NT;ti2++) for(int p=0;p<3;p++){
    int c=tri[ti2][p]; int q=nl[c]++;
    cti[c*NCOLE+q]=ti2; cpos[c*NCOLE+q]=p;
  }
  for(int c=0;c<K;c++) colstart[c]=c*NCOLE;
  colstart[K]=K*NCOLE;
  for(int c=0;c<K;c++) for(int j=0;j<nl[c];j++){
    int ti2=cti[c*NCOLE+j], p=cpos[c*NCOLE+j];
    colti[c*NCOLE+j]=ti2;
    colmul[c*NCOLE+j] = (p==0? V*V : (p==1? V : 1));
  }
  FILE *f=fopen(inp,"r"); if(!f) die("cannot open in");
  for(int r=0;r<N;r++) for(int j=0;j<K;j++){ int v; if(fscanf(f,"%d",&v)!=1) die("bad input"); A[r][j]=(unsigned char)v; }
  fclose(f);
  for(int ti2=0;ti2<NT;ti2++){
    int a=tri[ti2][0], b=tri[ti2][1], c=tri[ti2][2];
    for(int r=0;r<N;r++){
      int code=A[r][a]*V*V + A[r][b]*V + A[r][c];
      KEY[r*NT+ti2]=(unsigned char)code;
      CNT[ti2*VS+code]++;
    }
  }
  U=0;
  for(int i=0;i<NT*VS;i++) if(CNT[i]==0){ isunc[i]=1; if(U>=UNCCAP){UNCCAP*=2;unc=realloc(unc,(size_t)UNCCAP*sizeof(int));} unc[U++]=i; }
  memset(tabu,0xff,sizeof(tabu));
  int U0=U, best=U;
  unsigned char (*bestA)[64]=malloc(256*64);
  memcpy(bestA,A,256*64);
  long long steps=0;
  clock_t t0=clock();
  printf("ca_chain N=%d k=%d V=%d U0=%d seed=%s tenure=%d kick=%d\n",N,K,V,U0,argv[6],tenure,kick);
  fflush(stdout);
  static int cr[1024], cc[1024], cv[1024], cdU[1024], cprot[1024];
  static int hcr[1024], hcc[1024], hcv[1024];
  for(;;){
    steps++;
    if(U==0) break;
    if((steps&255)==0 && ((double)(clock()-t0)/CLOCKS_PER_SEC) > secs) break;
    int chosen=0, chosen_ti=0, chosen_code=0, nc=0;
    int starti = (int)(rnd()%(unsigned)U);
    for(int tries=0; tries<U && tries<6; tries++){
      int ui=(starti+tries)%U;
      int idx=unc[ui], ti2=idx/VS, code=idx%VS, ndist;
      int ncl=closing(ti2,code,hcr,hcc,hcv,&ndist);
      int m=0;
      for(int q=0;q<ncl;q++){
        if(tabu[hcr[q]][hcc[q]] > (int)steps) continue;
        cr[m]=hcr[q]; cc[m]=hcc[q]; cv[m]=hcv[q]; m++;
      }
      if(m>0){ chosen=1; chosen_ti=ti2; chosen_code=code; nc=m; break; }
    }
    if(!chosen){ memset(tabu,0xff,sizeof(tabu)); continue; }
    for(int q=0;q<nc;q++) evalc(cr[q],cc[q],cv[q],&cdU[q],&cprot[q]);
    int bi=0;
    for(int q=1;q<nc;q++) if(cdU[q]<cdU[bi] || (cdU[q]==cdU[bi] && cprot[q]>cprot[bi])) bi=q;
    int r=cr[bi], c=cc[bi], nv=cv[bi];
    apply(r,c,nv);
    tabu[r][c]=(int)steps+tenure;
    if(U<best){
      best=U; memcpy(bestA,A,256*64);
      printf("[%7.2fs] step %lld U=%d attack=(%d,%d,%d)->%d move=(%d,%d)->%d dU=%+d prot=%d\n",
             (double)(clock()-t0)/CLOCKS_PER_SEC, steps, U, tri[chosen_ti][0],tri[chosen_ti][1],tri[chosen_ti][2],
             chosen_code, r,c,nv, cdU[bi], cprot[bi]);
      fflush(stdout);
      if(best==0) break;
    }
    if(kick && (steps % kick)==0 && best==U){
      for(int z=0;z<2;z++){
        int rr=(int)(rnd()%(unsigned)N), ccc=(int)(rnd()%(unsigned)K), nv2=(int)(rnd()%(unsigned)V);
        if(nv2!=A[rr][ccc]) apply(rr,ccc,nv2);
      }
      memset(tabu,0xff,sizeof(tabu));
    }
  }
  write_arr(outp,bestA);
  printf("uncovered_tuples=%d steps=%lld U0=%d best=%d %.1fs\n",best,steps,U0,best,(double)(clock()-t0)/CLOCKS_PER_SEC);
  return 0;
}
