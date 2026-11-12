// Exact list sizes on the native circle code N = {deg P <= 16, P0 + P16 = 0} on mu32 in F_(p^2), p = 2^31-1,
// agreement 20. A 20-set S (complement C, 12 points) passes for word f = sum a_m x^m iff f is orthogonal to the 4 dual
// vectors of N vanishing on C: g_e = x^(e+1) Q_C (e=0,1,2) and g4 = 1 + x^16 + x R(x), R = interpolant of
// -(1+t^16)/t on C (deg <= 11). Pairing <f,g> = 32 * sum_m a_m g_{(32-m) mod 32}. Flag "rs" drops g4 (RS[mu32,17] check).
// Usage: ./native THREADS words.txt [rs]   ; words.txt lines: m re im m re im ...
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef uint64_t u64;
#define P 2147483647ULL
typedef struct { u64 r, i; } fe;
static inline u64 red(u64 x) { x = (x & P) + (x >> 31); x = (x & P) + (x >> 31); return x >= P ? x - P : x; }
static inline fe fmul(fe a, fe b) { u64 ac = red(a.r * b.r), bd = red(a.i * b.i); fe z = { red(ac + P - bd), red(a.r * b.i + a.i * b.r) }; return z; }
static inline fe fsub(fe a, fe b) { fe z = { red(a.r + P - b.r), red(a.i + P - b.i) }; return z; }
static inline fe fadd(fe a, fe b) { fe z = { red(a.r + b.r), red(a.i + b.i) }; return z; }
static inline int fz(fe a) { return a.r == 0 && a.i == 0; }
static u64 pw(u64 b, u64 e) { u64 r = 1; while (e) { if (e & 1) r = red(r * b); b = red(b * b); e >>= 1; } return r; }
static fe finv(fe a) { u64 n = red(red(a.r * a.r) + red(a.i * a.i)); u64 ni = pw(n, P - 2); fe z = { red(a.r * ni), red((P - a.i) % P * ni) }; return z; }
static fe Z[32], IZ[32], ID[32][32];
#define MAXW 4096
static int nw, wn[MAXW], wm[MAXW][16]; static fe wa[MAXW][16]; static int rsonly;
static u64 g_pass[MAXW], g_dist[MAXW], g_maxag[MAXW], g_hist[MAXW][33];
static int tj0[20000], tj1[20000], ntask, nexttask; static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;

// dual vectors for complement set given by list cs[12]; out g[4][32]
static void duals(const int *cs, fe g[4][32]) {
    fe q[13]; memset(q, 0, sizeof q); q[0].r = 1;
    for (int d = 0; d < 12; d++) { for (int t = d + 1; t >= 0; t--) { fe hi = t >= 1 ? q[t - 1] : (fe){0,0}, lo = t <= d ? fmul(q[t], Z[cs[d]]) : (fe){0,0}; q[t] = fsub(hi, lo); } }
    memset(g, 0, sizeof(fe) * 4 * 32);
    for (int e = 0; e < 3; e++) for (int j = 0; j <= 12; j++) g[e][e + 1 + j] = q[j];
    if (rsonly) return;
    fe R[12]; memset(R, 0, sizeof R);
    for (int a = 0; a < 12; a++) {
        int t = cs[a]; if (t & 1) continue;               // t^16 = (-1)^index; value 0 at odd points
        fe den = { 1, 0 }; for (int b = 0; b < 12; b++) if (b != a) den = fmul(den, ID[t][cs[b]]); // 1/prod(t-u)
        fe v = fmul(den, IZ[t]); v.r = (2 * (P - v.r)) % P; v.i = (2 * (P - v.i)) % P;   // -2/t / prod
        // Q_C/(x - t): synthetic division, degree 11
        fe qq[12]; qq[11] = q[12]; for (int j = 11; j >= 1; j--) qq[j - 1] = fadd(q[j], fmul(Z[t], qq[j]));
        for (int j = 0; j < 12; j++) R[j] = fadd(R[j], fmul(v, qq[j]));
    }
    g[3][0].r = 1; g[3][16].r = 1; for (int j = 0; j < 12; j++) g[3][j + 1] = fadd(g[3][j + 1], R[j]);
}

// ---- four-monomial census on the native code: supports m1<..<m4 in 16..31 (1820). Mode 0 = probe (counts), 1 = full dump.
static int S4[1820][4], ns4, MODE, TLIM; static FILE *outp, *outl;
static u64 g_r3, g_r2, g_r1, g_r0, g_pt, g_ptc, s_low[1820];
static int pidx[16][16];
#define SK (1u<<27)
static uint16_t *sk0,*sk1; static int THR=20;
static inline u64 mix(u64 x){x^=x>>33;x*=0xff51afd7ed558ccdULL;x^=x>>33;x*=0xc4ceb9fe1a85ec53ULL;x^=x>>33;return x;}
static int tj2[20000];
static inline fe det3(fe a0,fe a1,fe a2,fe b0,fe b1,fe b2,fe c0,fe c1,fe c2){
  return fadd(fadd(fmul(a0,fsub(fmul(b1,c2),fmul(b2,c1))),fmul(a1,fsub(fmul(b2,c0),fmul(b0,c2)))),fmul(a2,fsub(fmul(b0,c1),fmul(b1,c0))));}
static int passw(const int *mm,const fe *aa,uint32_t cmask){ int cs[12],n=0; for(int j=0;j<32;j++) if((cmask>>j)&1) cs[n++]=j; fe g[4][32]; duals(cs,g);
  for(int e=0;e<4;e++){ fe s={0,0}; for(int k=0;k<4;k++) s=fadd(s,fmul(aa[k],g[e][(32-mm[k])&31])); if(!fz(s)) return 0;} return 1;}
static void leaf(const int *cs,u64 *loc,u64 *low){
  fe g[4][32]; duals(cs,g); fe G[4][16]; for(int e=0;e<4;e++) for(int m=0;m<16;m++) G[e][m]=g[e][(16-m)&31]; // column m <-> monomial 16+m
  fe A[120],B[120]; int k=0;
  for(int i=0;i<16;i++) for(int j=i+1;j<16;j++){ A[k]=fsub(fmul(G[0][i],G[1][j]),fmul(G[0][j],G[1][i])); B[k]=fsub(fmul(G[2][i],G[3][j]),fmul(G[2][j],G[3][i])); k++; }
  uint32_t cm=0; for(int a=0;a<12;a++) cm|=1u<<cs[a]; int maxS=31; while((cm>>maxS)&1) maxS--;
  for(int s=0;s<ns4;s++){ int a=S4[s][0]-16,b=S4[s][1]-16,c=S4[s][2]-16,d=S4[s][3]-16;
    // Laplace on rows (0,1)|(2,3): det = A_ab B_cd - A_ac B_bd + A_ad B_bc + A_bc B_ad - A_bd B_ac + A_cd B_ab
    fe t=fsub(fmul(A[pidx[a][b]],B[pidx[c][d]]),fmul(A[pidx[a][c]],B[pidx[b][d]]));
    t=fadd(t,fmul(A[pidx[a][d]],B[pidx[b][c]])); t=fadd(t,fmul(A[pidx[b][c]],B[pidx[a][d]]));
    t=fsub(t,fmul(A[pidx[b][d]],B[pidx[a][c]])); t=fadd(t,fmul(A[pidx[c][d]],B[pidx[a][b]]));
    if(!fz(t)) continue;
    int col[4]={a,b,c,d}; fe kv[4]; int got=0;
    for(int om=0;om<4&&!got;om++){ int r[3],n=0; for(int i=0;i<4;i++) if(i!=om) r[n++]=i;
      for(int q=0;q<4;q++){ int cc[3],m=0; for(int i=0;i<4;i++) if(i!=q) cc[m++]=col[i];
        fe v=det3(G[r[0]][cc[0]],G[r[0]][cc[1]],G[r[0]][cc[2]],G[r[1]][cc[0]],G[r[1]][cc[1]],G[r[1]][cc[2]],G[r[2]][cc[0]],G[r[2]][cc[1]],G[r[2]][cc[2]]);
        if(q&1){ v.r=(P-v.r)%P; v.i=(P-v.i)%P; } kv[q]=v; }
      if(!(fz(kv[0])&&fz(kv[1])&&fz(kv[2])&&fz(kv[3]))) got=1; }
    if(!got){ loc[1]++; low[s]++;
      if(MODE==1){ pthread_mutex_lock(&mu); fprintf(outl,"%d %u\n",s,cm); pthread_mutex_unlock(&mu); }
      continue; }
    loc[0]++;
    if(fz(kv[0])||fz(kv[1])||fz(kv[2])||fz(kv[3])) continue;
    loc[2]++;
    if(MODE==0) continue;
    fe iv=finv(kv[0]); fe aa[4]={{1,0},fmul(kv[1],iv),fmul(kv[2],iv),fmul(kv[3],iv)}; int mm[4]={a+16,b+16,c+16,d+16};
    int canon=1,ag=20;
    for(int t2=0;t2<12&&canon;t2++){ int y=cs[t2]; if(y>maxS) break; uint32_t c2=(cm&~(1u<<y))|(1u<<maxS); if(passw(mm,aa,c2)) canon=0; }
    if(canon){ loc[3]++; u64 h=mix((u64)s*0x9e3779b97f4a7c15ULL^aa[1].r); h=mix(h^aa[1].i); h=mix(h^aa[2].r); h=mix(h^aa[2].i); h=mix(h^aa[3].r); h=mix(h^aa[3].i);
      uint32_t i0=h&(SK-1), i1=(h>>32)&(SK-1); uint16_t v0=__atomic_add_fetch(&sk0[i0],1,__ATOMIC_RELAXED), v1=__atomic_add_fetch(&sk1[i1],1,__ATOMIC_RELAXED);
      uint16_t v=v0<v1?v0:v1; if(v!=THR) continue; pthread_mutex_lock(&mu); fprintf(outp,"%d %llu %llu %llu %llu %llu %llu %d\n",s,(unsigned long long)aa[1].r,(unsigned long long)aa[1].i,(unsigned long long)aa[2].r,(unsigned long long)aa[2].i,(unsigned long long)aa[3].r,(unsigned long long)aa[3].i,ag); pthread_mutex_unlock(&mu); }
  }
}
static void rec(int *cs,int d,int s,u64 *loc,u64 *low){ if(d==12){leaf(cs,loc,low);return;} for(int j=s;j<=32-(12-d);j++){cs[d]=j;rec(cs,d+1,j+1,loc,low);} }
static void *worker(void *x){ (void)x; u64 loc[4]={0},*low=calloc(1820,8); int cs[12];
  for(;;){ pthread_mutex_lock(&mu); int k=nexttask++; pthread_mutex_unlock(&mu); if(k>=ntask) break; cs[0]=tj0[k]; cs[1]=tj1[k]; cs[2]=tj2[k]; rec(cs,3,tj2[k]+1,loc,low); }
  pthread_mutex_lock(&mu); g_r3+=loc[0]; g_r2+=loc[1]; g_pt+=loc[2]; g_ptc+=loc[3]; for(int s=0;s<1820;s++) s_low[s]+=low[s]; pthread_mutex_unlock(&mu); return NULL; }
int main(int argc,char **argv){ int nth=atoi(argv[1]); MODE=atoi(argv[2]); TLIM=atoi(argv[3]);
  if(MODE==1){ outp=fopen(argv[4],"w"); outl=fopen(argv[5],"w"); }
  fe w={906276279ULL,1179735656ULL}; Z[0].r=1; Z[0].i=0; for(int j=1;j<32;j++) Z[j]=fmul(Z[j-1],w);
  for(int j=0;j<32;j++){ IZ[j]=finv(Z[j]); for(int k=0;k<32;k++) if(k!=j) ID[j][k]=finv(fsub(Z[j],Z[k])); }
  int k=0; for(int i=0;i<16;i++) for(int j=i+1;j<16;j++){ pidx[i][j]=k; pidx[j][i]=k; k++; }
  for(int a=16;a<32;a++) for(int b=a+1;b<32;b++) for(int c=b+1;c<32;c++) for(int d=c+1;d<32;d++){S4[ns4][0]=a;S4[ns4][1]=b;S4[ns4][2]=c;S4[ns4][3]=d;ns4++;}
  sk0=calloc(SK,2); sk1=calloc(SK,2);
  for(int a=0;a<32;a++) for(int b=a+1;b<32;b++) for(int c=b+1;c<32;c++) if(32-c-1>=9 && ntask<TLIM){tj0[ntask]=a;tj1[ntask]=b;tj2[ntask]=c;ntask++;}
  pthread_t th[64]; for(int q=0;q<nth;q++) pthread_create(&th[q],0,worker,0); for(int q=0;q<nth;q++) pthread_join(th[q],0);
  if(MODE==1){ fclose(outp); fclose(outl); }
  printf("tasks=%d rank3(point) sets=%llu  full-support points=%llu  canonical=%llu  rank<=2 sets=%llu\n",ntask,(unsigned long long)g_r3,(unsigned long long)g_pt,(unsigned long long)g_ptc,(unsigned long long)g_r2);
  int shown=0; for(int s=0;s<ns4;s++) if(s_low[s]) shown++; printf("supports with rank<=2 sets: %d\n",shown);
  return 0; }
