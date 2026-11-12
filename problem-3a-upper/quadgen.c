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
static int tj0[600], tj1[600], ntask, nexttask; static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;

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

// Reads quad_low.txt (support index, complement mask); for each set computes the solution space of the 4x4 system on the
// support in reduced row-echelon form of its row space (canonical key). Prints: support dim(kernel) rref-rows.
static int S4[1820][4];
static u64 mx(u64 x){x^=x>>33;x*=0xff51afd7ed558ccdULL;x^=x>>33;x*=0xc4ceb9fe1a85ec53ULL;x^=x>>33;return x;}
static int passw4(const int *mm,const fe *aa,uint32_t cmask){ int cs[12],n=0; for(int j=0;j<32;j++) if((cmask>>j)&1) cs[n++]=j; fe g[4][32]; duals(cs,g);
  for(int e=0;e<4;e++){ fe s={0,0}; for(int k=0;k<4;k++) s=fadd(s,fmul(aa[k],g[e][(32-mm[k])&31])); if(!fz(s)) return 0;} return 1;}
int main(int argc,char **argv){
  fe w={906276279ULL,1179735656ULL}; Z[0].r=1; Z[0].i=0; for(int j=1;j<32;j++) Z[j]=fmul(Z[j-1],w);
  for(int j=0;j<32;j++){ IZ[j]=finv(Z[j]); for(int k=0;k<32;k++) if(k!=j) ID[j][k]=finv(fsub(Z[j],Z[k])); }
  int n=0; for(int a=16;a<32;a++) for(int b=a+1;b<32;b++) for(int c=b+1;c<32;c++) for(int d=c+1;d<32;d++){S4[n][0]=a;S4[n][1]=b;S4[n][2]=c;S4[n][3]=d;n++;}
  int s; unsigned cm; FILE *f=fopen(argv[1],"r");
  while(fscanf(f,"%d %u",&s,&cm)==2){
    int cs[12],k=0; for(int j=0;j<32;j++) if((cm>>j)&1) cs[k++]=j; fe g[4][32]; duals(cs,g);
    fe M[4][4]; for(int e=0;e<4;e++) for(int q=0;q<4;q++) M[e][q]=g[e][(32-S4[s][q])&31];
    int r=0, piv[4];
    for(int c=0;c<4&&r<4;c++){ int p=-1; for(int i=r;i<4;i++) if(!fz(M[i][c])){p=i;break;} if(p<0) continue;
      for(int q=0;q<4;q++){ fe t=M[r][q]; M[r][q]=M[p][q]; M[p][q]=t; }
      fe iv=finv(M[r][c]); for(int q=0;q<4;q++) M[r][q]=fmul(M[r][q],iv);
      for(int i=0;i<4;i++) if(i!=r&&!fz(M[i][c])){ fe t=M[i][c]; for(int q=0;q<4;q++) M[i][q]=fsub(M[i][q],fmul(t,M[r][q])); }
      piv[r++]=c; }
    // generic word on the kernel, seeded by the rref key so every set of the same space uses the same word
    u64 h=mx(s+1); for(int i=0;i<r;i++) for(int q=0;q<4;q++) h=mx(h^M[i][q].r^(M[i][q].i<<31));
    int isp[4]={0}; for(int i=0;i<r;i++) isp[piv[i]]=1; fe wv[4]={{0,0},{0,0},{0,0},{0,0}};
    for(int f=0;f<4;f++) if(!isp[f]){ h=mx(h+f); fe c={h%P,(h>>32)%P}; wv[f]=fadd(wv[f],c); for(int i=0;i<r;i++){ fe t=fmul(c,M[i][f]); wv[piv[i]]=fsub(wv[piv[i]],t);} }
    int canon=1; int maxS=31; while((cm>>maxS)&1) maxS--;
    for(int t=0;t<12&&canon;t++){ int y=cs[t]; if(y>maxS) break; if(passw4(S4[s],wv,(cm&~(1u<<y))|(1u<<maxS))) canon=0; }
    printf("%d %d %d",canon,s,4-r); for(int i=0;i<r;i++) for(int q=0;q<4;q++) printf(" %llu %llu",(unsigned long long)M[i][q].r,(unsigned long long)M[i][q].i); printf("\n");
  }
  return 0; }
