#define LINESONLY 1
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
static int nw, wn[MAXW], wm[MAXW][16]; static fe wa[MAXW][16]; static int rsonly=0;
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

// ---- three-monomial census on the native code: supports m1<m2<m3 in 16..31 (560 supports) ----
static int S3[560][3], ns3;
static FILE *outp; static u64 g_rk[4], g_line[560], g_all[560];
static inline fe d3(fe a0,fe a1,fe a2,fe b0,fe b1,fe b2,fe c0,fe c1,fe c2){
  return fadd(fadd(fmul(a0,fsub(fmul(b1,c2),fmul(b2,c1))),fmul(a1,fsub(fmul(b2,c0),fmul(b0,c2)))),fmul(a2,fsub(fmul(b0,c1),fmul(b1,c0))));}
static int passw(const int *mm, const fe *aa, uint32_t cmask){ int cs[12],n=0; for(int j=0;j<32;j++) if((cmask>>j)&1) cs[n++]=j; fe g[4][32]; duals(cs,g);
  for(int e=0;e<4;e++){ fe s={0,0}; for(int k=0;k<3;k++) s=fadd(s,fmul(aa[k],g[e][(32-mm[k])&31])); if(!fz(s)) return 0;} return 1;}
static void leaf(const int *cs, u64 *rk, u64 *ln, u64 *al){
  fe g[4][32]; duals(cs,g); fe G[4][32]; for(int e=0;e<4;e++) for(int m=16;m<32;m++) G[e][m]=g[e][(32-m)&31];
  uint32_t cm=0; for(int a=0;a<12;a++) cm|=1u<<cs[a]; int maxS=31; while((cm>>maxS)&1) maxS--;
  for(int s=0;s<ns3;s++){ int a=S3[s][0],b=S3[s][1],c=S3[s][2]; int r[4]={0,1,2,3}; int bad=0;
    for(int q=0;q<4&&!bad;q++){ int i0=q==0?1:0, i1=q<=1?2:1, i2=q<=2?3:2; // omit row q
      if(!fz(d3(G[i0][a],G[i0][b],G[i0][c],G[i1][a],G[i1][b],G[i1][c],G[i2][a],G[i2][b],G[i2][c]))) bad=1; }
    (void)r; if(bad){ continue; }
    // rank <= 2: find two independent rows -> kernel by cross product
    fe k[3]; int found=0;
    for(int i=0;i<4&&!found;i++) for(int j=i+1;j<4&&!found;j++){
      k[0]=fsub(fmul(G[i][b],G[j][c]),fmul(G[i][c],G[j][b])); k[1]=fsub(fmul(G[i][c],G[j][a]),fmul(G[i][a],G[j][c])); k[2]=fsub(fmul(G[i][a],G[j][b]),fmul(G[i][b],G[j][a]));
      if(!(fz(k[0])&&fz(k[1])&&fz(k[2]))) found=1; }
    if(!found){ int nz=0; for(int i=0;i<4;i++) if(!(fz(G[i][a])&&fz(G[i][b])&&fz(G[i][c]))) nz=1; if(nz){ rk[1]++; ln[s]++; for(int i=0;i<4;i++){ fe r0=G[i][a],r1=G[i][b],r2=G[i][c]; if(fz(r0)&&fz(r1)&&fz(r2)) continue; fe nr[3]={r0,r1,r2}; int p0=0; while(fz(nr[p0])) p0++; fe iv=finv(nr[p0]); for(int t=0;t<3;t++) nr[t]=fmul(nr[t],iv); pthread_mutex_lock(&mu); fprintf(outp,"L %d %d %d %llu %llu %llu %llu %llu %llu\n",a,b,c,(unsigned long long)nr[0].r,(unsigned long long)nr[0].i,(unsigned long long)nr[1].r,(unsigned long long)nr[1].i,(unsigned long long)nr[2].r,(unsigned long long)nr[2].i); pthread_mutex_unlock(&mu); break; } } else { rk[0]++; al[s]++; } continue; }
    rk[2]++; if(LINESONLY) continue;
    if(fz(k[0])||fz(k[1])||fz(k[2])) continue; // fewer monomials: covered by M-0908/M-0911
    fe iv=finv(k[0]); fe aa[3]={{1,0},fmul(k[1],iv),fmul(k[2],iv)}; int mm[3]={a,b,c};
    int canon=1, ag=20;
    for(int t=0;t<12;t++){ int y=cs[t]; uint32_t c2=(cm&~(1u<<y))|(1u<<maxS); if(passw(mm,aa,c2)){ ag++; if(y<maxS) canon=0; } }
    if(canon){ pthread_mutex_lock(&mu); fprintf(outp,"%d %d %d %llu %llu %llu %llu %d\n",a,b,c,(unsigned long long)aa[1].r,(unsigned long long)aa[1].i,(unsigned long long)aa[2].r,(unsigned long long)aa[2].i,ag); pthread_mutex_unlock(&mu); }
  }
}
static void rec(int *cs,int d,int s,u64 *rk,u64 *ln,u64 *al){ if(d==12){leaf(cs,rk,ln,al);return;} for(int j=s;j<=32-(12-d);j++){cs[d]=j;rec(cs,d+1,j+1,rk,ln,al);} }
static void *worker(void *x){ (void)x; u64 rk[4]={0},*ln=calloc(560,8),*al=calloc(560,8); int cs[12];
  for(;;){ pthread_mutex_lock(&mu); int k=nexttask++; pthread_mutex_unlock(&mu); if(k>=ntask) break; cs[0]=tj0[k]; cs[1]=tj1[k]; rec(cs,2,tj1[k]+1,rk,ln,al); }
  pthread_mutex_lock(&mu); for(int q=0;q<4;q++) g_rk[q]+=rk[q]; for(int s=0;s<560;s++){g_line[s]+=ln[s]; g_all[s]+=al[s];} pthread_mutex_unlock(&mu); return NULL; }
int main(int argc,char **argv){ int nth=atoi(argv[1]); outp=fopen(argv[2],"w");
  fe w={906276279ULL,1179735656ULL}; Z[0].r=1; Z[0].i=0; for(int j=1;j<32;j++) Z[j]=fmul(Z[j-1],w);
  for(int j=0;j<32;j++){ IZ[j]=finv(Z[j]); for(int k=0;k<32;k++) if(k!=j) ID[j][k]=finv(fsub(Z[j],Z[k])); }
  for(int a=16;a<32;a++) for(int b=a+1;b<32;b++) for(int c=b+1;c<32;c++){S3[ns3][0]=a;S3[ns3][1]=b;S3[ns3][2]=c;ns3++;}
  for(int a=0;a<32;a++) for(int b=a+1;b<32;b++) if(32-b-1>=10){tj0[ntask]=a;tj1[ntask]=b;ntask++;}
  pthread_t th[64]; for(int k=0;k<nth;k++) pthread_create(&th[k],0,worker,0); for(int k=0;k<nth;k++) pthread_join(th[k],0); fclose(outp);
  printf("rank2(point) sets=%llu  rank1(line) sets=%llu  rank0(all-word) sets=%llu\n",(unsigned long long)g_rk[2],(unsigned long long)g_rk[1],(unsigned long long)g_rk[0]);
  for(int s=0;s<ns3;s++) if(g_line[s]||g_all[s]) printf("support %d %d %d: line_sets=%llu all_sets=%llu\n",S3[s][0],S3[s][1],S3[s][2],(unsigned long long)g_line[s],(unsigned long long)g_all[s]);
  return 0; }
