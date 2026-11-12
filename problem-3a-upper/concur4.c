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

// concur4: per support, lines (2 equations) and planes (1 equation) in P^3 with multiplicities.
// Candidate words = full-support intersection points (line-line, line-plane, derived-line x line/plane);
// derived lines = plane-plane. Weight(w) = sum of mult of all spaces containing w. Prints heavy items.
typedef struct { int ne; fe e[2][4]; int mult; } sp;
static int rref(fe M[][4], int n, int *piv){ int r=0; for(int c=0;c<4&&r<n;c++){ int p=-1; for(int i=r;i<n;i++) if(!fz(M[i][c])){p=i;break;} if(p<0) continue;
  for(int q=0;q<4;q++){ fe t=M[r][q]; M[r][q]=M[p][q]; M[p][q]=t; } fe iv=finv(M[r][c]); for(int q=0;q<4;q++) M[r][q]=fmul(M[r][q],iv);
  for(int i=0;i<n;i++) if(i!=r&&!fz(M[i][c])){ fe t=M[i][c]; for(int q=0;q<4;q++) M[i][q]=fsub(M[i][q],fmul(t,M[r][q])); } piv[r++]=c; } return r; }
static int kerpt(fe M[][4], int n, fe *w){ int piv[4]; int r=rref(M,n,piv); if(r!=3) return 0; int isp[4]={0}; for(int i=0;i<3;i++) isp[piv[i]]=1; int f=0; while(isp[f]) f++;
  for(int q=0;q<4;q++) w[q]=(fe){0,0}; w[f]=(fe){1,0}; for(int i=0;i<3;i++){ fe t=M[i][f]; w[piv[i]]=(fe){(P-t.r)%P,(P-t.i)%P}; }
  for(int q=0;q<4;q++) if(fz(w[q])) return 0; fe iv=finv(w[0]); for(int q=0;q<4;q++) w[q]=fmul(w[q],iv); return 1; }
static int inside(const sp *k, const fe *w){ for(int i=0;i<k->ne;i++){ fe s={0,0}; for(int q=0;q<4;q++) s=fadd(s,fmul(k->e[i][q],w[q])); if(!fz(s)) return 0; } return 1; }
#define MAXS 400000
static sp L[MAXS]; static int nl;
#define HB (1u<<22)
static fe hw[HB][4]; static char hu[HB]; static int hlist[HB], nh;
static u64 hk(const fe *w){ u64 h=1469598103934665603ULL; for(int q=1;q<4;q++){ h=(h^w[q].r)*1099511628211ULL; h=(h^w[q].i)*1099511628211ULL; } return h; }
static void addpt(const fe *w){ u64 h=hk(w)&(HB-1); while(hu[h]){ if(!memcmp(hw[h],w,sizeof(fe)*4)) return; h=(h+1)&(HB-1);} hu[h]=1; memcpy(hw[h],w,sizeof(fe)*4); if(nh<(int)HB-1) hlist[nh++]=h; }
static int S4[1820][4]; static int THR;
static void process(int s){
  // derived lines from plane pairs
  int np0=nl; for(int i=0;i<np0;i++) if(L[i].ne==1) for(int j=i+1;j<np0;j++) if(L[j].ne==1){ fe M[4][4]; memcpy(M[0],L[i].e[0],sizeof M[0]); memcpy(M[1],L[j].e[0],sizeof M[1]); int piv[4]; if(rref(M,2,piv)!=2) continue;
      sp d; d.ne=2; memcpy(d.e[0],M[0],sizeof M[0]); memcpy(d.e[1],M[1],sizeof M[1]); d.mult=0; int dup=0; for(int k=np0;k<nl;k++) if(!memcmp(L[k].e,d.e,sizeof d.e)) dup=1; if(!dup && nl<MAXS) L[nl++]=d; }
  // derived-line weights
  for(int k=np0;k<nl;k++){ int wsum=0; for(int i=0;i<np0;i++){ int ok=1; if(L[i].ne==1){ /* plane contains line iff its eq is in span: test 2 kernel vecs */ }
      fe M[4][4]; memcpy(M[0],L[k].e[0],sizeof M[0]); memcpy(M[1],L[k].e[1],sizeof M[1]); for(int t=0;t<L[i].ne;t++) memcpy(M[2+t],L[i].e[t],sizeof M[0]); int piv[4]; ok=(rref(M,2+L[i].ne,piv)==2); if(ok) wsum+=L[i].mult; }
    if(wsum>=THR) printf("DLINE %d %d\n",s,wsum); }
  // candidate points
  memset(hu,0,HB); nh=0;
  for(int i=0;i<nl;i++) for(int j=i+1;j<nl;j++){ if(L[i].ne+L[j].ne<3) continue; fe M[4][4]; int n=0; for(int t=0;t<L[i].ne;t++) memcpy(M[n++],L[i].e[t],sizeof M[0]); for(int t=0;t<L[j].ne;t++) memcpy(M[n++],L[j].e[t],sizeof M[0]);
      fe w[4]; if(kerpt(M,n,w)) addpt(w); }
  int best=0;
  for(int a=0;a<nh;a++){ fe *w=hw[hlist[a]]; int wt=0; for(int i=0;i<np0;i++) if(inside(&L[i],w)) wt+=L[i].mult;
    if(wt>best) best=wt;
    if(wt>=THR){ printf("PT %d %d",s,wt); for(int q=1;q<4;q++) printf(" %llu %llu",(unsigned long long)w[q].r,(unsigned long long)w[q].i); printf("\n"); } }
  int ml=0; for(int i=0;i<np0;i++) if(L[i].mult>ml) ml=L[i].mult;
  printf("SUP %d spaces=%d points=%d maxpointweight=%d maxspace=%d\n",s,np0,nh,best,ml); fflush(stdout);
}
int main(int argc,char **argv){ THR=atoi(argv[2]); FILE *f=fopen(argv[1],"r"); char line[1024]; int cur=-1;
  static struct { int s; sp k; } *all; all=malloc(sizeof(*all)*1100000); int na=0;
  while(fgets(line,sizeof line,f)){ char *p=line; int m,s,d,nn; sscanf(p,"%d %d %d%n",&m,&s,&d,&nn); p+=nn; sp k; k.ne=4-d; k.mult=m;
    for(int i=0;i<k.ne;i++) for(int q=0;q<4;q++){ unsigned long long r,im; sscanf(p,"%llu %llu%n",&r,&im,&nn); p+=nn; k.e[i][q].r=r; k.e[i][q].i=im; }
    all[na].s=s; all[na].k=k; na++; }
  // sort by support (counting)
  static int cnt[1821]; for(int i=0;i<na;i++) cnt[all[i].s+1]++; for(int i=1;i<=1820;i++) cnt[i]+=cnt[i-1];
  int *ord=malloc(sizeof(int)*na); int pos[1821]; memcpy(pos,cnt,sizeof pos); for(int i=0;i<na;i++) ord[pos[all[i].s]++]=i;
  int sstart=atoi(argv[3]), sstep=atoi(argv[4]);
  for(int s=sstart;s<1820;s+=sstep){ if(cnt[s+1]==cnt[s]) continue; nl=0; for(int t=cnt[s];t<cnt[s+1];t++) L[nl++]=all[ord[t]].k; process(s); }
  (void)cur; return 0; }
