/* ca_focus_t: aimed local search for covering arrays of any strength t (2..6).
 * usage: ca_focus_t t N k v init.txt out.txt seed seconds [walk=300] [T=0.20] [tenure=6]
 * init may have N+1 rows (the row covering the fewest tuples alone is dropped), N rows, or N-1 rows (a random row is
 * added). Each step picks a random uncovered (column-set, value-tuple), scores every row by how many tuples it would
 * lose/gain if forced to that tuple (cells changed one by one, counts updated exactly), takes the best non-tabu row,
 * and accepts worse moves with probability exp(-delta/T). Every `walk` steps without improvement it restarts from the
 * best array. Writes the best array to out.txt and prints:
 *   CA(N;t,k,v): start=S uncovered_tuples=U moves=M seconds=X array_written=out.txt
 * The output still has to pass the frozen checker; this program certifies nothing. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

static int t, N, K, V, NC, VT;           /* NC = C(K,t), VT = V^t */
static int *A;                           /* N x K */
static int *cols;                        /* NC x t */
static int **ccombo, **cpos, *ccnt;      /* for each column: combos containing it and its position */
static int pw[8];
static int16_t *cnt;                     /* NC x VT */
static int *unc, *upos, nunc;            /* uncovered list of (combo*VT+code), upos[-1 if covered] */
static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static void unc_add(int id) { if (upos[id] < 0) { upos[id] = nunc; unc[nunc++] = id; } }
static void unc_del(int id) { int p = upos[id]; if (p >= 0) { int last = unc[--nunc]; unc[p] = last; upos[last] = p; upos[id] = -1; } }

static int code_of(int r, int c) { int s = 0; for (int j = 0; j < t; j++) s += A[r * K + cols[c * t + j]] * pw[t - 1 - j]; return s; }

/* change cell (r,col) to b, updating counts and the uncovered list; returns change in #uncovered */
static int set_cell(int r, int col, int b) {
    int a = A[r * K + col]; if (a == b) return 0;
    int d = 0;
    for (int i = 0; i < ccnt[col]; i++) {
        int c = ccombo[col][i], p = cpos[col][i];
        int oc = code_of(r, c), nc = oc + (b - a) * pw[t - 1 - p];
        int io = c * VT + oc, in = c * VT + nc;
        if (--cnt[io] == 0) { unc_add(io); d++; }
        if (cnt[in]++ == 0) { unc_del(in); d--; }
    }
    A[r * K + col] = b;
    return d;
}

static void rebuild(void) {
    memset(cnt, 0, sizeof(int16_t) * (size_t)NC * VT);
    for (int r = 0; r < N; r++) for (int c = 0; c < NC; c++) cnt[c * VT + code_of(r, c)]++;
    nunc = 0; for (int i = 0; i < NC * VT; i++) { upos[i] = -1; if (!cnt[i]) unc_add(i); }
}

int main(int argc, char **argv) {
    if (argc < 9) { fprintf(stderr, "usage: ca_focus_t t N k v init out seed seconds [walk T tenure]\n"); return 2; }
    t = atoi(argv[1]); N = atoi(argv[2]); K = atoi(argv[3]); V = atoi(argv[4]);
    const char *initf = argv[5], *outf = argv[6];
    rs = (uint64_t)atoll(argv[7]) * 2654435761ULL + 88172645463325252ULL; if (!rs) rs = 42;
    double secs = atof(argv[8]);
    int walk = argc > 9 ? atoi(argv[9]) : 300; double T = argc > 10 ? atof(argv[10]) : 0.20; int tenure = argc > 11 ? atoi(argv[11]) : 6;
    if (t < 2 || t > 6 || t > K) { fprintf(stderr, "bad t\n"); return 2; }
    pw[0] = 1; for (int i = 1; i <= t; i++) pw[i] = pw[i - 1] * V; VT = pw[t];
    /* enumerate combos */
    double nc = 1; for (int i = 0; i < t; i++) nc = nc * (K - i) / (i + 1); NC = (int)(nc + 0.5);
    if ((double)NC * VT > 4e8) { fprintf(stderr, "too large: %d combos x %d tuples\n", NC, VT); return 2; }
    cols = malloc(sizeof(int) * (size_t)NC * t);
    int idx[8]; for (int i = 0; i < t; i++) idx[i] = i;
    for (int c = 0; c < NC; c++) {
        for (int i = 0; i < t; i++) cols[c * t + i] = idx[i];
        int i = t - 1; while (i >= 0 && idx[i] == K - t + i) i--;
        if (i < 0) break; idx[i]++; for (int j = i + 1; j < t; j++) idx[j] = idx[j - 1] + 1;
    }
    ccnt = calloc(K, sizeof(int)); ccombo = malloc(sizeof(int *) * K); cpos = malloc(sizeof(int *) * K);
    int per = (int)(nc * t / K + 0.5);
    for (int k = 0; k < K; k++) { ccombo[k] = malloc(sizeof(int) * per); cpos[k] = malloc(sizeof(int) * per); }
    for (int c = 0; c < NC; c++) for (int j = 0; j < t; j++) { int k = cols[c * t + j]; ccombo[k][ccnt[k]] = c; cpos[k][ccnt[k]++] = j; }
    cnt = malloc(sizeof(int16_t) * (size_t)NC * VT); unc = malloc(sizeof(int) * (size_t)NC * VT); upos = malloc(sizeof(int) * (size_t)NC * VT);
    /* read init */
    FILE *f = fopen(initf, "r"); if (!f) { fprintf(stderr, "ca_focus_t: cannot open init file\n"); return 2; }
    int cap = (N + 2) * K, *raw = malloc(sizeof(int) * cap), n = 0, x;
    while (n < cap && fscanf(f, "%d", &x) == 1) raw[n++] = x;
    fclose(f);
    int rows = n / K; if (rows * K != n || rows < N - 1 || rows > N + 1) { fprintf(stderr, "init shape: %d values, k=%d\n", n, K); return 2; }
    for (int i = 0; i < n; i++) if (raw[i] < 0 || raw[i] >= V) { fprintf(stderr, "init symbol out of range\n"); return 2; }
    A = malloc(sizeof(int) * (size_t)(N + 1) * K);
    if (rows == N + 1) {          /* drop the row that alone covers the fewest tuples */
        int saveN = N; N = N + 1; memcpy(A, raw, sizeof(int) * n); rebuild();
        int best = 0, bs = 1 << 30;
        for (int r = 0; r < N; r++) { int s = 0; for (int c = 0; c < NC; c++) if (cnt[c * VT + code_of(r, c)] == 1) s++; if (s < bs) { bs = s; best = r; } }
        memmove(A + best * K, A + (best + 1) * K, sizeof(int) * (size_t)(N - best - 1) * K); N = saveN;
    } else {
        memcpy(A, raw, sizeof(int) * n);
        if (rows == N - 1) for (int c = 0; c < K; c++) A[(N - 1) * K + c] = rnd() % V;
    }
    rebuild();
    int start = nunc, bestU = nunc, *B = malloc(sizeof(int) * (size_t)N * K); memcpy(B, A, sizeof(int) * (size_t)N * K);
    long *tabu = calloc(N, sizeof(long)); long moves = 0, since = 0;
    int *save = malloc(sizeof(int) * t);
    clock_t t0 = clock();
    while (nunc > 0) {
        if ((moves & 255) == 0 && (double)(clock() - t0) / CLOCKS_PER_SEC > secs) break;
        moves++;
        int id = unc[rnd() % nunc], c = id / VT, code = id % VT, want[8];
        for (int j = 0; j < t; j++) want[j] = (code / pw[t - 1 - j]) % V;
        int br = -1, bd = 1 << 30, ties = 0;
        for (int r = 0; r < N; r++) {
            if (tabu[r] > moves) continue;
            int d = 0;
            for (int j = 0; j < t; j++) { int col = cols[c * t + j]; save[j] = A[r * K + col]; d += set_cell(r, col, want[j]); }
            for (int j = t - 1; j >= 0; j--) set_cell(r, cols[c * t + j], save[j]);
            if (d < bd) { bd = d; br = r; ties = 1; } else if (d == bd && rnd() % ++ties == 0) br = r;
        }
        if (br < 0) continue;
        if (bd > 0 && urand() >= exp(-bd / T)) { since++; goto chk; }
        for (int j = 0; j < t; j++) set_cell(br, cols[c * t + j], want[j]);
        tabu[br] = moves + tenure;
        if (nunc < bestU) { bestU = nunc; memcpy(B, A, sizeof(int) * (size_t)N * K); since = 0; } else since++;
    chk:
        if (since > walk * 50) { memcpy(A, B, sizeof(int) * (size_t)N * K); rebuild(); since = 0; }
    }
    if (nunc < bestU || nunc == 0) { bestU = nunc; memcpy(B, A, sizeof(int) * (size_t)N * K); }
    FILE *o = fopen(outf, "w");
    for (int r = 0; r < N; r++) { for (int k = 0; k < K; k++) fprintf(o, k ? " %d" : "%d", B[r * K + k]); fputc('\n', o); }
    fclose(o);
    printf("CA(%d;%d,%d,%d): start=%d uncovered_tuples=%d moves=%ld seconds=%.0f array_written=%s\n",
           N, t, K, V, start, bestU, moves, (double)(clock() - t0) / CLOCKS_PER_SEC, outf);
    return 0;
}
