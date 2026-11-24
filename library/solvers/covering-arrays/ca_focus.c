/* ca_focus.c - aimed search for covering arrays CA(N; 3, k, v), for the endgame (few uncovered tuples).
 *
 * Instead of changing random cells, every main move is aimed at a combination that is still missing:
 *   1. pick an uncovered (column-triple, tuple);
 *   2. for every row, measure what writing that tuple into the row would cost (how many other tuples
 *      it would uncover) and write it into the cheapest row that is not tabu;
 *   3. in between, take a short random walk over cells that cost nothing, so the slack moves around.
 * Cost = number of uncovered (triple, tuple) pairs; 0 means a covering array.
 *
 * usage: ca_focus N k v init.txt out.txt seed seconds [walk=300] [T=0.20] [tenure=6] [maxup=99] [slack=9999]
 *   maxup: an aimed move is skipped if it would raise the cost by more than this
 *   slack: when the cost drifts this far above the best seen, go back to the best array
 * The init file may differ from the target by one step, so a passing array can be pushed further:
 *   N+1 rows   -> the row that covers the fewest tuples on its own is dropped (row descent)
 *   N-1 rows   -> one random row is added
 *   k-1 columns-> one balanced random column is added (column extension)
 * Strength t is fixed at 3. Prints progress on stderr and one summary line on stdout.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int N, K, V, NT, V3;          /* NT = number of column triples, V3 = v^3 */
static int *A;                       /* N x K */
static int (*tc)[3];                 /* columns of each triple */
static int *colq, *colz, cpc;        /* per column: triples through it and the column's position */
static int16_t *cnt;                 /* NT x V3 coverage counts */
static int *unc, ucnt, *upos;        /* set of uncovered keys (q*V3+tuple) */
static int pw[3];
static uint64_t rs;

static inline uint64_t rnd64(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline int rnd(int n) { return (int)(rnd64() % (uint64_t)n); }
static inline double rnd01(void) { return (rnd64() >> 11) * (1.0 / 9007199254740992.0); }

static inline void uadd(int key) { upos[key] = ucnt; unc[ucnt++] = key; }
static inline void udel(int key) { int p = upos[key], last = unc[--ucnt]; unc[p] = last; upos[last] = p; upos[key] = -1; }

static inline void setcell(int i, int j, int nw) {
    int *row = A + (size_t)i * K, old = row[j];
    if (old == nw) return;
    const int *qs = colq + (size_t)j * cpc, *zs = colz + (size_t)j * cpc;
    for (int m = 0; m < cpc; m++) {
        int q = qs[m];
        const int *c = tc[q];
        int oi = row[c[0]] * pw[0] + row[c[1]] * pw[1] + row[c[2]] * pw[2];
        int ni = oi + (nw - old) * pw[zs[m]];
        int16_t *cq = cnt + (size_t)q * V3;
        if (--cq[oi] == 0) uadd(q * V3 + oi);
        if (cq[ni]++ == 0) udel(q * V3 + ni);
    }
    row[j] = nw;
}

static void die(const char *m) { fprintf(stderr, "ca_focus: %s\n", m); exit(2); }

int main(int argc, char **argv) {
    if (argc < 8) die("usage: ca_focus N k v init.txt out.txt seed seconds [walk] [T] [tenure]");
    N = atoi(argv[1]); K = atoi(argv[2]); V = atoi(argv[3]);
    const char *initf = argv[4], *outf = argv[5];
    rs = (uint64_t)atoll(argv[6]) * 2654435761ULL + 88172645463325252ULL; if (!rs) rs = 42;
    double secs = atof(argv[7]);
    int walk = argc > 8 ? atoi(argv[8]) : 300;
    double T = argc > 9 ? atof(argv[9]) : 0.20;
    int tenure = argc > 10 ? atoi(argv[10]) : 6;
    int maxup = argc > 11 ? atoi(argv[11]) : 99;
    int slack = argc > 12 ? atoi(argv[12]) : 9999;
    V3 = V * V * V; pw[0] = V * V; pw[1] = V; pw[2] = 1;
    NT = K * (K - 1) * (K - 2) / 6; cpc = (K - 1) * (K - 2) / 2;
    A = malloc(sizeof(int) * (size_t)N * K);
    tc = malloc(sizeof(int[3]) * (size_t)NT);
    colq = malloc(sizeof(int) * (size_t)K * cpc); colz = malloc(sizeof(int) * (size_t)K * cpc);
    cnt = calloc((size_t)NT * V3, sizeof(int16_t));
    unc = malloc(sizeof(int) * (size_t)NT * V3); upos = malloc(sizeof(int) * (size_t)NT * V3);
    int *fill = calloc((size_t)K, sizeof(int)), *best = malloc(sizeof(int) * (size_t)N * K);
    long *tabu = calloc((size_t)N, sizeof(long));
    if (!A || !tc || !colq || !colz || !cnt || !unc || !upos || !fill || !best || !tabu) die("out of memory");
    FILE *f = fopen(initf, "r"); if (!f) die("cannot open init file");
    int *raw = calloc((size_t)(N + 1) * K, sizeof(int)), rows_in = 0, cols_in = 0;
    if (!raw) die("out of memory");
    { char line[1 << 16];
      while (fgets(line, sizeof line, f)) {
          int c = 0, val, used; char *p = line;
          while (sscanf(p, "%d%n", &val, &used) == 1) {
              if (val < 0 || val >= V) die("bad symbol in init file");
              if (rows_in > N || c >= K) die("init file is too large for N and k");
              raw[(size_t)rows_in * K + c++] = val; p += used; }
          if (!c) continue;
          if (!cols_in) cols_in = c; else if (c != cols_in) die("rows of different length in init file");
          rows_in++; } }
    fclose(f);
    if ((cols_in != K && cols_in != K - 1) || rows_in < N - 1 || rows_in > N + 1) die("init file must be within one row and one column of the target");
    if (cols_in == K - 1) {                 /* add a balanced random last column */
        for (int i = 0; i < rows_in; i++) raw[(size_t)i * K + K - 1] = i % V;
        for (int i = rows_in - 1; i > 0; i--) { int j = rnd(i + 1), tmp = raw[(size_t)i * K + K - 1];
            raw[(size_t)i * K + K - 1] = raw[(size_t)j * K + K - 1]; raw[(size_t)j * K + K - 1] = tmp; }
        fprintf(stderr, "added column %d\n", K - 1);
    }
    if (rows_in == N - 1) { for (int j = 0; j < K; j++) raw[(size_t)rows_in * K + j] = rnd(V); rows_in++; fprintf(stderr, "added one random row\n"); }
    int drop = -1;
    int q = 0;
    for (int a = 0; a < K; a++) for (int b = a + 1; b < K; b++) for (int c = b + 1; c < K; c++) {
        tc[q][0] = a; tc[q][1] = b; tc[q][2] = c;
        colq[(size_t)a * cpc + fill[a]] = q; colz[(size_t)a * cpc + fill[a]++] = 0;
        colq[(size_t)b * cpc + fill[b]] = q; colz[(size_t)b * cpc + fill[b]++] = 1;
        colq[(size_t)c * cpc + fill[c]] = q; colz[(size_t)c * cpc + fill[c]++] = 2;
        q++;
    }
    if (rows_in == N + 1) {                 /* drop the row with the fewest uniquely covered tuples */
        for (int i = 0; i <= N; i++) { const int *row = raw + (size_t)i * K;
            for (q = 0; q < NT; q++) cnt[(size_t)q * V3 + row[tc[q][0]] * pw[0] + row[tc[q][1]] * pw[1] + row[tc[q][2]] * pw[2]]++; }
        int bu = 1 << 30;
        for (int i = 0; i <= N; i++) { const int *row = raw + (size_t)i * K; int u = 0;
            for (q = 0; q < NT; q++) u += cnt[(size_t)q * V3 + row[tc[q][0]] * pw[0] + row[tc[q][1]] * pw[1] + row[tc[q][2]] * pw[2]] == 1;
            if (u < bu) { bu = u; drop = i; } }
        memset(cnt, 0, sizeof(int16_t) * (size_t)NT * V3);
        fprintf(stderr, "dropped row %d (%d tuples only it covered)\n", drop, bu);
    }
    for (int i = 0, o = 0; i < rows_in; i++) { if (i == drop) continue; memcpy(A + (size_t)o++ * K, raw + (size_t)i * K, sizeof(int) * (size_t)K); }
    free(raw);
    for (int i = 0; i < N; i++) { const int *row = A + (size_t)i * K;
        for (q = 0; q < NT; q++) cnt[(size_t)q * V3 + row[tc[q][0]] * pw[0] + row[tc[q][1]] * pw[1] + row[tc[q][2]] * pw[2]]++; }
    ucnt = 0;
    for (int key = 0; key < NT * V3; key++) { upos[key] = -1; if (!cnt[key]) uadd(key); }
    int bestcost = ucnt, start = ucnt;
    memcpy(best, A, sizeof(int) * (size_t)N * K);
    time_t t0 = time(NULL), tlast = t0;
    long it = 0, focused = 0;
    double pacc1 = exp(-1.0 / T);
    while (ucnt > 0) {
        it++;
        if ((it & 15) == 0) { time_t now = time(NULL);
            if (now - t0 >= secs) break;
            if (now - tlast >= 10) { tlast = now; fprintf(stderr, "t=%lds cost=%d best=%d focused=%ld\n", (long)(now - t0), ucnt, bestcost, focused); } }
        /* random walk over cells: free and improving changes always, a one-step worsening rarely */
        for (int m = 0; m < walk && ucnt > 0; m++) {
            int i = rnd(N), j = rnd(K), old = A[(size_t)i * K + j], nw = rnd(V - 1); if (nw >= old) nw++;
            int c0 = ucnt; setcell(i, j, nw);
            int d = ucnt - c0;
            if (d > 0 && !(d == 1 && rnd01() < pacc1)) setcell(i, j, old);
        }
        if (ucnt == 0) break;
        if (ucnt < bestcost) { bestcost = ucnt; memcpy(best, A, sizeof(int) * (size_t)N * K); }
        /* aimed move: write one missing tuple into the row where it does the least damage */
        int key = unc[rnd(ucnt)]; q = key / V3; int tup = key % V3;
        int tv[3] = { tup / pw[0], (tup / pw[1]) % V, tup % V };
        int brow = -1, bc = 1 << 30, ties = 0;
        for (int i = 0; i < N; i++) {
            if (tabu[i] > it) continue;
            int *row = A + (size_t)i * K, o[3] = { row[tc[q][0]], row[tc[q][1]], row[tc[q][2]] };
            for (int z = 0; z < 3; z++) setcell(i, tc[q][z], tv[z]);
            int c = ucnt;
            for (int z = 2; z >= 0; z--) setcell(i, tc[q][z], o[z]);
            if (c < bc) { bc = c; brow = i; ties = 1; }
            else if (c == bc && rnd(++ties) == 0) brow = i;
        }
        if (brow < 0 || bc - ucnt > maxup) continue;
        for (int z = 0; z < 3; z++) setcell(brow, tc[q][z], tv[z]);
        tabu[brow] = it + tenure + rnd(tenure + 1);
        focused++;
        if (ucnt < bestcost) { bestcost = ucnt; memcpy(best, A, sizeof(int) * (size_t)N * K); }
        if (ucnt > bestcost + slack) {          /* drifted too far: return to the best array */
            for (int i = 0; i < N; i++) for (int j = 0; j < K; j++) setcell(i, j, best[(size_t)i * K + j]);
            memset(tabu, 0, sizeof(long) * (size_t)N);
        }
    }
    if (ucnt == 0) { bestcost = 0; memcpy(best, A, sizeof(int) * (size_t)N * K); }
    f = fopen(outf, "w"); if (!f) die("cannot write output");
    for (int i = 0; i < N; i++) { for (int j = 0; j < K; j++) fprintf(f, "%d%c", best[(size_t)i * K + j], j + 1 < K ? ' ' : '\n'); }
    fclose(f);
    printf("CA(%d;3,%d,%d): start=%d uncovered_tuples=%d focused_moves=%ld seconds=%ld array_written=%s\n",
           N, K, V, start, bestcost, focused, (long)(time(NULL) - t0), outf);
    return bestcost == 0 ? 0 : 1;
}
