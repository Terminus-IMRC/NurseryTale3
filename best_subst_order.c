/*
 * Copyright (c) 2019 Sugizaki Yukimasa (sugizaki@hpcs.cs.tsukuba.ac.jp)
 * All rights reserved.
 *
 * This software is licensed under a Modified (3-Clause) BSD License.
 * You should have received a copy of this license along with this
 * software. If not, contact the copyright holder above.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>


#define __stringify_1(x...) #x
#define __stringify(x...) __stringify_1(x)
#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)


#if ORDER != 5 && ORDER != 6
#error "Choose ORDER from: 5, 6"
#endif

typedef unsigned __int128 score_t;

#if   ORDER == 5
typedef uint32_t binsquare_t;
#define popcount_bs __builtin_popcount
#define DEPTH_MAX 14
#define SCORE_BASE_SHIFT 5
#define BSLIST_LEN 1188905
#elif ORDER == 6
typedef uint64_t binsquare_t;
#define popcount_bs __builtin_popcountll
#define DEPTH_MAX 23
#define SCORE_BASE_SHIFT 5
#define BSLIST_LEN 0xdeadbeaf
#endif


static void print_binsquare(const binsquare_t binsquare)
{
    int i;
    printf("(");
    for (i = sizeof(binsquare) * 8 - 1; i >= ORDER*ORDER; i --)
        printf("%d", !!(binsquare & (((binsquare_t) 1) << i)));
    printf(")");
    for (i = ORDER*ORDER-1; i >= 0; i --)
        printf("%d", !!(binsquare & (((binsquare_t) 1) << i)));
    printf("\n");
}

static void* bslist_load(void)
{
    const char *filename = "bslist." __stringify(ORDER);
    struct stat sb;
    size_t rets;
    void *p;
    FILE *fp;

    if (lstat(filename, &sb)) {
        fprintf(stderr, "lstat: %s: %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (sb.st_size != BSLIST_LEN * sizeof(binsquare_t)) {
        fprintf(stderr, "%s: size does not match\n", filename);
        exit(EXIT_FAILURE);
    }

    if (posix_memalign(&p, 65536, BSLIST_LEN * sizeof(binsquare_t))) {
        fprintf(stderr, "posix_memalign: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "fopen: %s: %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    rets = fread(p, sizeof(binsquare_t), BSLIST_LEN, fp);
    if (rets != BSLIST_LEN) {
        fprintf(stderr, "fread: %s: %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (fclose(fp)) {
        fprintf(stderr, "fclose: %s: %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    return p;
}

static binsquare_t *bslist;

#define ALONE_1(x) (((x) & ((x) - 1)) ? 0 : (x))
#define RIGHTMOST_0(x) (~(x) & ((x) + 1))
#define NEXT_BF(bf, filled) RIGHTMOST_0((filled) | (((bf) << 1) - 1))

static score_t best_order_recurse(binsquare_t filled, const unsigned depth)
{
    size_t i;

    binsquare_t bf, obv;
    binsquare_t best_sets_filled[ORDER*ORDER];
    size_t best_sets_len = 0;
    unsigned best_sets_score_base_best = 0; /* Score based on only this call. */

    for (bf = RIGHTMOST_0(filled); bf != 0;
            bf = NEXT_BF(bf, filled), obv = (binsquare_t) 0) {

        filled |= bf;

retry:
        for (i = 0; i < BSLIST_LEN; i ++) {
            binsquare_t tmp = ALONE_1(~filled & bslist[i]);
            if (tmp) {
                filled |= tmp;
                obv |= tmp;
                goto retry;
            }
        }

        const unsigned score_base = 1 + popcount_bs(obv);

        if (filled == (binsquare_t) -1) {
            assert(depth == DEPTH_MAX - 1);
            return (score_t) score_base << (SCORE_BASE_SHIFT * (DEPTH_MAX - 1 - depth));
        }

        if (score_base > best_sets_score_base_best) {
            best_sets_filled[0] = filled;
            best_sets_len = 1;
            best_sets_score_base_best = score_base;
        } else if (score_base == best_sets_score_base_best)
            best_sets_filled[best_sets_len++] = filled;

        filled ^= bf | obv;
    }

    assert(best_sets_len != 0);

    /*
     * If it is at the bottom of the recursion and has not yet returned,
     * it means that filling in depth=depth_max is failed.
     */
    if (depth >= DEPTH_MAX - 1)
        return 0;

    if (depth <= 10)
        printf("depth=%2u: len=%2zu, best=%2d\n",
                depth, best_sets_len, best_sets_score_base_best);

    score_t score_base_children_best = 0;
    for (i = 0; i < best_sets_len; i ++) {
        score_t score_base_children;
        filled = best_sets_filled[i];

        score_base_children = best_order_recurse(filled, depth + 1);
        if (score_base_children > score_base_children_best)
            score_base_children_best = score_base_children;
    }

    return (best_sets_score_base_best
            << (SCORE_BASE_SHIFT * (DEPTH_MAX - 1 - depth)))
                    | score_base_children_best;
}

int main(void)
{
    unsigned i;
    binsquare_t filled = ~((1 << (ORDER*ORDER)) - 1);

    bslist = bslist_load();

    score_t score = best_order_recurse(filled, 0);

    printf("Result:\n");
    for (i = 0; i < DEPTH_MAX; i ++, score >>= SCORE_BASE_SHIFT)
        printf("depth=%2u: score_base=%u\n",
                DEPTH_MAX - 1 - i, (uint8_t) (score & ((1 << SCORE_BASE_SHIFT) - 1)));

    return 0;
}
