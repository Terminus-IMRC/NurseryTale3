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
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>
#include <immintrin.h>
#include <omp.h>


#if defined(ORDER) && ORDER != 6
#error "This program is for ORDER=6 only"
#endif
#define ORDER 6

#if !defined(COEFF_MIN) || !defined(COEFF_MAX)
#error "Define COEFF_MIN and COEFF_MAX"
#endif

#if COEFF_MIN > COEFF_MAX
#error "COEFF_MIN must be smaller or equal to COEFF_MAX"
#endif

#ifndef PREFIX
#define PREFIX ""
#endif /* PREFIX */

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

/*
 * (1 << (8 * sizeof(square_elem_t) - 1)) must be larger than
 * maxabs(coeff_min, coeff_max) * (2*order + 2)
 */
typedef int_fast8_t square_elem_t;
typedef __m512i square_t;
typedef uint64_t binsquare_t;

#define BINSQUARE_MAP_SIZE (((size_t) 1) << (ORDER*ORDER - 3))

#define _mm512_extract_epi8(a, imm8) \
    _mm_extract_epi8(_mm512_extracti32x4_epi32(a, imm8 / (128 / 8)), imm8 % (128 / 8))

static void print_square(square_t square)
{
    printf("(%2d ", (int8_t) _mm512_extract_epi8(square, 63));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 62));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 61));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 60));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 59));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 58));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 57));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 56));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 55));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 54));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 53));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 52));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 51));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 50));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 49));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 48));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 47));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 46));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 45));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 44));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 43));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 42));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 41));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 40));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 39));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 38));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 37));
    printf("%2d) ", (int8_t) _mm512_extract_epi8(square, 36));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 35));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 34));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 33));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 32));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 31));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 30));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 29));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 28));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 27));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 26));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 25));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 24));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 23));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 22));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 21));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 20));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 19));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 18));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 17));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 16));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 15));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 14));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 13));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 12));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 11));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 10));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 9));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 8));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 7));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 6));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 5));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 4));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 3));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 2));
    printf("%2d ", (int8_t) _mm512_extract_epi8(square, 1));
    printf("%2d\n", (int8_t)_mm512_extract_epi8(square, 0));
}

#define get_add(line_id, c) \
    ({ \
        square_t __attribute__((aligned(32))) adds[ORDER*2+2] = { \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c,c,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c,c,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c,c,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c,c,0,0,0,0,0,0,0,0,0,0,0,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c,c,0,0,0,0,0,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c,c), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,0,0,c,0,0,0,0,0,0,c,0,0,0,0,0,0,c,0,0,0,0,0,0,c,0,0,0,0,0,0,c), \
                _mm512_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,0), \
        }; \
        adds[line_id]; \
    })

#define square_addsub_line(square, line_id, c) \
    square = _mm512_add_epi8(square, get_add(line_id, c))


static void* my_cvalloc(const size_t size)
{
    intptr_t p;
    const size_t align = 65536;
    const size_t size_suf = size + align - 1;
    p = (intptr_t) calloc(size_suf, 1);
    if (p == (intptr_t) NULL)
        return NULL;
    const size_t mod = ((intptr_t) p) % align;
    if (mod)
        p += align - mod;
    return (void*) p;
}

static void* binsquare_init(void)
{
     void *map;

     printf("Mapping %zu bytes\n", BINSQUARE_MAP_SIZE);

     /*
      * order size
      * 3     64B
      * 4     8KiB
      * 5     4MiB
      * 6     8GiB
      */

     map = my_cvalloc(BINSQUARE_MAP_SIZE);
     if (map == NULL) {
         fprintf(stderr, "my_cvalloc: %s\n", strerror(errno));
         exit(EXIT_FAILURE);
     }

     return map;
}

static void binsquare_finalize(void *map)
{
    FILE *fp;
    size_t rets;
    char str[0x100];

    snprintf(str, sizeof(str), "bsmap.%d.%d.%d", ORDER, COEFF_MIN, COEFF_MAX);

    fp = fopen(str, "wb");
    if (fp == NULL) {
        fprintf(stderr, "fopen: %s: %s\n", str, strerror(errno));
        exit(EXIT_FAILURE);
    }

    rets = fwrite(map, 1, BINSQUARE_MAP_SIZE, fp);
    if (rets != BINSQUARE_MAP_SIZE) {
        fprintf(stderr, "fwrite: %s: %s\n", str, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (fclose(fp)) {
        fprintf(stderr, "fclose: %s: %s\n", str, strerror(errno));
        exit(EXIT_FAILURE);
    }
}

int main(void)
{
    printf("Built on %s %s\n", __DATE__, __TIME__);
    printf("ORDER = %d\n", ORDER);
    printf("COEFF_{MIN,MAX} = {%d, %d}\n", COEFF_MIN, COEFF_MAX);

#define PUSH(n) square_orig_c##n = square

#define LOOP(n) for (c##n = COEFF_MIN; c##n <= COEFF_MAX; c##n ++)

#define ADD(n, c) square_addsub_line(square, n, c)

#define POP(n) square = square_orig_c##n

#define STA(n) PUSH(n); ADD(n, COEFF_MIN); LOOP(n) {
#define END(n) ADD(n, 1); } POP(n)

    uint64_t *map = binsquare_init();

#pragma omp parallel firstprivate(map)
    {
        square_t square;
        int c0, c1, c2, c3, c4, c5;

#define X(n) int c##n; square_t square_orig_c##n;
            X(6)
            X(7)
            X(8)
            X(9)
            X(10)
            X(11)
#undef X

#pragma omp for nowait collapse(6)
        for (c0 = 0; c0 <= COEFF_MAX; c0 ++) {
            for (c1 = COEFF_MIN; c1 <= COEFF_MAX; c1 ++) {
                for (c2 = COEFF_MIN; c2 <= COEFF_MAX; c2 ++) {
                    for (c3 = COEFF_MIN; c3 <= COEFF_MAX; c3 ++) {
                        for (c4 = COEFF_MIN; c4 <= COEFF_MAX; c4 ++) {
                            for (c5 = COEFF_MIN; c5 <= COEFF_MAX; c5 ++) {
                                square = _mm512_setzero_si512();
                                ADD(0, c0);
                                ADD(1, c1);
                                ADD(2, c2);
                                ADD(3, c3);
                                ADD(4, c4);
                                ADD(5, c5);
                                STA(6);
                                    STA(7);
                                        STA(8);
                                            STA(9);
                                                STA(10);
#if 0
                                                    STA(11);
                                                        const binsquare_t binsquare = _cvtmask32_u32(_mm256_cmpneq_epi8_mask(square, _mm256_setzero_si256()));
                                                        const size_t off = binsquare >> 6;
                                                        const uint64_t mask = ((uint64_t) 1) << (binsquare & ((binsquare_t) (64-1)));
                                                        if (unlikely(!(map[off] & mask))) {
                                                            map[off] |= mask;
                                                            innovative_count++;
                                                            //printf("inn: "); print_square(square);
                                                            //print_binsquare(binsquare);
                                                        }
                                                    END(11);
#else
                                                    PUSH(11);
                                                    ADD(11, COEFF_MIN);
                                                    __mmask64 mask = _mm512_cmpneq_epi8_mask(square, _mm512_setzero_si512());
                                                    for (c11 = COEFF_MIN+1; c11 <= COEFF_MAX; c11 ++) {
                                                        ADD(11, 1);
                                                        const binsquare_t binsquare = _cvtmask64_u64(mask);
                                                        mask = _mm512_cmpneq_epi8_mask(square, _mm512_setzero_si512());
                                                        const size_t off = binsquare >> 6;
                                                        const uint64_t hot = ((uint64_t) 1) << (binsquare & ((binsquare_t) (64-1)));
                                                        if (!(map[off] & hot)) {
#pragma omp critical
                                                            if (!(map[off] & hot))
                                                                map[off] |= hot;
                                                        }
                                                    }
                                                    const binsquare_t binsquare = _cvtmask64_u64(mask);
                                                    const size_t off = binsquare >> 6;
                                                    const uint64_t hot = ((uint64_t) 1) << (binsquare & ((binsquare_t) (64-1)));
                                                    if (!(map[off] & hot)) {
#pragma omp critical
                                                        if (!(map[off] & hot))
                                                            map[off] |= hot;
                                                    }
                                                    POP(11);
#endif
                                                END(10);
                                            END(9);
                                        END(8);
                                    END(7);
                                END(6);
                            }
                        }
                    }
                }
            }
        }

        printf("Final square:    ");
        print_square(square);
    }

    fflush(stdout);
    printf("Writing bsmap to file\n");
    binsquare_finalize(map);

    {
        char str[0x100];
        snprintf(str, sizeof(str), "./bsmap_gather bsmap.%d.%d.%d",
                ORDER, COEFF_MIN, COEFF_MAX);
        (void) execl("/bin/sh", "sh", "-c", str, NULL);
    }

    return 0;
}
