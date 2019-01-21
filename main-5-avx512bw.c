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


#if defined(ORDER) && ORDER != 5
#error "This program is for ORDER=5 only"
#endif
#define ORDER 5

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

typedef int_fast8_t square_elem_t;
typedef __m256i square_t;
typedef uint32_t binsquare_t;

#define BINSQUARE_MAP_SIZE (((size_t) 1) << (ORDER*ORDER - 3))

static void print_square(square_t square)
{
    printf("(%2d ", (int8_t) _mm256_extract_epi8(square, 31));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 30));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 29));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 28));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 27));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 26));
    printf("%2d) ", (int8_t) _mm256_extract_epi8(square, 25));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 24));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 23));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 22));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 21));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 20));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 19));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 18));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 17));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 16));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 15));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 14));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 13));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 12));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 11));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 10));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 9));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 8));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 7));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 6));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 5));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 4));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 3));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 2));
    printf("%2d ", (int8_t) _mm256_extract_epi8(square, 1));
    printf("%2d\n", (int8_t)_mm256_extract_epi8(square, 0));
}

#define get_add(line_id, c) \
    ({ \
        square_t __attribute__((aligned(32))) adds[ORDER*2+2] = { \
                _mm256_set_epi8(0,0,0,0,0,0,0,c,c,c,c,c,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
                _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), \
                _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c,0,0,0,0,0,0,0,0,0,0), \
                _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c,0,0,0,0,0), \
                _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,c,c,c,c,c), \
                _mm256_set_epi8(0,0,0,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0), \
                _mm256_set_epi8(0,0,0,0,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0), \
                _mm256_set_epi8(0,0,0,0,0,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0), \
                _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0), \
                _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c,0,0,0,0,c), \
                _mm256_set_epi8(0,0,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c,0,0,0,0,0,c), \
                _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,0,c,0,0,0,c,0,0,0,c,0,0,0,c,0,0,0,c,0,0,0,0), \
        }; \
        adds[line_id]; \
    })

#define square_addsub_line(square, line_id, c) \
    square = _mm256_add_epi8(square, get_add(line_id, c))


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

#if ORDER <= 5
#warning "Using in-memory index"
     map = my_cvalloc(BINSQUARE_MAP_SIZE);
     if (map == NULL) {
         fprintf(stderr, "my_cvalloc: %s\n", strerror(errno));
         exit(EXIT_FAILURE);
     }
#else
#warning "Using out-of-memory index"
     int fd;
     const char *filename = PREFIX "index-" STR(ORDER) ".db";
     int err;

     fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
     if (fd < 0) {
         fprintf(stderr, "open: %s\n", strerror(errno));
         exit(EXIT_FAILURE);
     }

     err = unlink(filename);
     if (err) {
         fprintf(stderr, "unlink: %s\n", strerror(errno));
         exit(EXIT_FAILURE);
     }

     err = ftruncate(fd, size);
     if (err) {
         fprintf(stderr, "ftruncate: %s\n", strerror(errno));
         exit(EXIT_FAILURE);
     }

     map = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
     if (map == MAP_FAILED) {
         fprintf(stderr, "mmap: %s\n", strerror(errno));
         exit(EXIT_FAILURE);
     }

     if (close(fd)) {
         fprintf(stderr, "close: %s\n", strerror(errno));
         exit(EXIT_FAILURE);
     }
#endif

     return map;
}

static void binsquare_finalize(void *map)
{
    FILE *fp;
    size_t rets;
    const int tid = omp_get_thread_num();
    char str[0x100];

    snprintf(str, sizeof(str), "bsmap.%d.%d.%d.%d",
            ORDER, COEFF_MIN, COEFF_MAX, tid);

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

#pragma omp parallel
    {
        square_t square;
        uint64_t *map = binsquare_init();
        int c0, c1, c2;

#define X(n) int c##n; square_t square_orig_c##n;
            X(3)
            X(4)
            X(5)
            X(6)
            X(7)
            X(8)
            X(9)
            X(10)
            X(11)
#undef X

#pragma omp for nowait collapse(3)
        for (c0 = COEFF_MIN; c0 <= COEFF_MAX; c0 ++) {
            for (c1 = COEFF_MIN; c1 <= COEFF_MAX; c1 ++) {
                for (c2 = COEFF_MIN; c2 <= COEFF_MAX; c2 ++) {
                    square = _mm256_setzero_si256();
                    ADD(0, c0);
                    ADD(1, c1);
                    ADD(2, c2);
                    STA(3);
                        STA(4);
                            STA(5);
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
                                                    __mmask32 mask = _mm256_cmpneq_epi8_mask(square, _mm256_setzero_si256());
                                                    for (c11 = COEFF_MIN+1; c11 <= COEFF_MAX; c11 ++) {
                                                        ADD(11, 1);
                                                        const binsquare_t binsquare = _cvtmask32_u32(mask);
                                                        mask = _mm256_cmpneq_epi8_mask(square, _mm256_setzero_si256());
                                                        const size_t off = binsquare >> 6;
                                                        const uint64_t hot = ((uint64_t) 1) << (binsquare & ((binsquare_t) (64-1)));
                                                        if (!(map[off] & hot))
                                                            map[off] |= hot;
                                                    }
                                                    const binsquare_t binsquare = _cvtmask32_u32(mask);
                                                    const size_t off = binsquare >> 6;
                                                    const uint64_t hot = ((uint64_t) 1) << (binsquare & ((binsquare_t) (64-1)));
                                                    if (!(map[off] & hot))
                                                        map[off] |= hot;
                                                    POP(11);
#endif
                                                END(10);
                                            END(9);
                                        END(8);
                                    END(7);
                                END(6);
                            END(5);
                        END(4);
                    END(3);
                }
            }
        }

        printf("Final square:    ");
        print_square(square);

        binsquare_finalize(map);
    }

    {
        char str[0x100];
        snprintf(str, sizeof(str), "./bsmap_gather bsmap.%d.%d.%d.*",
                ORDER, COEFF_MIN, COEFF_MAX);
        (void) execl("/bin/sh", "sh", "-c", str, NULL);
    }

    return 0;
}
