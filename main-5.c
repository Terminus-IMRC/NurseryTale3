/*
 * Copyright (c) 2019 Sugizaki Yukimasa (sugizaki@hpcs.cs.tsukuba.ac.jp)
 * All rights reserved.
 *
 * This software is licensed under a Modified (3-Clause) BSD License.
 * You should have received a copy of this license along with this
 * software. If not, contact the copyright holder above.
 */

#include <stdio.h>
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


#if !defined(COEFF_MIN) || !defined(COEFF_MAX)
#error "Define COEFF_MIN and COEFF_MAX"
#endif

#ifndef PREFIX
#define PREFIX ""
#endif /* PREFIX */

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

#define _STR(x) #x
#define STR(x) _STR(x)

typedef int_fast8_t square_elem_t;
typedef square_elem_t square_t [ORDER*ORDER];

#if   ORDER == 3 || ORDER == 4
typedef uint16_t binsquare_t;
#elif ORDER == 5
typedef uint32_t binsquare_t;
#elif ORDER == 6
typedef uint64_t binsquare_t;
#endif

static void *binsquare_map = NULL;
static const size_t binsquare_map_len = ((size_t) 1) << (ORDER*ORDER);
static uint32_t innovative_count = 0;

/*
 * 290
 * 4857
 * 55913
 * 557406
 */

/*
 * New records: order, depth, count
 * 3, 10, 399
 * 4, 12, 16143
 */

static void print_square(square_t square)
{
    int i;
    for (i = 0; i < ORDER*ORDER; i ++)
        printf("%2d%c", square[i], (i == ORDER*ORDER-1) ? '\n' : ' ');
}

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

static inline binsquare_t square_addsub_line(square_t square,
        const unsigned line_id, const square_elem_t c, binsquare_t binsquare)
{
    int i;
    assert(line_id <= ORDER*ORDER + 2);
    if (line_id < ORDER) {
        const unsigned offset = line_id * ORDER;
        binsquare &= ~(((((binsquare_t) 1) << ORDER) - 1) << offset);
        for (i = 0; i < ORDER; i ++) {
            const unsigned idx = offset + i;
            if ((square[idx] += c) != 0)
                binsquare |= (((binsquare_t) 1) << idx);
        }
    } else if (line_id < ORDER*2) {
        for (i = line_id - ORDER; i < ORDER*ORDER; i += ORDER) {
            if ((square[i] += c) != 0)
                binsquare |= (((binsquare_t) 1) << i);
            else
                binsquare &= ~(((binsquare_t) 1) << i);
        }
    } else if (line_id == ORDER*2 + 0) {
        for (i = 0; i < ORDER*ORDER; i += ORDER+1) {
            if ((square[i] += c) != 0)
                binsquare |= (((binsquare_t) 1) << i);
            else
                binsquare &= ~(((binsquare_t) 1) << i);
        }
    } else if (line_id == ORDER*2 + 1) {
        for (i = ORDER-1; i <= ORDER*(ORDER-1); i += ORDER-1) {
            if ((square[i] += c) != 0)
                binsquare |= (((binsquare_t) 1) << i);
            else
                binsquare &= ~(((binsquare_t) 1) << i);
        }
    } else {
        assert(0);
    }
    return binsquare;
}

static void binsquare_init(void)
{
     void *map;
     const size_t size = binsquare_map_len * sizeof(binsquare_t);

     printf("Mapping %zu bytes\n", size);

     /*
      * order size
      * 3     64B
      * 4     8KiB
      * 5     4MiB
      * 6     8GiB
      */

#if ORDER <= 5
#warning "Using in-memory index"
     if (posix_memalign(&map, 65536, size)) {
         fprintf(stderr, "posix_memalign: %s\n", strerror(errno));
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

     (void) memset(map, 0, size);
     binsquare_map = map;
}

static void binsquare_finalize(void)
{
    FILE *fp;
    size_t rets;

    fp = fopen("bsmap." STR(ORDER) "." STR(COEFF_MIN) "." STR(COEFF_MAX), "wb");
    if (fp == NULL) {
        fprintf(stderr, "fopen: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    rets = fwrite(binsquare_map, sizeof(binsquare_t), binsquare_map_len, fp);
    if (rets != binsquare_map_len) {
        fprintf(stderr, "fwrite: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (fclose(fp)) {
        fprintf(stderr, "fclose: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
}

int main(void)
{
    uint64_t *map;
    square_t square = {0};
    binsquare_t binsquare = (binsquare_t) 0;

#define X(n) int c##n; square_t square_orig_c##n; binsquare_t binsquare_orig_c##n;

    X(0)
    X(1)
    X(2)
    X(3)
    X(4)
    X(5)
    X(6)
    X(7)
    X(8)
    X(9)
    X(10)
    X(11)


    printf("Built on %s %s\n", __DATE__, __TIME__);
    printf("ORDER = %d\n", ORDER);
    printf("COEFF_{MIN,MAX} = {%d, %d}\n", COEFF_MIN, COEFF_MAX);
    print_square(square);

    binsquare_init();
    map = binsquare_map;


#define PUSH(n) \
    do { \
        (void) memcpy(square_orig_c##n, square, sizeof(square)); \
        binsquare_orig_c##n = binsquare; \
    } while (0)

#define LOOP(n) for (c##n = COEFF_MIN; c##n <= COEFF_MAX; c##n ++)

#define ADD(n, c) binsquare = square_addsub_line(square, n, c, binsquare)

#define POP(n) \
    do { \
        (void) memcpy(square, square_orig_c##n, sizeof(square)); \
        binsquare = binsquare_orig_c##n; \
    } while (0)

#define STA(n) PUSH(n); ADD(n, COEFF_MIN-1); LOOP(n) { ADD(n, 1)
#define END(n) } POP(n)

    STA(0);
        STA(1);
            STA(2);
                STA(3);
                    STA(4);
                        STA(5);
                            STA(6);
                                STA(7);
                                    STA(8);
                                        STA(9);
                                            STA(10);
                                                STA(11);
                                                    const size_t off = binsquare >> 6;
                                                    const uint64_t mask = ((uint64_t) 1) << (binsquare & ((binsquare_t) (64-1)));
                                                    if (!(map[off] & mask)) {
                                                        map[off] |= mask;
                                                        innovative_count++;
                                                        //printf("inn: "); print_square(square);
                                                        //print_binsquare(binsquare);
                                                    }
                                                END(11);
                                            END(10);
                                        END(9);
                                    END(8);
                                END(7);
                            END(6);
                        END(5);
                    END(4);
                END(3);
            END(2);
        END(1);
    END(0);


    printf("Final square:    ");
    print_square(square);
    printf("Final binsquare: ");
    print_binsquare(binsquare);
    printf("innovate_count = %" PRIu32 "\n", innovative_count);

    binsquare_finalize();
    return 0;
}
