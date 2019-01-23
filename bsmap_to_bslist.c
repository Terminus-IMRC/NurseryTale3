#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <inttypes.h>

#if !defined(ORDER)
#error "Define ORDER"
#endif
#if   ORDER == 5
typedef uint32_t binsquare_t;
#elif ORDER == 6
typedef uint64_t binsquare_t;
#else
#error "Unsupported ORDER specified"
#endif

static void bsmap_to_bslist(void)
{
    uint64_t u;
    size_t i, j, cnt = 0;
    for (i = 0; ; i ++) {
        const size_t rets = fread(&u, sizeof(u), 1, stdin);
        if (rets != 1) {
            if (feof(stdin))
                break;
            fprintf(stderr, "rets (%zu) != 1\n", rets);
            fprintf(stderr, "fread: %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < sizeof(u) * 8; j ++) {
            const __typeof__(u) m = ((__typeof__(u)) 1) << j;
            /* Exclude the null binsquare. */
            if (i == 0 && j == 0)
                continue;
            if (u & m) {
                const binsquare_t bs = (i * sizeof(u) * 8) | j;
                const size_t rets = fwrite(&bs, sizeof(bs), 1, stdout);
                if (rets != 1) {
                    fprintf(stderr, "fwrite: %s\n", strerror(errno));
                    exit(EXIT_FAILURE);
                }
                cnt ++;
            }
        }
    }
    fprintf(stderr, "%zu entries (%zu bytes) written\n",
            cnt, cnt * sizeof(binsquare_t));
    fprintf(stderr, "Note: This count excludes the null binsquare\n");
}

int main(void)
{
    if (isatty(STDIN_FILENO)) {
        fprintf(stderr, "stdin is a tty!\n");
        exit(EXIT_FAILURE);
    }

    if (isatty(STDOUT_FILENO)) {
        fprintf(stderr, "stdout is a tty!\n");
        exit(EXIT_FAILURE);
    }

    bsmap_to_bslist();

    return 0;
}
