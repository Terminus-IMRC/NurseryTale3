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

static uint64_t* gather(uint64_t *sum, const char *filename, size_t *sizep)
{
    size_t i;
    FILE *fp;
    struct stat sb;
    uint64_t *p;
    size_t rets;

    if (lstat(filename, &sb)) {
        fprintf(stderr, "lstat: %s: %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (*sizep != 0 && sb.st_size != *sizep) {
        fprintf(stderr, "size is different\n");
        exit(EXIT_FAILURE);
    }

    *sizep = sb.st_size;

    p = calloc(*sizep, 1);
    if (p == NULL) {
        fprintf(stderr, "calloc: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        fprintf(stderr, "fopen: %s: %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    rets = fread(p, 1, *sizep, fp);
    if (rets != *sizep) {
        fprintf(stderr, "fread: %s: %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (fclose(fp)) {
        fprintf(stderr, "close: %s: %s\n", filename, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (sum == NULL)
        return p;

    for (i = 0; i < *sizep / sizeof(*sum); i ++)
        sum[i] |= p[i];

    free(p);

    return NULL;
}

static uint32_t count_innovative(unsigned long long *sum, const size_t size)
{
    size_t i;
    uint32_t count = 0;
    for (i = 0; i < size / sizeof(*sum); i ++)
        count += __builtin_popcountll(sum[i]);
    return count;
}

int main(int argc, char *argv[])
{
    void *sum;
    int i;
    size_t size = 0;
    uint32_t count;

    if (argc <= 1) {
        fprintf(stderr, "error: Specify bsmap files\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "Gathering %d file(s)\n", argc - 1);

    sum = gather(NULL, argv[1], &size);
    for (i = 2; i < argc; i ++)
        gather(sum, argv[i], &size);

    count = count_innovative(sum, size);
    fprintf(stderr, "innovative_count = %" PRIu32 "\n", count);

    if (!isatty(STDOUT_FILENO)) {
        fprintf(stderr, "Writing gathered bsmap to stdout\n");
        size_t rets = fwrite(sum, 1, size, stdout);
        if (rets != size) {
            fprintf(stderr, "fwrite: stdout: %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
    } else
        printf("Redirect stdout to file to output gathered bsmap\n");

    free(sum);

    return 0;
}
