/* SPDX-License-Identifier: GPL-3.0-or-later */
#define _POSIX_C_SOURCE 200809L
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codeword_hooks.h"
#include "generations.h"

static int ncols, best_weight, verbose;
static int *columns;
static char *best;

int acm994_log(const char *format, ...) {
    if (!verbose) return 0;
    va_list args;
    va_start(args, format);
    int result = vfprintf(stderr, format, args);
    va_end(args);
    return result;
}

void acm994_use_columns(int *original_column) { columns = original_column; }

void acm994_swap_columns(int a, int b) {
    int saved = columns[a];
    columns[a] = columns[b];
    columns[b] = saved;
}

void acm994_consider_word(int weight, int words, const uint32_t *prefix,
                          const uint32_t *suffix) {
    if (weight <= 0 || weight >= best_weight) return;
    if (!columns || words != (ncols + 31) / 32) abort();
    for (int j = 0; j < ncols; ++j) {
        uint32_t bits = prefix[j / 32] ^ suffix[j / 32];
        best[columns[j]] = (bits >> (31 - j % 32)) & 1u;
    }
    best_weight = weight;
}

static void fail(const char *message) {
    fprintf(stderr, "acm994-codeword: %s\n", message);
    exit(EXIT_FAILURE);
}

static int dimension(char **text) {
    char *end;
    errno = 0;
    long value = strtol(*text, &end, 10);
    if (end == *text || errno || value <= 0 || value > INT_MAX - 32 ||
        (*end && !isspace((unsigned char)*end)))
        fail("expected two positive dimensions: k n");
    *text = end;
    return (int)value;
}

static int row_echelon(char *matrix, int k, int n, int *pivots) {
    int rank = 0;
    for (int j = 0; j < n && rank < k; ++j) {
        int r = rank;
        while (r < k && !matrix[(size_t)r * n + j]) ++r;
        if (r == k) continue;
        char *pivot = matrix + (size_t)rank * n;
        char *row = matrix + (size_t)r * n;
        for (int i = j; i < n; ++i) {
            char t = pivot[i];
            pivot[i] = row[i];
            row[i] = t;
        }
        for (r = rank + 1; r < k; ++r) {
            row = matrix + (size_t)r * n;
            if (row[j])
                for (int i = j; i < n; ++i) row[i] ^= pivot[i];
        }
        pivots[rank++] = j;
    }
    return rank;
}

int main(void) {
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length = getline(&line, &capacity, stdin);
    if (length < 0 || memchr(line, '\0', (size_t)length))
        fail("missing or malformed header");
    char *rest = line;
    int k = dimension(&rest);
    int n = dimension(&rest);
    while (isspace((unsigned char)*rest)) ++rest;
    if (*rest) fail("expected exactly two dimensions: k n");
    if ((int64_t)k * n > INT_MAX) fail("matrix exceeds CALGO 994 index limits");
    char *matrix = malloc((size_t)k * n);
    int *pivots = malloc((size_t)(k < n ? k : n) * sizeof(*pivots));
    best = malloc((size_t)n);
    if (!matrix || !pivots || !best) fail("out of memory");
    for (int r = 0; r < k; ++r) {
        length = getline(&line, &capacity, stdin);
        if (length > 0 && line[length - 1] == '\n') --length;
        if (length > 0 && line[length - 1] == '\r') --length;
        if (length != n) fail("expected one n-bit row per line");
        for (int j = 0; j < n; ++j) {
            if (line[j] != '0' && line[j] != '1') fail("rows must be binary");
            matrix[(size_t)r * n + j] = (char)(line[j] - '0');
        }
    }
    free(line);

    int rank = row_echelon(matrix, k, n, pivots);
    if (!rank) fail("zero code has no nonzero codeword");
    int64_t pairs = (int64_t)rank * (rank - 1) / 2;
    /* The saved-pair buffers are allocated via upstream signed-int counts. */
    if (pairs > (INT_MAX / (int)sizeof(uint32_t)) / ((n + 31) / 32))
        fail("saved-combination matrix exceeds CALGO 994 index limits");
    ncols = n;
    best_weight = n + 1;
    verbose = getenv("ACM994_VERBOSE") != NULL;
    int weight = compute_distance_of_matrix(matrix, rank, n);
    if (weight != best_weight || weight > n || weight <= 0)
        fail("missing or inconsistent search witness");

    /* Check the recovered witness independently in the original coordinates. */
    char *residual = malloc((size_t)n);
    if (!residual) fail("out of memory");
    memcpy(residual, best, (size_t)n);
    for (int r = 0; r < rank; ++r) {
        if (residual[pivots[r]])
            for (int j = pivots[r]; j < n; ++j)
                residual[j] ^= matrix[(size_t)r * n + j];
    }
    int actual_weight = 0;
    for (int j = 0; j < n; ++j) {
        if (residual[j]) fail("search witness is outside the input row span");
        actual_weight += best[j];
    }
    if (actual_weight != weight) fail("search witness weight mismatch");
    printf("RESULT %d ", weight);
    for (int j = 0; j < n; ++j) putchar('0' + best[j]);
    putchar('\n');
    free(residual);
    free(best);
    free(pivots);
    free(matrix);
    return ferror(stdout) ? EXIT_FAILURE : EXIT_SUCCESS;
}
