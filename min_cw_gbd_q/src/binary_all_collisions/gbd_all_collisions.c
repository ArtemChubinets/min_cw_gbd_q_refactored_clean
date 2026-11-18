/* Standalone fixed-S binary GBD with full buckets. */
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Each row has words little-endian uint64 limbs, higher padding zero.
 * Each key is a packed projection using key_bits low bits.
 * Returns 1 for a nonzero minimum, 0 for the zero subcode, -1 on error.
 * Counts nonzero candidate checks, and scans the entire right list. */
int gbd_all_collisions(const uint64_t *left, const uint32_t *left_keys, unsigned k1,
                       const uint64_t *right, const uint32_t *right_keys, unsigned k2,
                       unsigned key_bits, size_t words, uint64_t *best,
                       uint64_t *weight_out, uint64_t *scanned, uint64_t *checked) {
    if (!words || words > SIZE_MAX / 30 / sizeof(uint64_t) ||
        k1 > 30 || k2 > 30 || key_bits > 24 ||
        !best || !weight_out || !scanned || !checked ||
        (k1 && (!left || !left_keys)) || (k2 && (!right || !right_keys))) return -1;
    uint32_t keys = UINT32_C(1) << key_bits;
    uint32_t nl = UINT32_C(1) << k1, nr = UINT32_C(1) << k2;
    for (unsigned j = 0; j < k1; ++j) if (left_keys[j] >= keys) return -1;
    for (unsigned j = 0; j < k2; ++j) if (right_keys[j] >= keys) return -1;

    uint32_t *head = calloc(keys, sizeof(*head));
    uint32_t *next = malloc((size_t)nl * sizeof(*next));
    uint64_t *rv = calloc(words, sizeof(*rv));
    uint64_t *candidate = malloc(words * sizeof(*candidate));
    if (!head || !next || !rv || !candidate) {
        free(head); free(next); free(rv); free(candidate);
        return -1;
    }
    *scanned = 0;
    *checked = 0;
    *weight_out = 0;
    memset(best, 0, words * sizeof(*best));

    uint32_t key = 0;
    for (uint32_t i = 0; i < nl; ++i) {
        if (i) key ^= left_keys[__builtin_ctz(i)];
        next[i] = head[key];
        head[key] = i + 1; /* Zero is the empty-chain sentinel. */
    }

    key = 0;
    uint64_t best_weight = UINT64_MAX;
    for (uint32_t i = 0; i < nr; ++i) {
        if (i) {
            unsigned bit = __builtin_ctz(i);
            key ^= right_keys[bit];
            const uint64_t *row = right + bit * words;
            for (size_t j = 0; j < words; ++j) rv[j] ^= row[j];
        }
        for (uint32_t node = head[key]; node; node = next[node - 1]) {
            uint32_t index = node - 1;
            uint32_t mask = index ^ (index >> 1);
            memcpy(candidate, rv, words * sizeof(*rv));
            while (mask) {
                const uint64_t *row = left + __builtin_ctz(mask) * words;
                mask &= mask - 1;
                for (size_t j = 0; j < words; ++j) candidate[j] ^= row[j];
            }
            uint64_t w = 0;
            for (size_t j = 0; j < words; ++j) w += __builtin_popcountll(candidate[j]);
            if (!w) continue;
            ++*checked;
            if (w < best_weight) {
                best_weight = w;
                memcpy(best, candidate, words * sizeof(*best));
            }
        }
    }
    *scanned = nr;
    free(head); free(next); free(rv); free(candidate);
    if (best_weight == UINT64_MAX) return 0;
    *weight_out = best_weight;
    return 1;
}
