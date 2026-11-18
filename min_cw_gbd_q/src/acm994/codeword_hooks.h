/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef ACM994_CODEWORD_HOOKS_H
#define ACM994_CODEWORD_HOOKS_H

#include <stdint.h>

/* Private context for one call in the single-threaded CLI. */
void acm994_use_columns(int *original_column);
void acm994_swap_columns(int a, int b);
void acm994_consider_word(int weight, int words, const uint32_t *prefix,
                          const uint32_t *suffix);
int acm994_log(const char *format, ...);

#endif
