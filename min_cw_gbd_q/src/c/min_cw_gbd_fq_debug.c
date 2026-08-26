/*
 * DEBUG VERSION: min_cw_gbd_fq_debug.c
 * 
 * Added extensive debugging to find why the algorithm returns no solution.
 * Focus: GF(2) [3,2] case should return minimum weight 2.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>  // Added for debugging

typedef uint8_t field_elem_t;

typedef struct {
    int q;
    field_elem_t* add_table;
    field_elem_t* mul_table;  
    field_elem_t* neg_table;
} field_t;

typedef struct {
    field_elem_t* data;
    int length;
} vector_fq_t;

typedef struct hash_entry {
    field_elem_t* key;
    field_elem_t* message;
    int message_length;
    struct hash_entry* next;
} hash_entry_t;

typedef struct {
    hash_entry_t** buckets;
    int size;
    int key_length;
} hash_table_t;

// Field operations
static inline field_elem_t field_add(const field_t* F, field_elem_t a, field_elem_t b) {
    return F->add_table[a * F->q + b];
}

static inline field_elem_t field_neg(const field_t* F, field_elem_t a) {
    return F->neg_table[a];
}

static inline field_elem_t field_mul(const field_t* F, field_elem_t a, field_elem_t b) {
    return F->mul_table[a * F->q + b];
}

// Vector operations
static vector_fq_t* vector_create(int length) {
    vector_fq_t* v = malloc(sizeof(vector_fq_t));
    if (!v) return NULL;
    
    v->data = calloc(length, sizeof(field_elem_t));
    v->length = length;
    
    if (!v->data) {
        free(v);
        return NULL;
    }
    return v;
}

static void vector_free(vector_fq_t* v) {
    if (v) {
        free(v->data);
        free(v);
    }
}

static void vector_zero(vector_fq_t* v) {
    memset(v->data, 0, v->length * sizeof(field_elem_t));
}

static void vector_copy(vector_fq_t* dest, const vector_fq_t* src) {
    memcpy(dest->data, src->data, src->length * sizeof(field_elem_t));
}

static int vector_hamming_weight(const vector_fq_t* v) {
    int weight = 0;
    for (int i = 0; i < v->length; i++) {
        if (v->data[i] != 0) weight++;
    }
    return weight;
}

// Hash table operations  
static uint32_t hash_function(const field_elem_t* key, int key_length, int field_size) {
    uint32_t hash = 5381;
    for (int i = 0; i < key_length; i++) {
        hash = ((hash << 5) + hash) + key[i];
    }
    return hash;
}

static hash_table_t* hash_table_create(int size, int key_length) {
    hash_table_t* table = malloc(sizeof(hash_table_t));
    if (!table) return NULL;
    
    table->buckets = calloc(size, sizeof(hash_entry_t*));
    table->size = size;
    table->key_length = key_length;
    
    if (!table->buckets) {
        free(table);
        return NULL;
    }
    return table;
}

static void hash_table_free(hash_table_t* table) {
    if (!table) return;
    
    for (int i = 0; i < table->size; i++) {
        hash_entry_t* entry = table->buckets[i];
        while (entry) {
            hash_entry_t* next = entry->next;
            free(entry->key);
            free(entry->message);
            free(entry);
            entry = next;
        }
    }
    free(table->buckets);
    free(table);
}

static int hash_table_insert(hash_table_t* table, const field_elem_t* key, 
                           const field_elem_t* message, int message_length) {
    uint32_t hash = hash_function(key, table->key_length, 0);
    int bucket = hash % table->size;
    
    hash_entry_t* entry = malloc(sizeof(hash_entry_t));
    if (!entry) return -1;
    
    entry->key = malloc(table->key_length * sizeof(field_elem_t));
    entry->message = malloc(message_length * sizeof(field_elem_t));
    
    if (!entry->key || !entry->message) {
        free(entry->key);
        free(entry->message);
        free(entry);
        return -1;
    }
    
    memcpy(entry->key, key, table->key_length * sizeof(field_elem_t));
    memcpy(entry->message, message, message_length * sizeof(field_elem_t));
    entry->message_length = message_length;
    
    entry->next = table->buckets[bucket];
    table->buckets[bucket] = entry;
    
    return 0;
}

static hash_entry_t* hash_table_find(const hash_table_t* table, const field_elem_t* key) {
    uint32_t hash = hash_function(key, table->key_length, 0);
    int bucket = hash % table->size;
    
    hash_entry_t* entry = table->buckets[bucket];
    while (entry) {
        if (memcmp(entry->key, key, table->key_length * sizeof(field_elem_t)) == 0) {
            return entry;
        }
        entry = entry->next;
    }
    return NULL;
}

// Q-ary enumeration
static int qary_increment(field_elem_t* digits, int length, int q) {
    for (int i = 0; i < length; i++) {
        digits[i]++;
        if (digits[i] < q) {
            return 1; // Success
        }
        digits[i] = 0; // Carry
    }
    return 0; // Overflow
}

// Vector reconstruction  
static void reconstruct_vector(const vector_fq_t** generator_rows, int k, 
                              const field_elem_t* message, vector_fq_t* result, 
                              const field_t* F) {
    vector_zero(result);
    for (int i = 0; i < k; i++) {
        if (message[i] != 0) {
            for (int j = 0; j < result->length; j++) {
                field_elem_t term = field_mul(F, message[i], generator_rows[i]->data[j]);
                result->data[j] = field_add(F, result->data[j], term);
            }
        }
    }
}

// Extract positions
static void vector_extract_positions(const vector_fq_t* v, const int* positions, 
                                   int count, field_elem_t* output) {
    for (int i = 0; i < count; i++) {
        output[i] = v->data[positions[i]];
    }
}

// DEBUG version of the main algorithm
int gbd_minimum_weight_fq_debug(
    const field_t* F,
    const vector_fq_t** generator_rows,
    int k, int n,
    const int* filter_positions,
    int s,
    int max_weight,
    vector_fq_t* best_codeword,
    int* best_weight
) {
    printf("DEBUG: Starting GBD F_%d, k=%d, n=%d, s=%d\n", F->q, k, n, s);
    
    int q = F->q;
    int k1 = k / 2;
    int k2 = k - k1;
    
    printf("DEBUG: Split k1=%d, k2=%d\n", k1, k2);
    
    // Simple hash size
    int hash_size = 1024;
    hash_table_t* L1 = hash_table_create(hash_size, s);
    if (!L1) return -1;
    
    vector_fq_t* v1 = vector_create(n);
    vector_fq_t* v2 = vector_create(n);
    vector_fq_t* candidate = vector_create(n);
    field_elem_t* key1 = malloc(s * sizeof(field_elem_t));
    field_elem_t* key2 = malloc(s * sizeof(field_elem_t));
    field_elem_t* message1 = malloc(k1 * sizeof(field_elem_t));
    field_elem_t* message2 = malloc(k2 * sizeof(field_elem_t));
    
    if (!v1 || !v2 || !candidate || !key1 || !key2 || !message1 || !message2) {
        printf("DEBUG: Memory allocation failed\n");
        // Cleanup...
        return -1;
    }
    
    *best_weight = n + 1;
    int l1_entries = 0;
    
    // Phase 1: Build L1 
    printf("DEBUG: Phase 1 - Building L1 table\n");
    memset(message1, 0, k1 * sizeof(field_elem_t));
    
    do {
        // Check if zero message
        int is_zero = 1;
        for (int i = 0; i < k1; i++) {
            if (message1[i] != 0) {
                is_zero = 0;
                break;
            }
        }
        
        if (!is_zero) {
            // Reconstruct vector
            reconstruct_vector(generator_rows, k1, message1, v1, F);
            
            // Extract key
            vector_extract_positions(v1, filter_positions, s, key1);
            
            // Insert into hash table
            hash_table_insert(L1, key1, message1, k1);
            l1_entries++;
            
            if (l1_entries <= 5) { // Debug first few entries
                printf("DEBUG: L1[%d]: message=(", l1_entries);
                for (int i = 0; i < k1; i++) printf("%d,", message1[i]);
                printf(") -> key=(");
                for (int i = 0; i < s; i++) printf("%d,", key1[i]);
                printf(")\n");
            }
        }
        
    } while (qary_increment(message1, k1, q));
    
    printf("DEBUG: L1 built with %d entries\n", l1_entries);
    
    // Phase 2: Search L2
    printf("DEBUG: Phase 2 - Searching L2\n");
    memset(message2, 0, k2 * sizeof(field_elem_t));
    int l2_checks = 0;
    int collisions_found = 0;
    
    do {
        // Check if zero message
        int is_zero = 1;
        for (int i = 0; i < k2; i++) {
            if (message2[i] != 0) {
                is_zero = 0;
                break;
            }
        }
        
        if (!is_zero) {
            l2_checks++;
            
            // Reconstruct vector
            reconstruct_vector(generator_rows + k1, k2, message2, v2, F);
            
            // Extract negated key
            vector_extract_positions(v2, filter_positions, s, key2);
            for (int i = 0; i < s; i++) {
                key2[i] = field_neg(F, key2[i]);
            }
            
            // Look for collision
            hash_entry_t* collision = hash_table_find(L1, key2);
            if (collision) {
                collisions_found++;
                
                // Reconstruct candidate
                reconstruct_vector(generator_rows, k1, collision->message, candidate, F);
                
                // Add v2
                for (int j = 0; j < n; j++) {
                    candidate->data[j] = field_add(F, candidate->data[j], v2->data[j]);
                }
                
                int weight = vector_hamming_weight(candidate);
                
                printf("DEBUG: Collision %d: weight=%d\n", collisions_found, weight);
                
                if (weight > 0 && weight < *best_weight) {
                    *best_weight = weight;
                    vector_copy(best_codeword, candidate);
                    
                    printf("DEBUG: New best weight: %d\n", weight);
                    
                    if (weight <= max_weight) break;
                }
            }
        }
        
    } while (qary_increment(message2, k2, q));
    
    printf("DEBUG: L2 searched %d entries, found %d collisions\n", l2_checks, collisions_found);
    printf("DEBUG: Final best weight: %d\n", *best_weight);
    
    // Cleanup
    hash_table_free(L1);
    vector_free(v1); vector_free(v2); vector_free(candidate);
    free(key1); free(key2); free(message1); free(message2);
    
    return (*best_weight <= n) ? 0 : -1;
}

// Simple wrapper that tries different filter sets automatically
int gbd_adaptive_search_fq(
    const field_t* F,
    const vector_fq_t** generator_rows,
    int k, int n,
    vector_fq_t* best_codeword,
    int* best_weight
) {
    printf("DEBUG: gbd_adaptive_search_fq called\n");
    
    *best_weight = n + 1;
    int found_any = 0;
    
    // Try different filter configurations
    for (int s = 1; s <= k/2 && s <= n/2 && s <= 3; s++) {
        printf("DEBUG: Trying filter size s=%d\n", s);
        
        // Try different starting positions
        for (int start = 0; start <= n-s; start++) {
            int* filter = malloc(s * sizeof(int));
            if (!filter) continue;
            
            for (int i = 0; i < s; i++) {
                filter[i] = (start + i) % n;  // Wrap around
            }
            
            printf("DEBUG: Filter positions: ");
            for (int i = 0; i < s; i++) printf("%d ", filter[i]);
            printf("\n");
            
            vector_fq_t* candidate = vector_create(n);
            int candidate_weight;
            
            if (candidate && gbd_minimum_weight_fq_debug(F, generator_rows, k, n, 
                                                       filter, s, *best_weight - 1,
                                                       candidate, &candidate_weight) == 0) {
                if (candidate_weight < *best_weight) {
                    *best_weight = candidate_weight;
                    vector_copy(best_codeword, candidate);
                    found_any = 1;
                    printf("DEBUG: New best from filter s=%d, start=%d: weight=%d\n", s, start, candidate_weight);
                    
                    // Early exit if found optimal
                    if (candidate_weight <= 2) {
                        vector_free(candidate);
                        free(filter);
                        return 0;
                    }
                }
            }
            
            vector_free(candidate);
            free(filter);
        }
    }
    
    printf("DEBUG: Final adaptive result: weight=%d, found=%d\n", *best_weight, found_any);
    return found_any ? 0 : -1;
}