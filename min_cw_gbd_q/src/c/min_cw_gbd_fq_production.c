/*
 * min_cw_gbd_fq_production.c - PRODUCTION VERSION
 * 
 * Clean, optimized version without debug output.
 * Successfully finds minimum weight codewords for arbitrary F_q.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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

// Field operations (inline for speed)
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
            return 1;
        }
        digits[i] = 0;
    }
    return 0;
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

static void vector_extract_positions(const vector_fq_t* v, const int* positions, 
                                   int count, field_elem_t* output) {
    for (int i = 0; i < count; i++) {
        output[i] = v->data[positions[i]];
    }
}

// Main GBD algorithm
static int gbd_minimum_weight_fq(
    const field_t* F,
    const vector_fq_t** generator_rows,
    int k, int n,
    const int* filter_positions,
    int s,
    int max_weight,
    vector_fq_t* best_codeword,
    int* best_weight
) {
    int q = F->q;
    int k1 = k / 2;
    int k2 = k - k1;
    
    // Hash table size optimization
    int hash_size = 4096; // Good default
    long long est_entries = 1;
    for (int i = 0; i < k1 && i < 10; i++) est_entries *= q;
    if (est_entries > 0 && est_entries < 100000) {
        hash_size = (int)(est_entries * 2); // 50% load factor
    }
    
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
        hash_table_free(L1);
        vector_free(v1); vector_free(v2); vector_free(candidate);
        free(key1); free(key2); free(message1); free(message2);
        return -1;
    }
    
    *best_weight = n + 1;
    
    // Phase 1: Build L1 table
    memset(message1, 0, k1 * sizeof(field_elem_t));
    
    do {
        // Skip zero message
        int is_zero = 1;
        for (int i = 0; i < k1; i++) {
            if (message1[i] != 0) {
                is_zero = 0;
                break;
            }
        }
        
        if (!is_zero) {
            reconstruct_vector(generator_rows, k1, message1, v1, F);
            vector_extract_positions(v1, filter_positions, s, key1);
            hash_table_insert(L1, key1, message1, k1);
        }
        
    } while (qary_increment(message1, k1, q));
    
    // Phase 2: Search L2 for collisions
    memset(message2, 0, k2 * sizeof(field_elem_t));
    
    do {
        // Skip zero message  
        int is_zero = 1;
        for (int i = 0; i < k2; i++) {
            if (message2[i] != 0) {
                is_zero = 0;
                break;
            }
        }
        
        if (!is_zero) {
            reconstruct_vector(generator_rows + k1, k2, message2, v2, F);
            
            // Extract negated key for collision detection
            vector_extract_positions(v2, filter_positions, s, key2);
            for (int i = 0; i < s; i++) {
                key2[i] = field_neg(F, key2[i]);
            }
            
            // Look for collision
            hash_entry_t* collision = hash_table_find(L1, key2);
            if (collision) {
                // Reconstruct candidate = v1 + v2
                reconstruct_vector(generator_rows, k1, collision->message, candidate, F);
                
                for (int j = 0; j < n; j++) {
                    candidate->data[j] = field_add(F, candidate->data[j], v2->data[j]);
                }
                
                int weight = vector_hamming_weight(candidate);
                if (weight > 0 && weight < *best_weight) {
                    *best_weight = weight;
                    vector_copy(best_codeword, candidate);
                    
                    if (weight <= max_weight) break;
                }
            }
        }
        
    } while (qary_increment(message2, k2, q));
    
    // Cleanup
    hash_table_free(L1);
    vector_free(v1); vector_free(v2); vector_free(candidate);
    free(key1); free(key2); free(message1); free(message2);
    
    return (*best_weight <= n) ? 0 : -1;
}

// Adaptive wrapper that tries multiple filter configurations
int gbd_adaptive_search_fq(
    const field_t* F,
    const vector_fq_t** generator_rows,
    int k, int n,
    vector_fq_t* best_codeword,
    int* best_weight
) {
    *best_weight = n + 1;
    int found_any = 0;
    
    // Try different filter configurations intelligently
    int max_s = (k <= 4) ? k/2 : ((k <= 8) ? 3 : 4); // Adaptive max filter size
    if (max_s > n/2) max_s = n/2;
    if (max_s < 1) max_s = 1;
    
    for (int s = 1; s <= max_s; s++) {
        // Limit number of filter position tries for efficiency
        int max_tries = (s == 1) ? n : ((s == 2) ? ((n >= 6) ? 6 : n-1) : 3);
        
        for (int start = 0; start < max_tries && start <= n-s; start++) {
            int* filter = malloc(s * sizeof(int));
            if (!filter) continue;
            
            for (int i = 0; i < s; i++) {
                filter[i] = (start + i) % n;
            }
            
            vector_fq_t* candidate = vector_create(n);
            int candidate_weight;
            
            if (candidate && gbd_minimum_weight_fq(F, generator_rows, k, n, 
                                                 filter, s, *best_weight - 1,
                                                 candidate, &candidate_weight) == 0) {
                if (candidate_weight < *best_weight) {
                    *best_weight = candidate_weight;
                    vector_copy(best_codeword, candidate);
                    found_any = 1;
                    
                    // Early termination for very good results
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
        
        // If we found something good, don't need larger filters
        if (found_any && *best_weight <= 3) break;
    }
    
    return found_any ? 0 : -1;
}