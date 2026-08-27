# Быстрый тест алгоритмического сравнения
import time
import itertools
from sage.all import *

print("TEST АЛГОРИТМИЧЕСКОЕ СРАВНЕНИЕ (быстрый тест)")
print("=" * 50)

def brute_force_algorithm(G):
    """Чистый брутфорс алгоритм."""
    F = G.base_ring()
    k = G.nrows()

    best_weight = G.ncols() + 1
    best_vector = None

    count = 0
    for message in itertools.product(F, repeat=k):
        if all(x == 0 for x in message):
            continue
        count += 1
        codeword = sum(message[i] * G[i] for i in range(k))
        weight = codeword.hamming_weight()

        if weight < best_weight:
            best_weight = weight
            best_vector = codeword

        # Ограничиваем для больших полей
        if count > 1000:
            break

    return best_vector, count

def gbd_algorithm_basic(G):
    """Базовый GBD алгоритм без оптимизаций."""
    F = G.base_ring()
    n = G.ncols()
    k = G.nrows()

    # Простое разделение пополам
    k1 = k // 2
    k2 = k - k1
    G1 = G.matrix_from_rows(range(k1))
    G2 = G.matrix_from_rows(range(k1, k))

    # Размер фильтра = половина позиций
    s = min(k1, n // 2)
    filter_positions = list(range(s))

    best_weight = n + 1
    best_vector = None
    operations = 0

    # Строим таблицу L1
    L1 = {}
    for m1 in itertools.product(F, repeat=k1):
        if all(x == 0 for x in m1):
            continue
        operations += 1
        v1 = sum(m1[i] * G1[i] for i in range(k1))
        key = tuple(v1[j] for j in filter_positions)
        L1[key] = v1

        if operations > 500:  # Ограничиваем
            break

    # Ищем коллизии в L2
    for m2 in itertools.product(F, repeat=k2):
        if all(x == 0 for x in m2):
            continue
        operations += 1
        v2 = sum(m2[i] * G2[i] for i in range(k2))
        neg_key = tuple(-v2[j] for j in filter_positions)

        if neg_key in L1:
            candidate = L1[neg_key] + v2
            weight = candidate.hamming_weight()
            if weight < best_weight:
                best_weight = weight
                best_vector = candidate

        if operations > 1000:  # Ограничиваем
            break

    return best_vector, operations

# Тестируем на малых кодах
test_codes = [
    ("GF(2)", 6, 3),
    ("GF(3)", 6, 3),
    ("GF(4)", 6, 3),
    ("GF(5)", 6, 3),
]

for field_name, n, k in test_codes:
    print(f"\n--- {field_name} [{n},{k}] ---")

    # Создаём код
    set_random_seed(42)
    F = eval(field_name)
    G = random_matrix(F, k, n, algorithm='echelonizable', rank=k)
    C = LinearCode(G)

    # Эталон: Sage GAP
    start = time.time()
    sage_dist = C.minimum_distance()
    sage_time = time.time() - start

    # Брутфорс алгоритм
    start = time.time()
    bf_vec, bf_ops = brute_force_algorithm(G)
    bf_time = time.time() - start
    bf_weight = bf_vec.hamming_weight() if bf_vec else None

    # GBD алгоритм
    start = time.time()
    gbd_vec, gbd_ops = gbd_algorithm_basic(G)
    gbd_time = time.time() - start
    gbd_weight = gbd_vec.hamming_weight() if gbd_vec else None

    print(f"Sage GAP:    {sage_time:.4f}s, distance={sage_dist}")
    print(f"Brute Force: {bf_time:.4f}s, weight={bf_weight}, ops={bf_ops}")
    print(f"GBD Basic:   {gbd_time:.4f}s, weight={gbd_weight}, ops={gbd_ops}")

    # Проверка корректности
    correct_bf = "OK" if bf_weight == sage_dist else "FAIL"
    correct_gbd = "OK" if gbd_weight == sage_dist else "FAIL"
    print(f"Корректность: BF={correct_bf}, GBD={correct_gbd}")

    # Сравнение эффективности
    if bf_time > 0 and gbd_time > 0:
        speedup = float(bf_time) / float(gbd_time)
        ops_ratio = float(gbd_ops) / float(bf_ops) if bf_ops > 0 else 0
        print(f"GBD speedup: {speedup:.1f}x (time), ops ratio: {ops_ratio:.1f}")

print(f"\n{'='*50}")
print("ВЫВОДЫ:")
print("- Сравниваем алгоритмы, а не реализации")
print("- GBD должен быть быстрее на средних/больших кодах")
print("- Для малых кодов brute force может выигрывать")