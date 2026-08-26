# Отладка C-библиотеки шаг за шагом
import ctypes
import numpy as np
from sage.all import *

print("🔧 ОТЛАДКА C-БИБЛИОТЕКИ")
print("=" * 30)

# 1. Проверка загрузки
try:
    lib = ctypes.CDLL('./min_cw_gbd_fq.so')
    print("✓ Библиотека загружена")
    
    # 2. Проверка функций
    if hasattr(lib, 'gbd_adaptive_search_fq'):
        print("✓ Функция gbd_adaptive_search_fq найдена")
    else:
        print("✗ Главная функция НЕ найдена")
        
    # 3. Простой тест с минимальным кодом
    print("\n🧪 Тест простейшего случая: GF(2) [3,2]")
    
    F2 = GF(2)
    G = matrix(F2, [[1,0,1], [0,1,1]])
    C = LinearCode(G)
    
    print("Генераторная матрица:")
    print(G)
    print("Все кодовые слова:", [c for c in C])
    print("Веса:", [c.hamming_weight() for c in C])
    print("Sage minimum_distance:", C.minimum_distance())
    
    # 4. Проверим, почему C-код не срабатывает
    print("\n🔍 Отладка причин...")
    
    try:
        from gbd_complete_wrapper import gbd_minimum_weight_c
        result = gbd_minimum_weight_c(G, verbose=True)
        print("C result:", result)
    except Exception as e:
        print("C function error:", e)
        import traceback
        traceback.print_exc()
        
except Exception as e:
    print(f"✗ Ошибка библиотеки: {e}")
    
print("\n📋 ДИАГНОЗ:")
print("• Проверить сигнатуры C-функций")
print("• Убедиться что все структуры правильно определены")  
print("• Добавить детальную отладку в Python wrapper")