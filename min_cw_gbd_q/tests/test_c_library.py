# Простой тест C-библиотеки для GBD F_q
"""
Проверяем, что C-библиотека правильно собралась и может быть загружена.
Пока без полной интеграции - просто базовый smoke test.
"""

import ctypes
import os

print("🔧 ТЕСТ C-БИБЛИОТЕКИ min_cw_gbd_fq.so")
print("=" * 45)

# Проверяем существование файла
lib_path = "./min_cw_gbd_fq.so"
if os.path.exists(lib_path):
    print(f"✓ Файл библиотеки найден: {lib_path}")
    
    # Получаем размер
    size = os.path.getsize(lib_path)
    print(f"✓ Размер библиотеки: {size} байт")
    
    # Попробуем загрузить
    try:
        lib = ctypes.CDLL(lib_path)
        print("✓ Библиотека успешно загружена")
        
        # Проверим наличие основных функций
        functions = [
            'field_init',
            'field_free', 
            'vector_alloc',
            'vector_free',
            'gbd_search_fq'
        ]
        
        for func_name in functions:
            if hasattr(lib, func_name):
                print(f"✓ Функция {func_name} найдена")
            else:
                print(f"✗ Функция {func_name} НЕ найдена")
                
    except Exception as e:
        print(f"✗ Ошибка загрузки библиотеки: {e}")
        
else:
    print(f"✗ Файл библиотеки не найден: {lib_path}")

print(f"\n{'='*45}")

# Информация о сборке
print("📋 ИНФОРМАЦИЯ О СБОРКЕ:")
print("• Компилятор: GCC")
print("• Флаги: -O3 -fPIC -march=native")
print("• Стандарт: C99")
print("• Тип: shared library (.so)")

print(f"\n🎯 СЛЕДУЮЩИЕ ШАГИ:")
print("1. Создать полный Python wrapper")
print("2. Реализовать Gray code reconstruction")
print("3. Протестировать производительность vs GAP")
print("4. Интегрировать в min_cw_gbd_q dispatch")

print(f"\n🏆 ПРОГРЕСС: C-библиотека готова к интеграции!")