# Простой тест оптимизированной версии  
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("Тест оптимизированной min_cw_gbd_q")
print("=" * 40)

try:
    from min_cw_gbd_q import min_cw_gbd_q
    print("✓ Импорт успешен")
    
    # Простой бинарный тест  
    F2 = GF(2)
    G = matrix(F2, 3, 6, [
        [1,0,0,1,1,0],
        [0,1,0,1,0,1], 
        [0,0,1,0,1,1]
    ])
    C = LinearCode(G)
    
    print(f"\nТест binary [6,3]:")
    result_vec = min_cw_gbd_q(C)
    weight = result_vec.hamming_weight() if result_vec is not None else None
    distance = C.minimum_distance()
    
    print(f"Наш результат: weight={weight}, vector={result_vec}")
    print(f"Sage результат: distance={distance}") 
    print(f"Совпадение: {'✓' if weight == distance else '✗'}")
    
    # Простой тернарный тест
    F3 = GF(3)
    G3 = matrix(F3, 2, 5, [
        [1,0,1,2,1],
        [0,1,2,1,2]
    ])
    C3 = LinearCode(G3)
    
    print(f"\nТест ternary [5,2]:")
    result_vec3 = min_cw_gbd_q(C3)
    weight3 = result_vec3.hamming_weight() if result_vec3 is not None else None
    distance3 = C3.minimum_distance()
    
    print(f"Наш результат: weight={weight3}, vector={result_vec3}")
    print(f"Sage результат: distance={distance3}")
    print(f"Совпадение: {'✓' if weight3 == distance3 else '✗'}")
    
except Exception as e:
    import traceback
    print(f"✗ Ошибка: {e}")
    print(traceback.format_exc())