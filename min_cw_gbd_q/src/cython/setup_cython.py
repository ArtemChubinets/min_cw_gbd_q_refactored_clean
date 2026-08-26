"""
Setup script для сборки Cython extension gbd_fast.

Использование:
    python setup_cython.py build_ext --inplace
"""
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import os

# Пути к Sage headers (могут отличаться)
sage_include_dirs = []

# Стандартные пути где может быть Sage
possible_sage_paths = [
    "/usr/lib/python3.14/site-packages/sage", 
    "/usr/local/lib/python*/site-packages/sage",
    "/opt/sagemath/lib/python*/site-packages/sage"
]

# Найдем Sage включения
for path_pattern in possible_sage_paths:
    import glob
    for path in glob.glob(path_pattern):
        if os.path.exists(path):
            sage_include_dirs.append(path)
            break

extensions = [
    Extension(
        "gbd_fast",
        sources=["gbd_fast.pyx"],
        include_dirs=[
            np.get_include(),
        ] + sage_include_dirs,
        language="c++",
        extra_compile_args=[
            "-O3",           # Максимальная оптимизация 
            "-ffast-math",   # Быстрая математика
            "-march=native", # Оптимизация под текущий процессор
            "-std=c++11",    # C++11 стандарт
        ],
        extra_link_args=["-O3"],
    )
]

setup(
    name="gbd_fast",
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "language_level": 3,
            "boundscheck": False,    # Отключить проверки границ
            "wraparound": False,     # Отключить отрицательные индексы
            "cdivision": True,       # C-деление (быстрее)
            "initializedcheck": False, # Не проверять инициализацию
        }
    ),
    zip_safe=False,
)

print("""
Cython extension собран!

Для тестирования:
    from gbd_fast import pack_key_q_fast, linear_combination_fast, gbd_core_optimized
    
Для интеграции в min_cw_gbd_q:
    1. Скопировать gbd_fast.*.so в min_cw_gbd_q/
    2. Добавить import в core.py:
        try:
            from .gbd_fast import gbd_core_optimized
            USE_CYTHON = True
        except ImportError:
            USE_CYTHON = False
    3. Использовать в gbd_search_q при USE_CYTHON=True
""")