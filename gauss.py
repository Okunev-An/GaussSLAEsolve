
""" Решение системы линейных алгебраических
уравнений методом Гаусса"""
import numpy as np
import copy

def forward(M, ext=False):
    """Прямой проход метода Гаусса.
     Поочередно избавляемся от переменных"""
    for i in range(0, min(vars, strings)):
        if i > len(M):
            continue
        if not M[i][i]:
            ind = i
            for x in range(i + 1, len(M)):
                if M[x][i]:
                    ind = x
            if ind == i:
                continue
            M[i], M[x] = M[x], M[i].copy()
        for j in range(i + 1, strings):
            if any(abs(M[j])>0.00001):
                M[j] -= M[i] * (M[j][i] / M[i][i])
    return M

def remove_null(M):
    """Удаляем нулевые строки из матрицы"""
    to_del = []
    for i in range(strings):
        if all(abs(M[i]) < 0.0001):
            to_del.append(i)
    M = np.delete(M, to_del, axis=0)
    return M

M = []
"""Программа получает на вход 2 числа : кол-во N
строк и количество M переменных.
Последующие N строк содержат M+1 число """
strings, vars = map(int, input().split())
for i in range(strings):
    M.append([float(x) for x in input().split()])
Mext = np.array(M)
# Расширенная матрица из коэффициентов
M = copy.deepcopy(Mext)[:, :vars]
# Матрица только из левых частей СЛАУ


"""Проводим прямой проход метода Гаусса в обеих 
матрицах, удаляем нулевые строки"""
ext = remove_null(forward(Mext))
simp = remove_null(forward(M))
if len(simp) != len(ext):
    """Система не имеет решений, если 
    получились матрицы разного ранга"""
    print('NO')
elif len(ext) < vars:
    """Если ранг одинаковый, но меньше числа переменных -
    - система неопределена, бесконечное множ. решений"""
    print('INF')
else:
    """В противном случае - есть ровно одно решение, 
    которое мы и будем находить"""
    print('YES')
    answ = []
    """Обратный проход метода Гаусса.
    Проходим по всем строкам получившейся диагональной 
    матрицы(расширенной) в обратном порядке. Вычисляем в i-той 
    строке переменную, домножаем i-ю ячейку
    всех строк на получившееся значение. Список переменных выводим в 
    обратном порядке
    """
    for i in range(len(ext) - 1, -1, -1):
        x = (ext[i][-1] - sum(ext[i][i + 1: -1])) / ext[i][i]
        answ.append(x)
        for j in ext:
            j[i] *= x
    answ = answ[::-1]
    print(*answ)

