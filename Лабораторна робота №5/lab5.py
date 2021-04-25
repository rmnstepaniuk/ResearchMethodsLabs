import math, random, numpy as np
from scipy.stats import f
from sklearn import linear_model as lm

def average(list):
    average = 0
    for element in list:
        average += element / len(list)
    return average

def dispersion(list):
    list_average = average(list)
    dispersion = 0
    for element in list:
        dispersion += (element - list_average)**2 / len(list)
    return dispersion

def coef_b(x, y):
    for i in range(len(x)):
        x[i].insert(0, 1)
    skm = lm.LinearRegression(fit_intercept = False)
    skm.fit(x, y)
    b = skm.coef_
    return b

def cochrane_criteria():
    global k, N

    gp_denominator = 0
    for disp in dispersion_list:
        gp_denominator += disp
    gp = max(dispersion_list) / gp_denominator
    
    f1 = k - 1
    f2 = N
    gt = 0.3346

    if gp < gt: return True
    else: return False

def students_criteria(b):
    global k, N
    sb = average(dispersion_list)
    
    s_beta_2 = sb / (N * k)
    s_beta = math.sqrt(s_beta_2)

    beta = [sum(average_list[j] * plan_matrix[j][s] for j in range(N)) / N for s in range(m)]

    t = [abs(beta[i]) / s_beta for i in range(m)]

    f3 = (k - 1) * N
    tt = 2.042

    student_check = {}
    for i in range(m):
        if (t[i] > tt): student_check[i] = b[i]
        else: 
            student_check[i] = 0
            b[i] = 0

    return student_check

def fisher_criteria():
    global k, N
    d = 0
    for key in students_criteria:
        if students_criteria[key] != 0: d += 1
    f1 = k - 1
    f2 = N
    f3 = (k - 1) * N
    f4 = N - d

    s2_ad = sum((regression_equation[i] - average_list[i])**2 for i in range(N))

    if (f4 == 0): s2_ad *= k / 10**-12
    else: s2_ad *=  k / f4

    s2_b = average(dispersion_list)

    fp = s2_ad / s2_b

    if fp > f.ppf(q=0.95, dfn=f4, dfd=f3):
        print("\nРівняння регресії неадекватно оригіналу при рівні значимості 0.05")
        return True
    else:
        print("\nРівняння регресії адекватно оригіналу при рівні значимості 0.05")
        return False

x_min = [-1, -10, -8]
x_max = [6, 5, 2]
y_min = 200 + average(x_min)
y_max = 200 + average(x_max)

k = 3
N = 15
m = 10
l = 1.215

# Матриця планування ОЦКП для k = 3
plan_matrix = [
    [-1,    -1,     -1],
    [-1,    -1,      1],
    [-1,     1,     -1],
    [-1,     1,      1],
    [ 1,    -1,     -1],
    [ 1,    -1,      1],
    [ 1,     1,     -1],
    [ 1,     1,      1],
    [-l,     0,      0],
    [ l,     0,      0],
    [ 0,    -l,      0],
    [ 0,     l,      0],
    [ 0,     0,     -l],
    [ 0,     0,      l],
    [ 0,     0,      0]
]

# Матриця планування ОЦКП із натуралізованими значеннями для k = 3
x0 = [(x_min[i] + x_max[i]) / 2 for i in range(k)]
delta_x = [x_max[i] - x0[i] for i in range(k)]
print(f"x\u2080: {x0}\nΔx: {delta_x}")

natur_matrix = [[plan_matrix[i][j] * delta_x[j] + x0[j] for j in range(k)] for i in range(N)]

for i in range(N):
    plan_matrix[i].append(plan_matrix[i][0] * plan_matrix[i][1])
    plan_matrix[i].append(plan_matrix[i][0] * plan_matrix[i][2])
    plan_matrix[i].append(plan_matrix[i][1] * plan_matrix[i][2])
    plan_matrix[i].append(plan_matrix[i][0] * plan_matrix[i][1] * plan_matrix[i][2])
    plan_matrix[i].append(plan_matrix[i][0] ** 2)
    plan_matrix[i].append(plan_matrix[i][1] ** 2)
    plan_matrix[i].append(plan_matrix[i][2] ** 2)

for i in range(N):
    natur_matrix[i].append(natur_matrix[i][0] * natur_matrix[i][1])
    natur_matrix[i].append(natur_matrix[i][0] * natur_matrix[i][2])
    natur_matrix[i].append(natur_matrix[i][1] * natur_matrix[i][2])
    natur_matrix[i].append(natur_matrix[i][0] * natur_matrix[i][1] * natur_matrix[i][2])
    natur_matrix[i].append(natur_matrix[i][0] ** 2)
    natur_matrix[i].append(natur_matrix[i][1] ** 2)
    natur_matrix[i].append(natur_matrix[i][2] ** 2)

print("\nМатриця планування ОЦКП із натуралізованими значеннями:")
for line in natur_matrix:
    print(line)

y_list = [[random.randint(int(y_min), int(y_max)) for _ in range(k)] for __ in range(N)]

print("\nФункції відгуку:")
for line in y_list:
    print(line)

dispersion_list = [dispersion(y_list[i]) for i in range(N)]

# Знайдемо середні значення функцій відгуку за рядками
average_list = [average(y_list[i]) for i in range(N)]

# Перевірка однорідності дисперсії за критерієм Кохрена
cochrane_criteria = cochrane_criteria()
if cochrane_criteria: print("\nДисперсія однорідна")
else: 
    print("\nДисперсія неоднорідна")
    exit()

# Знайдемо коефіцієнти рівняння регресії методом найменших квадратів

b0 =    sum(average_list) / N
b = coef_b(natur_matrix, average_list)

print("\nКоефіцієнти рівняння регресії:")
for i in range(len(b)):
    print(f"b{i} = {round(b[i], 3)}")

# Перевірка значущості коефіцієнтів за критерієм Стьюдента
significant_coefficients = 0
students_criteria = students_criteria(b)

for key in students_criteria:
    if students_criteria[key] != 0:
        significant_coefficients += 1

# Перевірка адекватності за критерієм Фішера
regression_equation = [
            b0 + 
            b[0] * natur_matrix[i][0] + 
            b[1] * natur_matrix[i][1] + 
            b[2] * natur_matrix[i][2] + 
            b[3] * natur_matrix[i][3] + 
            b[4] * natur_matrix[i][4] + 
            b[5] * natur_matrix[i][5] + 
            b[6] * natur_matrix[i][6] +
            b[7] * natur_matrix[i][7] +
            b[8] * natur_matrix[i][8] +
            b[9] * natur_matrix[i][9] for i in range(N)]

fisher_criteria()
