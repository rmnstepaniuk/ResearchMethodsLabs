import random, math, numpy as np
from scipy.stats import f

def main(m):
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

    def cochrane_criteria():
        global m, N

        gp_denominator = 0
        for disp in dispersion_list:
            gp_denominator += disp
        gp = max(dispersion_list) / gp_denominator
        

        f1 = m - 1
        f2 = N
        gt = 0.7679

        if gp < gt: return True
        else: return False

    def students_criteria(b):
        global m, N
        sb = average(dispersion_list)
        
        s_beta_2 = sb / (N * m)

        s_beta = math.sqrt(s_beta_2)

        beta = [sum(average_list[i] * plan_matrix[j][i] for i in range(N)) for j in range(N)]

        t = [abs(beta[i]) / s_beta for i in range(N)]

        f3 = (m - 1) * N
        tt = 2.306

        student_check = {}
        for i in range(N):
            if (t[i] > tt): student_check[i] = b[i]
            else: 
                student_check[i] = 0
                b[i] = 0

        return student_check

    def fisher_criteria():
        global m, N
        d = 0
        for key in students_criteria:
            if students_criteria[key] != 0: d += 1
        f1 = m - 1
        f2 = N
        f3 = (m - 1) * N
        f4 = N - d

        s2_ad = sum((regression_equation[i] - average_list[i])**2 for i in range(N))

        if (f4 == 0): s2_ad *= m / 10**-12
        else: s2_ad *=  m / f4

        s2_b = average(dispersion_list)

        fp = s2_ad / s2_b

        if fp > f.ppf(q=0.95, dfn=f4, dfd=f3):
            print("\nРівняння регресії неадекватно оригіналу при рівні значимості 0.05")
            return True
        else:
            print("\nРівняння регресії адекватно оригіналу при рівні значимості 0.05")
            return False

    x_min = [-20, -15, -15]
    x_max = [15, 35, -10]
    y_min = 200 + average(x_min)
    y_max = 200 + average(x_max)

    print(f"\nМінімальне значення ф-ції відгуку: {int(y_min)}")
    print(f"Максимальне значення ф-ції відгуку: {int(y_max)}")

    print(f"\nМінімальні значення х: {x_min}")
    print(f"Максимальні значення х: {x_max}")

    plan_matrix = [
                    [1, -1, -1, -1, 1, 1, 1, -1],
                    [1, -1, -1, 1, 1, -1, -1, 1],
                    [1, -1, 1, -1, -1, 1, -1, 1],
                    [1, -1, 1, 1, -1, -1, 1, -1],
                    [1, 1, -1, -1, -1, -1, 1, 1],
                    [1, 1, -1, 1, -1, 1, -1, -1],
                    [1, 1, 1, -1, 1, -1, -1, -1],
                    [1, 1, 1, 1, 1, 1, 1, 1]]

    print(f"\nМатриця планування експерименту:")
    for line in plan_matrix:
        print(line)

    x0 = [plan_matrix[i][0] for i in range(N)]
    x1 = [plan_matrix[i][1] for i in range(N)]
    x2 = [plan_matrix[i][2] for i in range(N)]
    x3 = [plan_matrix[i][3] for i in range(N)]
    x12 = [plan_matrix[i][4] for i in range(N)]
    x13 = [plan_matrix[i][5] for i in range(N)]
    x23 = [plan_matrix[i][6] for i in range(N)]
    x123 = [plan_matrix[i][7] for i in range(N)]

    experiment_matrix = [
                    [-20, -15, -15, -20*(-15), -20*(-15), -15*(-15), -20*(-15)*(-15)],
                    [-20, -15, -10, -20*(-15), -20*(-10), -15*(-10), -20*(-15)*(-10)],
                    [-20, 35, -15, -20*35, -20*(-15), 35*(-15), -20*35*(-15)],
                    [-20, 35, -10, -20*35, -20*(-10), 35*(-10), -20*35*(-10)],
                    [15, -15, -15, 15*(-15), 15*(-15), -15*(-15), 15*(-15)*(-15)],
                    [15, -15, -10, 15*(-15), 15*(-10), -15*(-10), 15*(-15)*(-10)],
                    [15, 35, -15, 15*35, 15*(-15), 35*(-15), 15*35*(-15)],
                    [15, 35, -10, 15*35, 15*(-10), 35*(-10), 15*35*(-10)]]
    print(f"\nЕкспериментальна матриця:")
    for line in experiment_matrix:
        print(line)

    y_list = [[random.randint(int(y_min), int(y_max)) for _ in range(m)] for __ in range(N)]

    print(f"\nФункції відгуку:")
    for line in y_list:
        print(line)

    dispersion_list = [
                        dispersion(y_list[0]),
                        dispersion(y_list[1]),
                        dispersion(y_list[2]),
                        dispersion(y_list[3]),
                        dispersion(y_list[4]),
                        dispersion(y_list[5]),
                        dispersion(y_list[6]),
                        dispersion(y_list[7])]

    # Рівняння регресії
    # y = b0 + b1*x1 + b2*x2 + b3*x3 + b12*x1*x2 + b13*x1*x3 + b23*x2*x3 + b123*x1*x2*x3

    # Знайдемо середні значення функцій відгуку за рядками

    average_list = [average(y_list[0]),
                    average(y_list[1]),
                    average(y_list[2]),
                    average(y_list[3]),
                    average(y_list[4]),
                    average(y_list[5]),
                    average(y_list[6]),
                    average(y_list[7])]


    # Перевірка однорідності дисперсії за критерієм Кохрена
    cochrane_criteria = cochrane_criteria()
    if cochrane_criteria: print("\nДисперсія однорідна")
    else: 
        print("\nДисперсія неоднорідна")
        exit()


    # Знайдемо коефіцієнти рівняння регресії методом найменших квадратів
    b0 =    sum(average_list[i] for i in range(N)) / N
    b1 =    sum(average_list[i] * plan_matrix[i][1] for i in range(N)) / N
    b2 =    sum(average_list[i] * plan_matrix[i][2] for i in range(N)) / N
    b3 =    sum(average_list[i] * plan_matrix[i][3] for i in range(N)) / N
    b12 =   sum(average_list[i] * plan_matrix[i][4] for i in range(N)) / N
    b13 =   sum(average_list[i] * plan_matrix[i][5] for i in range(N)) / N
    b23 =   sum(average_list[i] * plan_matrix[i][6] for i in range(N)) / N
    b123 =  sum(average_list[i] * plan_matrix[i][7] for i in range(N)) / N
    b = [b0, b1, b2, b3, b12, b13, b23, b123]

    # Перевірка значущості коефіцієнтів за критерієм Стьюдента
    significant_coefficients = 0
    students_criteria = students_criteria(b)

    for key in students_criteria:
        if students_criteria[key] != 0:
            significant_coefficients += 1

    # Перевірка адекватності за критерієм Фішера
    regression_equation = [
                b[0] + 
                b[1] * experiment_matrix[i][0] + 
                b[2] * experiment_matrix[i][1] + 
                b[3] * experiment_matrix[i][2] + 
                b[4] * experiment_matrix[i][3] + 
                b[5] * experiment_matrix[i][4] + 
                b[6] * experiment_matrix[i][5] + 
                b[7] * experiment_matrix[i][6] for i in range(N)]

    if fisher_criteria() == True:
        if m <= 20:
            m += 1
            main(m)
            exit()
        else: print("Рівняння регресії неадекватно оригіналу при m <= 20")
            
    else:
        
        print(f"\nСередні значення y: {average_list}")
        
        print(f"\nКоефіцієнти: {b}")

        print(f"\nКількість вагомих коефіцієнтів: {significant_coefficients}")

        for i in range(N):
            if b[i] == 0: print(f"Видалимо з рівняння регресії невагомий коефіцінт - b{i}")

        print(f"\nРівняння регресії:")
        print(f"{round(b[0], 3)} + {round(b[1], 3)}*x1 + {round(b[2], 3)}*x2 + {round(b[3], 3)}*x3 + {round(b[4], 3)}*x1*x2 + {round(b[5], 3)}*x1*x13 + {round(b[6], 3)}*x2*x3 + {round(b[7], 3)}*x1*x2*x3")
m = 3
N = 8

main(m)
