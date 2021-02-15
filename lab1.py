# Методи наукових досліджень 
#
# Степанюк Роман Вікторович ІВ-91 ФІОТ
#
# Варіант 125 ( Yет <-- )

import random

matrix = [[random.randint(0, 20) for i in range(3)] for i in range(8)]
print("Значення факторів у точках експерименту:")
for i in matrix:
    print(i)

a0 = random.randint(0, 20)
a1 = random.randint(0, 20)
a2 = random.randint(0, 20)
a3 = random.randint(0, 20)

Y_list = []
for line in matrix:
    Y = a0 + a1*line[0] + a2*line[1] + a3*line[2]
    Y_list.append(Y)
print(f"\nФункції відгуку у кожній точці експерименту:\n{Y_list}\n")

x0_1_set = {matrix[i][0] for i in range(8)}
x0_2_set = {matrix[i][1] for i in range(8)}
x0_3_set = {matrix[i][2] for i in range(8)}

x0_1 = (max(x0_1_set) + min(x0_1_set)) / 2
dx_1 = x0_1 - min(x0_1_set)

x0_2 = (max(x0_2_set) + min(x0_2_set)) / 2
dx_2 = x0_2 - min(x0_2_set)

x0_3 = (max(x0_3_set) + min(x0_3_set)) / 2
dx_3 = x0_3 - min(x0_3_set)

print(f"Нульовий рівень для першого фактора:\nX0 = {x0_1}\ndx = {dx_1}\n\nНульовий рівень для другого фактора:\nX0 = {x0_2}\ndx = {dx_2}\n\nНульовий рівень для третього фактора:\nX0 = {x0_3}\ndx = {dx_3}\n")

x0_list = [x0_1, x0_2, x0_3]
dx_list = [dx_1, dx_2, dx_3]

normalization = []
print("Значення факторів у точках експерименту після нормалізації:")
for i in range(8):
    normalization.append([])
    for j in range(3):
        normalization[i].append(round(((matrix[i][j] - x0_list[j]) / dx_list[j]), 5))
        if j == 2:
            print(normalization[i])

Yet = a0 + a1*x0_1 + a2*x0_2 + a3*x0_3

print(f"\nФ-ція відгуку від нульових рівнів факторів:\nYет = {Yet}")

diff_list = []  # у цей масив зберігаються значення різниці (Yi - Yет) для подальшого знаходження найменшої додатної з цих різниць
                # це потрібно безпосередньо для виконання завдання варіанту

for Y in Y_list:
    diff_list.append(Y - Yet)
print(f"\nРізниця ф-цій відгуку і ф-ції відгуку від нульових рівнів факторів:\n{diff_list}")


min_d = diff_list[0]
for d in diff_list:
    if d < 0:
        continue
    else:
        if min_d < 0:
            min_d = d
        elif d < abs(min_d):
            min_d = d
print(f"\nЗначення функції відгуку, яке найблище до значення еталонної ф-ції відгуку:\n{min_d + Yet}")

for i in range(len(Y_list)):
    if (min_d + Yet == Y_list[i]):
        print(f"Точка плану, що задовольняє критерій оптимальності (Yет <--):\n{matrix[i]}")    
