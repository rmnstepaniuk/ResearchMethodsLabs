import random

matrix = [[random.randint(0, 20) for i in range(3)] for i in range(8)]

# for i in matrix:
#     print(i)

a0 = random.randint(0, 20)
a1 = random.randint(0, 20)
a2 = random.randint(0, 20)
a3 = random.randint(0, 20)

# print("\n")

Y_list = []
for line in matrix:
    Y = a0 + a1*line[0] + a2*line[1] + a3*line[2]
    Y_list.append(Y)
    # print(f"Y = {Y}")
print(Y_list)

x0_1_set = {matrix[i][0] for i in range(8)}
x0_2_set = {matrix[i][1] for i in range(8)}
x0_3_set = {matrix[i][2] for i in range(8)}

# print(f"\n{x0_1_set}\n{x0_2_set}\n{x0_3_set}\n")

x0_1 = (max(x0_1_set) + min(x0_1_set)) / 2
dx_1 = x0_1 - min(x0_1_set)

x0_2 = (max(x0_2_set) + min(x0_2_set)) / 2
dx_2 = x0_2 - min(x0_2_set)

x0_3 = (max(x0_3_set) + min(x0_3_set)) / 2
dx_3 = x0_3 - min(x0_3_set)

print(f"X0 = {x0_1}\ndx = {dx_1}\n\nX0 = {x0_2}\ndx = {dx_2}\n\nX0 = {x0_3}\ndx = {dx_3}\n")

x0_list = [x0_1, x0_2, x0_3]
dx_list = [dx_1, dx_2, dx_3]

normalization = [[round((matrix[i][j]-x0_list[j])/dx_list[j], 5) for j in range(3)] for i in range(8)]
for line in normalization:
    print(line)

Yet = a0 + a1*x0_1 + a2*x0_2 + a3*x0_3
print(f"\nY = {a0} + {a1}*{x0_1}+ {a2}*{x0_2} + {a3}*{x0_3}\nYэт = {Yet}")