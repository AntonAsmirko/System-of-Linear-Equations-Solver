#!/usr/bin/env python3

import numpy as np
import sys
import random

# command line args description:
# first two args - shape of a matrix
# second two args - a range of values
# fifth arg is number of tests we want to make


def write_test(test, f):
    f.write(str(len(test[0])))
    f.write(' ')
    f.write(str(len(test[0][0])))
    f.write('\n')

    for i in test[0]:
        for j in i:
            f.write("{:.7f}".format(j))
            f.write(' ')
        f.write('\n')

    for i in test[1]:
        f.write("{:.7f}".format(i[0]))
        f.write(' ')

    f.write('\n')


def make_matrix(x_axes, y_axes, lower_bound, upper_bound):

    return (np.random.uniform(low=lower_bound, high=upper_bound, size=(x_axes, y_axes)), np.random.uniform(low=lower_bound, high=upper_bound, size=(x_axes, 1)))


def make_test(matrix):

    return (np.concatenate((matrix[0], matrix[1]), axis=1), np.linalg.solve(matrix[0], matrix[1]))


def random_matrix(x_axes, y_axes, lower_bound, upper_bound):

    return ([[random.uniform(lower_bound, upper_bound) for _ in range(0, y_axes)] for __ in range(0, x_axes)])


def make_linear_dependent_matrix(x_axes, y_axes, lower_bound, upper_bound):
    probably_linear_independent_matrix = random_matrix(
        x_axes, y_axes, lower_bound, upper_bound)

    rand_row = random.randint(0, x_axes - 1)
    second_rand_row = rand_row

    while(second_rand_row == rand_row):
        second_rand_row = random.randint(0, x_axes - 1)

    random_coef = random.random()
    copy_rand_row = probably_linear_independent_matrix[rand_row]

    copy_rand_row = [i * random_coef for i in copy_rand_row]

    print(rand_row)
    print(second_rand_row)
    print(random_coef)
    print()

    probably_linear_independent_matrix[second_rand_row] = copy_rand_row

    return probably_linear_independent_matrix

# Tests modes description:
# pli - probably linear independent
# ld - linear dependent


def main():
    if len(sys.argv) != 7:
        print('wrong number of arguments passed to script')
    print(sys.argv)

    f = open('/home/antonasmirko/MyProjects/C/SLAE-solver/tests/tests.txt', 'w')

    args = [int(sys.argv[2]), int(sys.argv[3]), float(
        sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])]

    mode = sys.argv[1]

    matrix_type = None

    if mode == 'pli':
        matrix_type = make_matrix
    elif mode == 'ld':
        matrix_type = make_linear_dependent_matrix

    if matrix_type == None:
        print('wrong mode\n')
        raise ValueError('wrong mode passed to test generator!')

    for i in range(args[4]):

        matrix = matrix_type(
            args[0], args[1], args[2], args[3])
        if mode == 'pli':
            test = make_test(matrix)
            write_test(test, f)
        elif mode == 'ld':
            test = (matrix, [])
            write_test(test, f)
    f.close()


if __name__ == '__main__':
    main()
