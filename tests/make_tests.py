#!/usr/bin/env python3

import numpy as np
import sys

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
            f.write(str(j))
            f.write(' ')
        f.write('\n')
    for i in test[1]:
        f.write(str(i[0]))
        f.write(' ')
    f.write('\n')


def make_matrix(x_axes, y_axes, lower_bound, upper_bound):
    return (np.random.uniform(low=lower_bound, high=upper_bound, size=(x_axes, y_axes)), np.random.uniform(low=lower_bound, high=upper_bound, size=(x_axes, 1)))


def make_test(matrix):
    return (np.concatenate((matrix[0], matrix[1]), axis=1), np.linalg.solve(matrix[0], matrix[1]))


def main():
    if len(sys.argv) != 6:
        print('wrong number of arguments passed to script')
    print(sys.argv)
    f = open('/home/antonasmirko/MyProjects/C/SLAE-solver/tests/tests.txt', 'w')
    args = [int(sys.argv[1]), int(sys.argv[2]), float(
        sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])]
    for i in range(args[4]):
        matrix = make_matrix(
            args[0], args[1], args[2], args[3])
        test = make_test(matrix)
        write_test(test, f)
    f.close()


if __name__ == '__main__':
    main()
