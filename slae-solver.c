#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#define ld double
#define DBL_EPSILON 1E-6

/* utils */

int float_equals(ld a, ld b, ld eps)
{
    ld diff = fabs(a - b);
    a = fabs(a);
    b = fabs(b);

    ld lar = (b > a) ? b : a;

    if (lar < eps)
    {
        if (diff <= eps)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        if (diff <= lar * eps)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

/* some structs for slae */

typedef enum answer_type
{
    no_solutions,
    one_solution,
    inf_solution,
    not_solved_yet
} answer_type;

typedef struct SLAE
{
    int size;
    ld *coeficients;
    answer_type at;
} SLAE;

int make_SLAE(int system_size, SLAE *system)
{
    system->coeficients = (ld *)malloc(sizeof(ld) * system_size * (system_size + 1));
    if (system->coeficients == NULL)
    {
        printf("%s\n", "An error occured while attempt to allocate memory for matrix of system's coefficients");
        return 1;
    }
    system->size = system_size;
    system->at = not_solved_yet;

    return 0;
}

/* read input */

int read_input(char *in_file, SLAE *input_system)
{
    FILE *file_in;
    if (!(file_in = fopen(in_file, "r")))
    {
        printf("%s", "error while opening input file\n");
        return 1;
    }
    int N;
    fscanf(file_in, "%d", &N);
    if (make_SLAE(N, input_system))
    {
        fclose(file_in);
        return 1;
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            fscanf(file_in, "%lf", &(input_system->coeficients[i * (input_system->size + 1) + j]));
        }
    }

    fclose(file_in);
    return 0;
}

/* write output */

int write_output(char *out_file, SLAE *solution)
{
    FILE *out;
    if (!(out = fopen(out_file, "w")))
    {
        return 1;
    }
    switch (solution->at)
    {
    case no_solutions:
        fprintf(out, "%s", "no solution\n");
        break;
    case inf_solution:
        fprintf(out, "%s", "many solutions\n");
        break;
    case one_solution:
        for (int i = 0; i < solution->size; i++)
        {
            fprintf(out, "%.10lf\n", solution->coeficients[i * (solution->size + 1) + solution->size]);
        }
        break;
    default:
        printf("%s", "strange behavior, troubles with solution cathegory\n");
        fclose(out);
        return 2;
    }
    fclose(out);
    return 0;
}

/* slae solver */

size_t max_elt_pos(SLAE *slae, size_t column)
{
    ld max = slae->coeficients[column * (slae->size + 1) + column];
    size_t max_pos = column;
    for (size_t i = column; i < slae->size; i++)
    {
        if (!float_equals(slae->coeficients[i * (slae->size + 1) + column], 0.0, DBL_EPSILON) &&
            (slae->coeficients[i * (slae->size + 1) + column] > max || float_equals(max, 0.0, DBL_EPSILON)))
        {
            max = slae->coeficients[i * (slae->size + 1) + column];
            max_pos = i;
        }
    }
    return max_pos;
}

void vector_times_scalar(ld *vector, size_t vector_size, ld scalar)
{
    for (size_t i = 0; i < vector_size; i++)
    {
        vector[i] = vector[i] * scalar;

        if (float_equals(vector[i], 0.0, DBL_EPSILON))
        {
            vector[i] = 0.0;
        }
        else if (float_equals(vector[i], 1.0, DBL_EPSILON))
        {
            vector[i] = 1.0;
        }
    }
}

void vector_subtruction(ld *vector1, ld *vector2, size_t vector_size)
{
    for (size_t i = 0; i < vector_size; i++)
    {
        if (float_equals(vector1[i], vector2[i], DBL_EPSILON))
        {
            vector1[i] = 0.0;
        }
        else
        {
            vector1[i] = vector1[i] - vector2[i];
        }
    }
}

int eliminate_column(ld *matrix, size_t matrix_size, size_t column, size_t pivot_index)
{
    const size_t row_len = (matrix_size + 1);
    ld scalar = 1.0 / matrix[pivot_index * row_len + column];
    vector_times_scalar(matrix + pivot_index * row_len, row_len, scalar);
    for (size_t i = 0; i < matrix_size; i++)
    {
        if (i == pivot_index)
        {
            continue;
        }
        else
        {
            ld koefficient = matrix[i * row_len + column];
            ld *tmp = (ld *)malloc(sizeof(ld) * (matrix_size + 1));
            if (tmp == NULL)
            {
                printf("%s\n", "Trouble with allocation of memory in fn 'eliminate_column'");
                return 1;
            }
            memcpy(tmp, matrix + pivot_index * row_len, row_len * sizeof(ld));
            vector_times_scalar(tmp, row_len, koefficient);
            vector_subtruction(matrix + i * row_len, tmp, row_len);
            free(tmp);
        }
    }
    return 0;
}

int pivoting(SLAE *slae)
{
    int was_skip = 0;

    for (size_t i = 0; i < slae->size; i++)
    {
        size_t max_pos = max_elt_pos(slae, i);
        if (float_equals(slae->coeficients[max_pos * (slae->size + 1) + i], 0.0, DBL_EPSILON))
        {
            was_skip = 1;
            continue;
        }
        if (max_pos != i)
        {
            ld *tmp = (ld *)malloc(sizeof(ld) * (slae->size + 1));
            if (tmp == NULL)
            {
                printf("%s\n", "Trouble with allocation of memory in fn 'pivoting'");
                return -1;
            }
            memcpy(tmp, slae->coeficients + max_pos * (slae->size + 1), (slae->size + 1) * sizeof(ld));
            memcpy(slae->coeficients + max_pos * (slae->size + 1), slae->coeficients + i * (slae->size + 1), (slae->size + 1) * sizeof(ld));
            memcpy(slae->coeficients + i * (slae->size + 1), tmp, (slae->size + 1) * sizeof(ld));
            free(tmp);
        }
        if (eliminate_column(slae->coeficients, slae->size, i, i))
        {
            return -1;
        }
    }
    return was_skip;
}

answer_type check_solution(SLAE *slae)
{
    for (size_t i = 0; i < slae->size; i++)
    {
        ld result = 0.0;
        for (size_t j = 0; j < slae->size; j++)
        {
            result += slae->coeficients[i * (slae->size + 1) + j] * slae->coeficients[j * (slae->size + 1) + slae->size];
        }
        if (!float_equals(result, slae->coeficients[i * (slae->size + 1) + slae->size], DBL_EPSILON))
        {
            return no_solutions;
        }
    }
    return inf_solution;
}

int solve_SLAE(SLAE *slae)
{
    int was_skip = pivoting(slae);
    if (was_skip == -1)
    {
        return 1;
    }
    if (was_skip)
    {
        slae->at = check_solution(slae);
    }
    else
    {
        slae->at = one_solution;
    }
    return 0;
}

/* main */

int main(int argc, char **argv)
{

    if (argc != 3)
    {
        printf("%s%d%s", "wrong number of arguments passed to program, got ", argc - 1, ", expected 2\n");
        exit(1);
    }

    SLAE system;
    if (read_input(argv[1], &system) == 1)
    {
        exit(1);
    }

    if (solve_SLAE(&system))
    {
        free(system.coeficients);
        exit(1);
    }

    if (write_output(argv[2], &system))
    {
        printf("%s", "an error occured while attempt to write solution to file\n");
        free(system.coeficients);
        exit(1);
    }

    free(system.coeficients);

    return 0;
}