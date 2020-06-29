#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#define ld long double
#define DBL_EPSILON 1E-7

/* utils */

int float_equals(ld LeftNumber, ld RightNumber, ld Epsilon)
{
    ld diff = fabs(LeftNumber - RightNumber);
    LeftNumber = fabs(LeftNumber);
    RightNumber = fabs(RightNumber);

    ld Largest = (RightNumber > LeftNumber) ? RightNumber : LeftNumber;

    if (Largest < Epsilon)
    {
        if (diff <= Epsilon)
            return 1;
        return 0;
    }
    else
    {
        if (diff <= Largest * Epsilon)
            return 1;
        return 0;
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
    ld **coeficients;
    answer_type at;
} SLAE;

int make_SLAE(int system_size, SLAE *system)
{
    system->coeficients = (ld **)malloc(sizeof(ld *) * system_size);
    for (int i = 0; i < system_size; i++)
    {
        system->coeficients[i] = (ld *)malloc(sizeof(ld) * (system_size + 1));
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
        return 1;
    }
    int N;
    fscanf(file_in, "%d", &N);
    make_SLAE(N, input_system);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            fscanf(file_in, "%lf", &(input_system->coeficients[i][j]));
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
            fprintf(out, "%lf\n", solution->coeficients[i][solution->size]);
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

void matrix_change_row(ld **matrix, ld *new_row, int pos)
{
    ld *to_erase = matrix[pos];
    matrix[pos] = new_row;
    free(to_erase);
}

size_t max_elt_pos(SLAE *slae, size_t column)
{
    ld max = slae->coeficients[column][column];
    size_t max_pos = column;
    for (size_t i = column; i < slae->size; i++)
    {
        if (!float_equals(slae->coeficients[i][column], 0.0, DBL_EPSILON) && (slae->coeficients[i][column] > max || float_equals(max, 0.0, DBL_EPSILON)))
        {
            max = slae->coeficients[i][column];
            max_pos = i;
        }
    }
    return max_pos;
}

void row_swap(ld **arr, int first_row_index, int second_row_index)
{
    ld *tmp = arr[first_row_index];
    arr[first_row_index] = arr[second_row_index];
    arr[second_row_index] = tmp;
}

ld *vector_times_scalar(ld *vector, size_t vector_size, ld scalar)
{
    ld *result = (ld *)malloc(sizeof(ld) * vector_size);
    for (size_t i = 0; i < vector_size; i++)
    {
        result[i] = vector[i] * scalar;

        if (float_equals(result[i], 0.0, DBL_EPSILON))
        {
            result[i] = 0.0;
        }
        else if (float_equals(result[i], 1.0, DBL_EPSILON))
        {
            result[i] = 1.0;
        }
    }
    return result;
}

ld *vector_subtruction(ld *vector1, ld *vector2, size_t vector_size)
{
    ld *result = (ld *)malloc(sizeof(ld) * vector_size);
    for (size_t i = 0; i < vector_size; i++)
    {
        if (float_equals(vector1[i], vector2[i], DBL_EPSILON))
        {
            result[i] = 0.0;
        }
        else
        {
            result[i] = vector1[i] - vector2[i];
        }
    }
    return result;
}

void eliminate_column(ld **matrix, size_t matrix_size, size_t column, size_t pivot_index)
{
    ld scalar = 1.0 / matrix[pivot_index][column];
    ld *norm_pivot = vector_times_scalar(matrix[pivot_index], matrix_size + 1, scalar);
    matrix_change_row(matrix, norm_pivot, pivot_index);
    for (size_t i = 0; i < matrix_size; i++)
    {
        if (i == pivot_index)
        {
            continue;
        }
        else
        {
            ld koefficient = matrix[i][column];
            ld *to_subtruct = vector_times_scalar(matrix[pivot_index], matrix_size + 1, koefficient);
            ld *ith_row_minus_to_subtruct = vector_subtruction(matrix[i], to_subtruct, matrix_size + 1);
            matrix_change_row(matrix, ith_row_minus_to_subtruct, i);
            free(to_subtruct);
        }
    }
}

/*
    Ensures that there are no zeros
        in main diagonal of matrix.
    We will choose for pivot max by absolute value
        element in column.
 */
int pivoting(SLAE *slae)
{
    int was_skip = 0;
    for (size_t i = 0; i < slae->size; i++)
    {
        size_t max_pos = max_elt_pos(slae, i);
        if (float_equals(slae->coeficients[max_pos][i], 0.0, DBL_EPSILON))
        {
            was_skip = 1;
            continue;
        }
        if (max_pos != i)
        {
            row_swap(slae->coeficients, max_pos, i);
        }
        eliminate_column(slae->coeficients, slae->size, i, i);
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
            result += slae->coeficients[i][j] * slae->coeficients[j][slae->size];
        }
        if (!float_equals(result, slae->coeficients[i][slae->size], DBL_EPSILON))
        {
            return no_solutions;
        }
    }
    return inf_solution;
}

void solve_SLAE(SLAE *slae)
{
    int was_skip = pivoting(slae);
    if (was_skip)
    {
        slae->at = check_solution(slae);
    }
    else
    {
        slae->at = one_solution;
    }
}

/* free utils */

void free_SLE_from_stack(SLAE *sle)
{
    for (int i = 0; i < sle->size; i++)
    {
        free(sle->coeficients[i]);
    }
    free(sle->coeficients);
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
        printf("%s", "error while opening input file\n");
        exit(1);
    };
    solve_SLAE(&system);

    if (write_output(argv[2], &system))
    {
        printf("%s", "an error occured while attempt to write solution to file\n");
        exit(1);
    }

    free_SLE_from_stack(&system);

    return 0;
}