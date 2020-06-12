#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#define ld long double

/* utils */

int float_equals(ld a, ld b, ld epsilon)
{
    if (fabs(a - b) < epsilon)
    {
        return 1;
    }
    else
    {
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
            fscanf(file_in, "%Lf", &(input_system->coeficients[i][j]));
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
            fprintf(out, "%Lf\n", solution->coeficients[i][solution->size]);
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

int max_elt_pos(SLAE *slae, int column)
{
    ld max = slae->coeficients[column][column];
    int max_pos = column;
    for (int i = column; i < slae->size; i++)
    {
        if (slae->coeficients[i][column] != 0 && (slae->coeficients[i][column] > max || max == 0))
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

ld *vector_times_scalar(ld *vector, int vector_size, ld scalar)
{
    ld *result = (ld *)malloc(sizeof(ld) * vector_size);
    for (int i = 0; i < vector_size; i++)
    {
        result[i] = vector[i] * scalar;
    }
    return result;
}

ld *vector_subtruction(ld *vector1, ld *vector2, int vector_size)
{
    ld *result = (ld *)malloc(sizeof(ld) * vector_size);
    for (int i = 0; i < vector_size; i++)
    {
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

void eliminate_column(ld **matrix, int matrix_size, int column, int pivot_index)
{
    ld *norm_pivot = vector_times_scalar(matrix[pivot_index], matrix_size + 1, 1.0 / matrix[pivot_index][column]);
    matrix_change_row(matrix, norm_pivot, pivot_index);
    for (int i = 0; i < matrix_size; i++)
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
    for (int i = 0; i < slae->size; i++)
    {
        int max_pos = max_elt_pos(slae, i);
        if (slae->coeficients[max_pos][i] == 0.0)
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
    for (int i = 0; i < slae->size; i++)
    {
        ld result = 0.0;
        for (int j = 0; j < slae->size; j++)
        {
            result += slae->coeficients[i][j] * slae->coeficients[j][slae->size];
        }
        if (!float_equals(result, slae->coeficients[i][slae->size], 0.000001))
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