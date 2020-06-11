#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define ld long double

/* utils */

// function to print squere matrix with its free coeficients
void print_squere_matrix(ld **matrix, int matrix_size, FILE *file)
{
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size + 1; j++)
        {
            if (j == matrix_size)
            {
                fprintf(file, "%Lf\n", matrix[i][j]);
            }
            else
            {
                fprintf(file, "%Lf ", matrix[i][j]);
            }
        }
    }
}

/* some structs for slae */

typedef struct SLAE
{
    int size;
    ld **coeficients;
} SLAE;

typedef enum answer_type
{
    no_solutions,
    one_solution,
    inf_solution
} answer_type;

typedef struct SLAE_solution
{
    int SLAE_size;
    answer_type at;
    ld *solution;

} SLAE_solution;

int make_SLAE(int system_size, SLAE *system)
{
    system->coeficients = (ld **)malloc(sizeof(ld *) * system_size);
    for (int i = 0; i < system_size; i++)
    {
        system->coeficients[i] = (ld *)malloc(sizeof(ld) * (system_size + 1));
    }

    system->size = system_size;

    return 0;
}

int make_SLAE_solution(int size, SLAE_solution *slae_solution)
{
    slae_solution = (SLAE_solution *)malloc(sizeof(SLAE_solution));
    slae_solution->SLAE_size = size;
    slae_solution->solution = (ld *)malloc(sizeof(ld) * size);

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

/* slae solver */

void matrix_change_row(ld **matrix, const ld *new_row, const int pos)
{
    ld *to_erase = matrix[pos];
    matrix[pos] = new_row;
    free(to_erase);
}

int max_elt_pos(const SLAE *slae, const int column)
{
    ld max = slae->coeficients[column][column];
    int max_pos = column;
    for (int i = column; i < slae->size; i++)
    {
        if (slae->coeficients[i][column] > max)
        {
            max = slae->coeficients[i][column];
            max_pos = i;
        }
    }
    return max_pos;
}

void row_swap(const ld **arr, const int first_row_index, const int second_row_index)
{
    ld *tmp = arr[first_row_index];
    arr[first_row_index] = arr[second_row_index];
    arr[second_row_index] = tmp;
}

ld *vector_times_scalar(const ld *vector, const int vector_size, const ld scalar)
{
    ld *result = (ld *)malloc(sizeof(ld) * vector_size);
    for (int i = 0; i < vector_size; i++)
    {
        result[i] = vector[i] * scalar;
    }
    return result;
}

ld *vector_subtruction(const ld *vector1, const ld *vector2, const int vector_size)
{
    ld *result = (ld *)malloc(sizeof(ld) * vector_size);
    for (int i = 0; i < vector_size; i++)
    {
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

void eliminate_column(const ld **matrix, const int matrix_size, const int column, const int pivot_index)
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
    Enshures that there are no zeros
        in main diagonal of matrix.
    We will choose for pivot max by absolute value
        element in column. 
 */
void pivoting(SLAE *slae)
{
    for (int i = 0; i < slae->size; i++)
    {
        int max_pos = max_elt_pos(slae, i);
        if (max_pos != i)
        {
            row_swap(slae->coeficients, max_pos, i);
        }
        eliminate_column(slae->coeficients, slae->size, i, i);
    }
}

// SLAE_solution *solve_SLAE(SLAE *slae)
// {

// }

void test_pivoting(SLAE *slae)
{
    FILE *test = fopen("test.txt", "w");
    pivoting(slae);
    print_squere_matrix(slae->coeficients, slae->size, test);
    fclose(test);
}

/* main */

int main(int argc, char **argv)
{

    if (argc != 3)
    {
        printf("%s%d%s", "wrong number of arguments passed to program, got ", argc - 1, ", expected 2");
        exit(1);
    }

    SLAE system;
    if (read_input(argv[1], &system) == 1)
    {
        printf("%s", "error while opening input file");
        exit(1);
    };

    test_pivoting(&system);

    return 0;
}