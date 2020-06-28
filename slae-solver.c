#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define ld float
#define DUMMY_EPS 0.000000000001f
#define ABS_TOL 1.192092896e-07
#define REL_TOL 1.0f
#define FLT_EPSILON_C 1.192092896e-07

/* utils */

void print_squere_matrix(ld **matrix, size_t matrix_size, FILE *file)
{
    for (size_t i = 0; i < matrix_size; i++)
    {
        for (size_t j = 0; j < matrix_size + 1; j++)
        {
            if (j == matrix_size)
            {
                fprintf(file, "%.16f\n", matrix[i][j]);
            }
            else
            {
                fprintf(file, "%.16f ", matrix[i][j]);
            }
        }
    }
}

void matrix_copy(ld **matrix, size_t matrix_size, ld **dest)
{
    for (size_t i = 0; i < matrix_size; i++)
    {
        for (size_t j = 0; j < matrix_size + 1; j++)
        {
            dest[i][j] = matrix[i][j];
        }
    }
}

int float_equals(ld a, ld b, ld abs_tol, ld rel_tol)
{
    ld diff = fabs(a - b);
    return (diff <= fabs(rel_tol * b)) ||
           (diff <= fabs(rel_tol * a)) ||
           (diff <= abs_tol);
}

ld **allocate_2D_double_array(size_t x_shape, size_t y_shape)
{
    ld **result = (ld **)malloc(sizeof(ld *) * x_shape);
    for (size_t i = 0; i < x_shape; i++)
    {
        result[i] = (ld *)malloc(sizeof(ld) * y_shape);
    }

    return result;
}

void append_int_to_query(char *query, size_t *path_len, int val)
{
    char *tmp = (char *)malloc(sizeof(char) * 50);
    sprintf(tmp, "%d", val);
    for (size_t i = 0; i < 50 && tmp[i] != '\0'; i++, (*path_len)++)
    {
        query[*path_len] = tmp[i];
    }
    query[*path_len] = ' ';
    (*path_len)++;
    free(tmp);
}

void append_str_to_query(char *query, size_t *path_len, char *str)
{
    for (size_t j = 0; *path_len < 2000 && str[j] != '\0'; (*path_len)++, j++)
    {
        query[*path_len] = str[j];
    }
    query[*path_len] = ' ';
    (*path_len)++;
}

char *make_query(char *mode, size_t x_shape, size_t y_shape, char *lower_bound, char *upper_bound, size_t count_tests)
{
    char *query = (char *)malloc(sizeof(char) * 2000);
    size_t path_len = 21;
    strcpy(query, "./tests/make_tests.py");
    query[path_len] = ' ';
    path_len++;
    append_str_to_query(query, &path_len, mode);
    append_int_to_query(query, &path_len, x_shape);
    append_int_to_query(query, &path_len, y_shape);
    append_str_to_query(query, &path_len, lower_bound);
    append_str_to_query(query, &path_len, upper_bound);
    append_int_to_query(query, &path_len, count_tests);

    return query;
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
    size_t size;
    ld **coeficients;
    answer_type at;
} SLAE;

typedef struct test_sla
{
    SLAE *system;
    ld *solution;
    ld **original_system;
} test_sla;

int make_SLAE(size_t system_size, SLAE *system)
{

    system->coeficients = allocate_2D_double_array(system_size, system_size + 1);
    system->size = system_size;
    system->at = not_solved_yet;

    return 0;
}

int make_test_sla(test_sla *test, SLAE *system, ld *solution)
{
    test->system = system;
    test->solution = solution;
    test->original_system = allocate_2D_double_array(system->size, system->size + 1);
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
            fscanf(file_in, "%f", &(input_system->coeficients[i][j]));
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
        for (size_t i = 0; i < solution->size; i++)
        {
            fprintf(out, "%.16f\n", solution->coeficients[i][solution->size]);
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

void matrix_change_row(ld **matrix, ld *new_row, size_t pos)
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
        if (slae->coeficients[i][column] != 0 && (slae->coeficients[i][column] > max || max == 0))
        {
            max = slae->coeficients[i][column];
            max_pos = i;
        }
    }
    return max_pos;
}

void row_swap(ld **arr, size_t first_row_index, size_t second_row_index)
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
        if (float_equals(result[i], 0.0f, ABS_TOL, FLT_EPSILON_C))
        {
            result[i] = 0.0f;
        }
        else if (float_equals(result[i], 1.0f, ABS_TOL, FLT_EPSILON_C))
        {
            result[i] = 1.0f;
        }
    }
    return result;
}

ld *vector_subtruction(ld *vector1, ld *vector2, size_t vector_size)
{
    ld *result = (ld *)malloc(sizeof(ld) * vector_size);
    for (size_t i = 0; i < vector_size; i++)
    {
        if (float_equals(vector1[i], vector2[i], ABS_TOL, FLT_EPSILON_C))
        {
            result[i] = 0.0f;
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
    ld scalar = 1.0f / matrix[pivot_index][column];
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
        if (float_equals(slae->coeficients[max_pos][i], 0.0f, ABS_TOL, FLT_EPSILON_C))
        {
            was_skip = 1;
            continue;
        }
        if (max_pos != i)
        {
            row_swap(slae->coeficients, max_pos, i);
        }
        eliminate_column(slae->coeficients, slae->size, i, i);
        print_squere_matrix(slae->coeficients, slae->size, stdout);
    }
    return was_skip;
}

answer_type check_solution(SLAE *slae)
{
    for (size_t i = 0; i < slae->size; i++)
    {
        ld result = 0.0f;
        for (size_t j = 0; j < slae->size; j++)
        {
            result += slae->coeficients[i][j] * slae->coeficients[j][slae->size];
        }
        if (!float_equals(result, slae->coeficients[i][slae->size], ABS_TOL, FLT_EPSILON_C))
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

/* utill for testing pivoting */

void test_pivoting(SLAE *slae)
{
    FILE *test = fopen("test.txt", "w");
    pivoting(slae);
    print_squere_matrix(slae->coeficients, slae->size, test);
    fclose(test);
}

/* free utils */

void free_SLE_from_stack(SLAE *sle)
{
    for (size_t i = 0; i < sle->size; i++)
    {
        free(sle->coeficients[i]);
    }
    free(sle->coeficients);
}

/* read test */

int read_test(FILE *tests, SLAE *system, test_sla *test, int require_solution)
{
    int x_size, y_size;
    fscanf(tests, "%d%d", &x_size, &y_size);
    for (size_t i = 0; i < x_size; i++)
    {
        for (size_t j = 0; j < y_size; j++)
        {
            fscanf(tests, "%f", &(system->coeficients[i][j]));
        }
    }
    ld *solution = (ld *)malloc(sizeof(ld) * x_size);
    if (require_solution)
    {
        for (size_t i = 0; i < x_size; i++)
        {
            fscanf(tests, "%f", &(solution[i]));
        }
    }

    make_test_sla(test, system, solution);
    matrix_copy(system->coeficients, system->size, test->original_system);
    return 0;
}

/* util for writing tests with wrong outcome */

int write_wrong_one_solution(test_sla *test, FILE *out)
{
    fprintf(out, "%s\n", "TEST:");
    print_squere_matrix(test->original_system, test->system->size, out);
    fprintf(out, "%s\n", "EXPECTED:");
    for (size_t i = 0; i < test->system->size; i++)
    {
        fprintf(out, "%.16f ", test->solution[i]);
    }
    fprintf(out, "\n%s\n", "GOT:");
    for (size_t i = 0; i < test->system->size; i++)
    {
        fprintf(out, "%.16f ", test->system->coeficients[i][test->system->size]);
    }
    fprintf(out, "%s", "\n");
    return 0;
}

/* run tests */

int check_correctness_one_solution(test_sla *test)
{
    size_t size_x = test->system->size;
    for (size_t i = 0; i < size_x; i++)
    {
        if (!float_equals(test->system->coeficients[i][size_x], test->solution[i], ABS_TOL, FLT_EPSILON_C))
        {
            return 0;
        }
    }
    return 1;
}

void run_tests(char *mode, size_t x_shape, size_t y_shape, char *lower_bound, char *upper_bound, int count_tests)
{
    char *query = make_query(mode, x_shape, y_shape, lower_bound, upper_bound, count_tests);
    system(query);
    free(query);
    FILE *tests;
    FILE *tests_results;
    int tests_passed = 0;
    if (!(tests = fopen("./tests/tests.txt", "r")))
    {
        printf("%s", "an error occurred while attempt to open file tests.txt");
    }

    if (!(tests_results = fopen("testsresults.txt", "w")))
    {
        printf("%s", "an error occurred while attempt to open file testsresults.txt");
    }

    SLAE *system;
    test_sla *test;

    for (size_t i = 0; i < count_tests; i++)
    {
        SLAE *system = (SLAE *)malloc(sizeof(SLAE));
        make_SLAE(x_shape, system);
        test_sla *test = (test_sla *)malloc(sizeof(test_sla));

        if (strcmp(mode, "ld") == 0)
        {
            read_test(tests, system, test, 0);
        }
        else if (strcmp(mode, "pli") == 0)
        {
            read_test(tests, system, test, 1);
        }

        solve_SLAE(system);

        if (strcmp(mode, "ld") == 0)
        {
            if (test->system->at == inf_solution)
            {
                printf("%s%d%s\n", "Test № ", i + 1, " passed successfully");
                tests_passed++;
            }
            else
            {
                fprintf(tests_results, "%d\n", test->system->size);

                print_squere_matrix(test->original_system, test->system->size, tests_results);
            }
        }
        else if (strcmp(mode, "pli") == 0)
        {
            if (!check_correctness_one_solution(test))
            {
                write_wrong_one_solution(test, tests_results);
            }
            else
            {
                printf("%s%d%s\n", "Test № ", i + 1, " passed successfully");
                printf("%s\n", "YOUR ANSWER:");

                for (size_t j = 0; j < x_shape; j++)
                {
                    printf("%.16f ", system->coeficients[j][system->size]);
                }

                printf("\n%s\n", "SYSTEMS ANSWER:");

                for (size_t j = 0; j < x_shape; j++)
                {
                    printf("%.16f ", test->solution[j]);
                }
                printf("%s", "\n");
                tests_passed++;
            }
        }

        free(system);
        free(test);
    }
    fprintf(tests_results, "%d %d\n", tests_passed, count_tests);
    fclose(tests);
    fclose(tests_results);
}

/* main */

// int main(int argc, char **argv)
// {
//     if (argc != 3)
//     {
//         printf("%s%d%s", "wrong number of arguments passed to program, got ", argc - 1, ", expected 2\n");
//         exit(1);
//     }

//     SLAE system;
//     if (read_input(argv[1], &system) == 1)
//     {
//         printf("%s", "error while opening input file\n");
//         exit(1);
//     };
//     solve_SLAE(&system);

//     printf("\n");

//     print_squere_matrix(system.coeficients, system.size, stdout);

//     if (write_output(argv[2], &system))
//     {
//         printf("%s", "an error occured while attempt to write solution to file\n");
//         exit(1);
//     }

//     test_pivoting(&system);

//     free_SLE_from_stack(&system);

//     return 0;
// }

int my_getnbr(char *str)
{
    int result;
    int puiss;

    result = 0;
    puiss = 1;
    while (('-' == (*str)) || ((*str) == '+'))
    {
        if (*str == '-')
            puiss = puiss * -1;
        str++;
    }
    while ((*str >= '0') && (*str <= '9'))
    {
        result = (result * 10) + ((*str) - '0');
        str++;
    }
    return (result * puiss);
}

/* 
run_tests args description:
mode - ld (linear dependent), pli (pribably linear independent)
x_shape
y_shape
lower_bound
upper_bound
count_tests
*/

int main(int argc, char **argv)
{
    if (argc != 7)
    {
        printf("%s\n", "wrong number of arguments been passed to programm, expected 6 args");
        exit(1);
    }

    run_tests(argv[1], my_getnbr(argv[2]), my_getnbr(argv[3]), argv[4], argv[5], my_getnbr(argv[6]));

    return 0;
}