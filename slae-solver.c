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

int pivoting(SLAE *slae)
{
    const size_t row_len = (slae->size + 1);
    int was_skip = 0;

    for (size_t i = 0; i < slae->size; i++)
    {
        ld max_koef = slae->coeficients[i * row_len + i];
        size_t max_pos = i;

        for (size_t j = i; j < slae->size; j++)
        {
            if (!float_equals(slae->coeficients[j * row_len + i], 0.0, DBL_EPSILON) &&
                (slae->coeficients[j * row_len + i] > max_koef || float_equals(max_koef, 0.0, DBL_EPSILON)))
            {
                max_koef = slae->coeficients[j * row_len + i];
                max_pos = j;
            }
        }

        if (float_equals(slae->coeficients[max_pos * row_len + i], 0.0, DBL_EPSILON))
        {
            was_skip = 1;
            continue;
        }
        if (max_pos != i)
        {
            ld *tmp = (ld *)malloc(sizeof(ld) * row_len);
            if (tmp == NULL)
            {
                printf("%s\n", "Trouble with allocation of memory in fn 'pivoting'");
                return -1;
            }
            memcpy(tmp, slae->coeficients + max_pos * row_len, row_len * sizeof(ld));
            memcpy(slae->coeficients + max_pos * row_len, slae->coeficients + i * row_len, row_len * sizeof(ld));
            memcpy(slae->coeficients + i * row_len, tmp, row_len * sizeof(ld));
            free(tmp);
        }

        ld scalar = 1.0 / slae->coeficients[i * row_len + i];
        vector_times_scalar(slae->coeficients + i * row_len, row_len, scalar);
        for (size_t j = 0; j < slae->size; j++)
        {
            if (j == i)
            {
                continue;
            }
            else
            {
                ld koefficient = slae->coeficients[j * row_len + i];
                ld *tmp = (ld *)malloc(sizeof(ld) * row_len);
                if (tmp == NULL)
                {
                    printf("%s\n", "Trouble with allocation of memory in fn 'pivoting'");
                    return -1;
                }
                memcpy(tmp, slae->coeficients + i * row_len, row_len * sizeof(ld));
                vector_times_scalar(tmp, row_len, koefficient);

                for (size_t k = 0; k < row_len; k++)
                {
                    if (float_equals((slae->coeficients + j * row_len)[k], tmp[k], DBL_EPSILON))
                    {
                        (slae->coeficients + j * row_len)[k] = 0.0;
                    }
                    else
                    {
                        (slae->coeficients + j * row_len)[k] -= tmp[k];
                    }
                }

                free(tmp);
            }
        }
    }
    return was_skip;
}

answer_type check_solution(SLAE *slae)
{
    const size_t row_len = (slae->size + 1);

    for (size_t i = 0; i < slae->size; i++)
    {
        ld result = 0.0;
        for (size_t j = 0; j < slae->size; j++)
        {
            result += slae->coeficients[i * row_len + j] * slae->coeficients[j * row_len + slae->size];
        }
        if (!float_equals(result, slae->coeficients[i * row_len + slae->size], DBL_EPSILON))
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