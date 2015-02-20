#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <float.h>
#include "defines.h"

#define F temp
#define MEM 04

double roundd (double);
void swap_ptr (double **, double **);


// Work with memory

int get_mem (symplex *sym, char map) {

    if (map & A) {
        sym->a = (double *) calloc(sym->m * sym->n, sizeof(double));
        if (sym->a != NULL)
            sym->mem |= A;
        else
            return NOT_ENOUGH_MEMORY;
    }
    if (map & B) {
        sym->b = (double *) calloc(sym->m, sizeof(double));
        if (sym->b != NULL)
            sym->mem |= B;
        else
            return NOT_ENOUGH_MEMORY;
    }
    if (map & C) {
        sym->c = (double *) calloc(sym->n, sizeof(double));
        if (sym->c != NULL)
            sym->mem |= C;
        else
            return NOT_ENOUGH_MEMORY;
    }
    if (map & D) {
        sym->d = (double *) calloc(sym->n, sizeof(double));
        if (sym->d != NULL)
            sym->mem |= D;
        else
            return NOT_ENOUGH_MEMORY;
    }
    if (map & X) {
        sym->x = (double *) calloc(sym->n, sizeof(double));
        if (sym->x != NULL)
            sym->mem |= X;
        else
            return NOT_ENOUGH_MEMORY;
    }
    return SUCCESS;
}


void clean_mem (symplex *sym) {

    if (sym->mem & A)
        free(sym->a);
    if (sym->mem & B)
        free(sym->b);
    if (sym->mem & C)
        free(sym->c);
    if (sym->mem & D)
        free(sym->d);
    if (sym->mem & X)
        free(sym->x);
    sym->mem = 0;
}


void copy_symplex (symplex *src, symplex *dest, int limits, char map) {
    int i, j;

    if (limits == 0)
        limits = src->n;
    for(i=0; i < src->m; i++){
        if (map & A)
            for(j=0; j < limits; j++)
                *(dest->a + i*dest->n + j) = *(src->a + i*src->n + j);
        if (map & B)
            *(dest->b + i) = *(src->b + i);
        if (map & BASE)
            dest->base[i] = src->base[i];
    }
    if (map & SIGNS) {
        for (i=0; i < src->m; i++)
            dest->signs[i] = src->signs[i];
        dest->signs[DIMENSION-1] = src->signs[DIMENSION-1];
    }
    if (map & C)
        for(i=0; i < limits; i++)
            *(dest->c + i) = *(src->c + i);
    if (map & D)
        for(i=0; i < limits; i++)
            *(dest->d + i) = *(src->d + i);
    if (map & X)
        for(i=0; i < limits; i++)
            *(dest->x + i) = *(src->x + i);
}



// FUNCTIONS *******************************************************************

int get_first_plan(symplex *obj){
    unsigned int i;
    symplex sym;

    sym.mem = 0;
    sym.m = obj->m;
    sym.n = obj->n + obj->m;
    if (get_mem(&sym, A | C | D | X) == SUCCESS) {
        copy_symplex(obj, &sym, 0, A);
        for (i=0; i<sym.m; i++) {
            *(sym.a + obj->n * (i+1) + i * (obj->m + 1)) = 1;
            *(sym.c + obj->n + i) = -1;
            *(sym.x + obj->n + i) = *(obj->b + i);
            sym.base[i] = obj->n + i;
        }
        get_deltas(&sym);
        switch (i = iteration(&sym)) {
            case SUCCESS:
                copy_symplex(&sym, obj, obj->n, 071);
                clean_mem(&sym);
                return SUCCESS;
            case ITERATION_FAIL:
                clean_mem(&sym);
                return ITERATION_FAIL;
            default:
                clean_mem(&sym);
				// TODO: log this case to stderr
                return i;
		}
    }
    return NOT_ENOUGH_MEMORY;
}



void get_deltas (symplex *sym) {
    unsigned int i, j;

    for (i=0; i < sym->n; i++){
        for(j=0; j < sym->m; j++)
            *(sym->d + i) += *(sym->c + sym->base[j]) * *(sym->a + j*sym->n + i);
        *(sym->d + i) -= *(sym->c + i);
    }
}


void make_standard(symplex *sym){
    unsigned int i,j;

    for (i = 0; i < sym->m; i++)
        if (sym->signs[i] & BIGGER) {
            for (j = 0; j < sym->n; j++)
                 *(sym->a + i * sym->n + j) *= -1;
            *(sym->b + i) *= -1;
            sym->signs[i] = LOWER;
        }
    if (sym->signs[DIMENSION-1] & MIN) {
        for (j = 0; j < sym->n; j++)
            *(sym->c + j) *= -1;
        sym->signs[DIMENSION-1] = MAX;
    }
}


int make_canonic(symplex *sym){
    unsigned int i, j;
    symplex temp;

    for (i = j = 0; i < sym->m; i++)
        if (sym->signs[i] & LOWER)
            j++;

    if (j != 0) {
        temp.m = sym->m;
        temp.n = sym->n + j;
        temp.mem = 0;
        if (get_mem(&temp, A | B | C | D | X )){
            printf("\nError: Can't get enough memory\n");
            clean_mem(&temp);
            return NOT_ENOUGH_MEMORY;
        } else
            copy_symplex(sym, &temp, 0, A | B | C | SIGNS);
        for (i = 0; i < temp.m; i++) {
            if (temp.signs[i] & LOWER) {
                *(temp.a + temp.n * (i+1) - j) = 1;
                *(temp.c + temp.n - j) = 0;
                j--;
                temp.signs[i] = EQUAL;
            }
        }

        clean_mem(sym);
        *sym = temp;

	for (i=0; i < sym->m; i++)
	    if (*(sym->b + i) < 0) {
            *(sym->b + i) *= -1;
            for (j=0; j < sym->n; j++)
                *(sym->a + i*sym->n + j) *= -1;
	    }
    }

    return SUCCESS;
}


// This function implement SYMPLEX METHOD
int iteration(symplex *mx) {
    char state;
    unsigned int i, j, r, s, n, m;
    double min, temp;
    short_symplex one, two, buffer;
	// We use <one> for analysing current state and <two> is a place where we will store our computations, then we'll swap pointers so two-->one and one-->two

    one.a = mx->a;
    one.d = mx->d;
    m = mx->m;
    n = mx->n;
    state = 0;

    do {
        min=0;
        state &= MEM;															//Оставляем флаг выделения памяти, если он присутствует
        for (j=0; j<n; j++)                                                     //Первая стадия итерации
            if (*(one.d+j) < 0) {                                               //Ищем отрицательные оценки замещения
                if ((state & ~MEM) == 0)                                        //Устанавливаем флаг наличия отрицательных оценок замещения
                    state |= ITERATION_FAIL;
                for (i=0; i<m; i++)
                    if (*(one.a+i*n+j) > 0){                                    //Считаем минимум если коэф. замещения неотрицателен
                        temp = *(mx->x + mx->base[i]) / *(one.a + i*n + j);     //Текущий епсилон
                        if (!(state & ITERATION_NEXT) || temp < min) {          //Первый епсилон удовлетворяющий условию
                            min=temp;                                           //Переопределение ведущего элемента
                            r=i;                                                //Выводимый из базиса вектор
                            s=j;                                                //Вводимый в базис вектор
                            state &= ~ITERATION_FAIL;                           //Флаг необходимости следующей итерации
                            state |= ITERATION_NEXT;
                        }
                    }
            }

        F = find_value(mx->x);
        print_main(one.a, mx->x, one.d, m, n, mx->base, NULL, r, s, F, (state & ITERATION_NEXT ? PRINTBASE | PRINTELEMENT : PRINTBASE));
        if (state & ITERATION_NEXT)
            putchar('\n');
        putchar('\n');

        if (state & ITERATION_NEXT) {                                           //Вторая стадия итерации если выполняется условие 3
            if (!(state & MEM)) {                                               //Выделяем память для второй стадии
                two.a = buffer.a = (double *) calloc(m*n, sizeof(double));
                two.d = buffer.d = (double *) calloc(n, sizeof(double));
                if (two.a && two.d == NULL)
                    return NOT_ENOUGH_MEMORY;
                state |= MEM;													//В конце необходимо будет очистить память для временных вычислений
            }
            for (i=0; i<m; i++)                                                 //Вычисление новых коэффициентов замещения
                for (j=0; j<n; j++)
                    if (i == r)
                        *(two.a+i*n+j) = *(one.a+r*n+j) / *(one.a+r*n+s);
                    else
                        *(two.a+i*n+j) = *(one.a+i*n+j) - *(one.a+r*n+j) * (*(one.a+i*n+s) / *(one.a+r*n+s));

            for (j=0; j<n; j++)                                                 //Вычисление новых оценок замещения
                *(two.d+j) = roundd(*(one.d+j) - *(one.d+s) * *(two.a+r*n+j));

            *(mx->x + mx->base[r])=0;                                           //Вычисляем новый вектор X
            mx->base[r] = s;
            for (i=0; i<m; i++)
                if (i == r)
                    *(mx->x + mx->base[r]) = min;
                else
                    *(mx->x + mx->base[i]) = *(mx->x + mx->base[i]) - *(one.a + i*n + s) * min;

            //Делаем вычисленный блок исходным блоком для следующей итерации, для этого меняем указатели.
            swap_ptr(&one.a, &two.a);
            swap_ptr(&one.d, &two.d);

        }
    } while (state & ITERATION_NEXT);

    if (state & MEM) {
        //Копируем результат в исходный блок памяти
        if (two.a == mx->a) {
            for (i=0; i < n*m; i++)
                *(mx->a+i) = *(one.a+i);
            for (i=0; i<n; i++)
                *(mx->d+i) = *(one.d+i);
        }
        //Освобождаем память для временных результатов
        free(buffer.a);
        free(buffer.d);
        state &= ~MEM;
    }

    return state;		// {SUCCESS || ITERATION_FAIL}
}


double roundd (double d) {
    return (fabs(d) < DBL_EPSILON ? 0 : d);
}

int is_integer (double d) {
	return d-floor(d) < FLT_EPSILON;
}

void swap_ptr (double **a, double **b){
    double *temp;

    temp = *a;
    *a = *b;
    *b = temp;
}
