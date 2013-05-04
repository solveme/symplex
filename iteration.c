#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "defines.h"

#define F temp

double roundd (double);
void swap_ptr (double **, double **);

int iteration(symplex *mx) {
    enum { C2 = 01, C3 = 02, MEM = 04 };
    char state;
    unsigned int i, j, r, s, n, m;
    double min, temp;
    short_symplex one, two, buffer;

    one.a = mx->a;
    one.d = mx->d;
    m = mx->m;
    n = mx->n;
    state = 0;

    do {
        min=0;
        state &= MEM;
        for (j=0; j<n; j++)                                                     //Первая стадия итерации
            if (*(one.d+j) < 0) {                                               //Ищем отрицательные оценки замещения
                if ((state & ~MEM) == 0)                                        //Флаг наличия отрицательных оценок замещения
                    state |= C2;
                for (i=0; i<m; i++)
                    if (*(one.a+i*n+j) > 0){                                    //Считаем минимум если коэф. замещения неотрицателен
                        temp = *(mx->x + mx->base[i]) / *(one.a + i*n + j);     //Текущий епсилон
                        if (!(state & C3) || temp < min) {                      //Первый епсилон удовлетворяющий условию
                            min=temp;                                           //Переопределение ведущего элемента
                            r=i;                                                //Выводимый из базиса вектор
                            s=j;                                                //Вводимый в базис вектор
                            state &= ~C2;                                       //Флаг необходимости следующей итерации
                            state |= C3;
                        }
                    }
            }

        F = find_value(mx->x, mx->c, mx->n);
        print_main(one.a, mx->x, one.d, m, n, mx->base, NULL, r, s, F, (state & C3 ? PRINTBASE | PRINTELEMENT : PRINTBASE));
        if (state & C3)
            putchar('\n');
        putchar('\n');

        if (state & C3) {                                                       //Вторая стадия итерации если выполняется условие 3
            if (!(state & MEM)) {                                               //Выделяем память для второй стадии
                two.a = buffer.a = (double *) calloc(m*n, sizeof(double));
                two.d = buffer.d = (double *) calloc(n, sizeof(double));
                if (two.a && two.d == NULL)
                    return 2;
                state |= MEM;
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
    } while (state & C3);

    if (state & MEM) {
        if (two.a == mx->a) {                                                   //Копируем результат в исходный блок памяти
            for (i=0; i < n*m; i++)
                *(mx->a+i) = *(one.a+i);
            for (i=0; i<n; i++)
                *(mx->d+i) = *(one.d+i);
        }
        free(buffer.a);                                                         //Освобождаем память для временных результатов
        free(buffer.d);
        state &= ~MEM;
    }

    return state;
}


double roundd (double d) {
    return (fabs(d) < DBL_EPSILON ? 0 : d);
}


void swap_ptr (double **a, double **b){
    double *temp;

    temp = *a;
    *a = *b;
    *b = temp;
}
