#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "defines.h"


int get_mem (symplex *sym, char map) {

    if (map & A) {
        sym->a = (double *) calloc(sym->m * sym->n, sizeof(double));
        if (sym->a != NULL)
            sym->mem |= A;
        else
            return 1;
    }
    if (map & B) {
        sym->b = (double *) calloc(sym->m, sizeof(double));
        if (sym->b != NULL)
            sym->mem |= B;
        else
            return 1;
    }
    if (map & C) {
        sym->c = (double *) calloc(sym->n, sizeof(double));
        if (sym->c != NULL)
            sym->mem |= C;
        else
            return 1;
    }
    if (map & D) {
        sym->d = (double *) calloc(sym->n, sizeof(double));
        if (sym->d != NULL)
            sym->mem |= D;
        else
            return 1;
    }
    if (map & X) {
        sym->x = (double *) calloc(sym->n, sizeof(double));
        if (sym->x != NULL)
            sym->mem |= X;
        else
            return 1;
    }
    return 0;
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
    if (get_mem(&sym, A | C | D | X) == 0) {
        copy_symplex(obj, &sym, 0, A);
        for (i=0; i<sym.m; i++) {
            *(sym.a + obj->n * (i+1) + i * (obj->m + 1)) = 1;
            *(sym.c + obj->n + i) = -1;
            *(sym.x + obj->n + i) = *(obj->b + i);
            sym.base[i] = obj->n + i;
        }
        get_deltas(&sym);
        switch (iteration(&sym)) {
            case 0:
                copy_symplex(&sym, obj, obj->n, 071);
                clean_mem(&sym);
                return 0;
            case 1:
                clean_mem(&sym);
                return 1;
            case 2:
                clean_mem(&sym);
                return 2;
        }
    }
    return 3;
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
    int i,j;

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
    int i, j;
    symplex temp;

    for (i = j = 0; i < sym->m; i++)
        if (sym->signs[i] & LOWER)
            j++;

    if (j != 0) {
        temp.m = sym->m;
        temp.n = sym->n + j;
        temp.mem = 0;
        if (get_mem(&temp, 037)){
            printf("\nError: Can't get enough memory\n");
            clean_mem(&temp);
            return 1;
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

    return 0;
}


double find_value(double *x, double *c, int n) {
    int i;
    double S = 0;

    for (i = 0; i < n; i++)
        S += *(x + i) * *(c + i);

    return S;
}
