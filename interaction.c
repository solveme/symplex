#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "defines.h"

#define NOLAST 0100
#define FIXED 01

int get_length (double);

void print_number (double, int, int, char);

void print_border (char *, int, char);

void print_symplex(symplex *sym, char map){
    int base[DIMENSION] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    char *s;
    int *bas;
    double *c, *b;

    if (map & PRINTTARGET)
		s = sym->signs;
    else 
		s = NULL;

    if (map & PRINTBASE)
        bas = sym->base;
    else
        bas = base;

    if (map & PRINTLAST) {
        c = sym->d;
        b = sym->x;
    } else {
        c = sym->c;
        b = sym->b;
    }

    print_main (sym->a, b, c, sym->m, sym->n, bas, s, 0, 0, 0, map);
}


void print_main(double *a, double *x, double *d, int m, int n, int base[], char signs[], int r, int s, double F, char map) {
    int length, i, j, t;
    char border[20], *lena;

    length = 0;
    lena = (char *) calloc(m*n + 2*n, sizeof(char));

    for (i=0; i<m; i++) {
        for (j=0; j<n; j++)
            if ((*(lena+ i*n + j) = get_length(*(a + i*n + j))) > length)
                length = *(lena+ i*n +j);
        if ((*(lena+ m*n + i) = get_length(*(x + i))) > length)
            length = *(lena+ m*n + i);
    }
    for (i=0; i<n; i++)
        if ((*(lena+ m*n + n + i) = get_length(*(d + i))) > length)
            length = *(lena+ m*n + n + i);


    if (length < DBL_DIG/3)
        length = DBL_DIG/3;
    else if (length < DBL_DIG)
        length += 2;

    for (i = 0; i < length+3; i++)
        border[i] = '-';
    border[length+3] = '\0';


    print_border(border, n, map);

    if (map & PRINTBASE)
        printf("| BS |%*c ", length+2, map & BINX ? 'B' : 'X');
    putchar('|');
    for (i = 0; i < n; i++)
        printf("%*s%u |", length-(i+1)/10+1, "A", i+1);
    if (map & PRINTTARGET)
        printf("    |");
    if ((map & PRINTBASE) != 1)
        printf("%*c |", length+2, map & BINX ? 'B' : 'X');
    putchar('\n');
    print_border(border, n, map);

    for (i = 0; i < m; i++) {
        if (map & PRINTBASE) {
            printf("| %*u |", 2, base[i]+1);
            print_number(*(x + base[i]), length, *(lena + m*n + i), FIXED);
            putchar(' ');
        }
        putchar('|');
        for (j = 0; j < n; j++) {
            print_number(*(a + i*n + j), length, *(lena + i*n + j), FIXED);
            if ((map & PRINTELEMENT) && (r == i) && (s == j))
                putchar('*');
            else
                putchar(' ');
            putchar('|');
        }
        if (map & PRINTTARGET)
            switch (signs[i]) {
                case LOWER:
                    printf(" <= |");
                    break;
                case BIGGER:
                    printf(" >= |");
                    break;
                default:
                    printf("  = |");
            }
        if ((map & PRINTBASE) != 1) {
            print_number(*(x + base[i]), length, *(lena + m*n + i), FIXED);
            printf(" |");
        }
        putchar('\n');
    }

    print_border(border, n, map);

    if (map & PRINTBASE) {
        printf("| F= |");
        t = get_length(F);
        print_number(F, length, t, FIXED);
        putchar(' ');
    }

    putchar('|');
    for (i = 0; i < n; i++) {
        print_number(*(d + i), length, *(lena + m*n + n + i), FIXED);
        printf(" |");
    }
    if (map & PRINTTARGET)
        printf(" -> |%-*s|", length+3, (signs[DIMENSION-1] & MAX ? " max " : " min "));
    putchar('\n');

    print_border(border, n, map);
    //print_border(border, n, (map & PRINTTARGET ? map : map | NOLAST));
    free(lena);
}


void print_border (char *border, int n, char map){
    int i;

    if (map & PRINTBASE)
        printf("+----");
    putchar('+');
    for (i = 0; i < n; i++)
        printf("%s+", border);
    if (map & PRINTTARGET)
        printf("----+");
    if (!(map & NOLAST))
        printf("%s+\n", border);
}


void print_xvector (double *x, int n){
    char lena[DIMENSION];
    int i, length;

    length = 0;
    for (i=0; i<n; i++)
        if ((lena[i] = get_length(*(x+i))) > length)
            length = lena[i];

    if (length < DBL_DIG/2)
        length = DBL_DIG/2;
    else if (length < DBL_DIG)
        length += 2;

    for (i = 0; i < n; i++) {
        printf(" x[%u] = ", i+1);
        print_number(*(x + i), length, lena[i], 0);
        putchar(';');
    }
}


int get_length (double val){
    double temp;
    int length;
    char s[2 * DBL_DIG];

    modf(fabs(val), &temp);
    sprintf(s, "%.f", temp);
    length = strlen(s);
    if (val < 0)
        length++;
    return length;
}


void print_number (double num, int a, int b, char map) {
    char s[20];
    int i;

    sprintf(s, "%*.*f", a, a-b, num);
    i = strlen(s);
    while (s[i-1] == '0')
        i--;
    if(s[i-1] == '.')
        i--;
    s[i] = '\0';
    if (map & FIXED)
        a += 2;
    else
        a = strlen(s);
    if (fabs(atof(s)) == 0)
        printf("%*i", a, (int) atof(s));
    else
        printf("%*s", a, s);
}


int get_objective(symplex *obj, char map) {
    const int DIM = 010;
    unsigned int i, j;
    char s[DIMENSION];

    if (map & DIM) {
        printf("\nEnter the number of variables (n = ) ");
        scanf("%u", &obj->n);
        printf("\nEnter the number of constraints (m = ) ");
        scanf("%u", &obj->m);
    }
    if (map & A & obj->mem) {
        printf("\nFilling the matrix of linear constraints (A)");
        for (i=0; i < obj->m; i++) {
            printf("\nEnter %u line of A\n", i+1);
            for (j=0; j < obj->n; j++) {
                //printf("a[%u,%u] = ", i+1, j+1);
                scanf("%lf", obj->a+i*obj->n+j);
            }
            scanf("%s", s);
            if (strstr(s, ">=") != NULL)
                obj->signs[i] = BIGGER;
            else if (strstr(s, "<=") != NULL)
                obj->signs[i] = LOWER;
            else if (strstr(s, "=") != NULL)
                obj->signs[i] = EQUAL;
            sscanf(strstr(s, "=")+1, "%lf", obj->b+i);
        }
    }

/*
    if (map & B & obj->mem) {
        printf("\nEnter the B vector of linear constraints\n");
        for (i=0; i < obj->m; i++){
            printf("b[%u] = ", i+1);
            scanf("%lf", obj->b+i);
        }
    }
*/


    if (map & C & obj->mem) {
        printf("\nEnter the target function (C)\n");
        for (i=0; i < obj->n; i++)
            //printf("c[%u] = ", i+1);
            scanf("%lf", obj->c+i);
        printf("\n1 - maximize\n2 - minimize\nEnter the objective: ");
        scanf("%u", &j);
        obj->signs[DIMENSION-1] = (j == 1 ? MAX : MIN );
    }

    return 0;
}



int edit_objective(symplex *obj) {
    unsigned int i,j;
    char c;

    printf("\nEditing mode\n a, b, c - for editing appropriate coefficients.\n f - for changing type of target function.\n s - for changing type of linear constraints.\n\n");
    do {
        print_symplex(obj, PRINTTARGET | BINX);
        printf("\n\nCommand: ");
        switch (c = getchar()){
            case 'a':
                printf("\nEnter row: ");
                scanf("%u", &i);
                if ((i > 0) && (i <= obj->m)) {
                    printf("\nEnter col: ");
                    scanf("%u", &j);
                    if ((j > 0) && (j <= obj->n)) {
                        printf("\nEnter new value of A[%u,%u] = ", i, j);
                        scanf("%lf", obj->a+(i-1)*obj->n+(j-1));
                    } else
                        printf("\nIndex out of bounds.\n");
                } else
                    printf("\nIndex out of bounds.\n");
                break;
            case 'b':
                printf("\nEnter row: ");
                scanf("%u", &i);
                if ((i > 0) && (i <= obj->m)) {
                    printf("\nEnter new value of B[%u] = ", i);
                    scanf("%lf", obj->b+(i-1));
                } else
                    printf("\nIndex out of bounds.\n");
                break;
            case 'c':
                printf("\nEnter col: ");
                scanf("%u", &j);
                if ((j > 0) && (j <= obj->n)) {
                    printf("\nEnter new value of C[%u] = ", j);
                    scanf("%lf", obj->c+(j-1));
                } else
                    printf("\nIndex out of bounds.\n");
                break;
            case 'f':
                printf("\nEnter new type of target function:\n 1 - max\n 2 - min\n Your choice: ");
                scanf("%u", &j);
                obj->signs[DIMENSION-1] = (j == 1 ? MAX : MIN );
                break;
            case 's':
                printf("\nEnter number of constraint, which sign you want to change: ");
                scanf("%u", &i);
                printf("\nSelect type of constraint:\n 1. >=\n 2. <=\n 3. =\nYour choice: ");
                scanf("%u", &j);
                switch (j) {
                    case 1: obj->signs[i-1] = BIGGER; break;
                    case 2: obj->signs[i-1] = LOWER; break;
                    case 3: obj->signs[i-1] = EQUAL; break;
                }
                break;
            default:
                if (c != 'q')
                    printf("Invalid command %c\n", c);
        }
        getchar();
    } while (c != 'q');

    return 0;
}



