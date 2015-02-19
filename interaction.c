#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "defines.h"

#define FIXED 01
#define FRACTIONALS 040

const int _base[DIMENSION] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

int get_length (double, char *);

void print_value (double, int, int, char);

void print_border (char *, int, char);

int get_value(double *F);

int getSign();

int clear_stdin();


void print_symplex(symplex *sym, char map){
    double *c, *b;

    if (map & PRINTLAST) {
		// TODO: does this part executed at least once
		printf("\nPRINTLAST\n");
        c = sym->d;
        b = sym->x;
    } else {
        c = sym->c;
        b = sym->b;
    }

    print_main (sym->a, b, c, sym->m, sym->n, 
		map & PRINTBASE ? sym->base : _base, 
		map & PRINTTARGET ? sym->signs : NULL, 
		0, 0, 0, map);
}


void print_main(double *a, double *x, double *d, int m, int n, const int base[], char signs[], int r, int s, double F, char map) {
    int length, i, j, f_length;
    char border[2 * DBL_DIG], *lena;

    length = 0;
    lena = (char *) calloc(m*n + 2*n, sizeof(char));			// There will be stored integer parts lengths of all coefficients 
	
	// Calculating length of coefficients and
	// Finding longest value, so all cells will have width equal to length of this value
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++)
            if ((*(lena+ i*n + j) = get_length(*(a + i*n + j), &map)) > length)
                length = *(lena+ i*n +j);
				//printf("\nLENG = %i\n", length);
        if ((*(lena+ m*n + i) = get_length(*(x + base[i]), &map)) > length)
            length = *(lena+ m*n + i);
    }
	
    for (i=0; i<n; i++)
        if ((*(lena+ m*n + n + i) = get_length(*(d + i), &map)) > length)
            length = *(lena+ m*n + n + i);
	
	if ((f_length = get_length(F, NULL)) > length)
		length = f_length;
	
	//printf("\nLENG = %i\n", length);
	
	// if length is too short we increase it a little (by increasing <length>, we increase precision of printing)
	if ( !(map & FRACTIONALS) ) {
		if (length < ((int)floor(log10(n))+2))
			length = (int)floor(log10(n))+2;
    } else if (length < DBL_DIG/3)
        length = DBL_DIG/3;
    else if (length < DBL_DIG)
        length += 2;
	
	// Calculating border size
    for (i = 0; i < length+3; i++)
        border[i] = '-';
    border[length+3] = '\0';

	// Printing header
    print_border(border, n, map);
    if (map & PRINTBASE)
        printf("| BS |%*c ", length+2, map & BINX ? 'B' : 'X');				// If it's iteration printing, print <Base> and <X> cols at start.
    putchar('|');
    for (i = 0; i < n; i++)
        printf("%*s%u |", length-(int)floor(log10(i+1))+1, "A", i+1);		// <A> vectors cols.
    if (map & PRINTTARGET)
        printf("    |");													// Gap for equation signs
    if (!(map & PRINTBASE))
        printf("%*c |", length+2, map & BINX ? 'B' : 'X');					// If it's objective printing, print <B> col at the end.
    putchar('\n');
    print_border(border, n, map);

	// Printing values
    for (i = 0; i < m; i++) {
		// If it's iteration printing, print <Base> and <X> cols at start.
        if (map & PRINTBASE) {
            printf("| %*u |", 2, base[i]+1);
            print_value(*(x + base[i]), length, *(lena + m*n + i), FIXED);
            putchar(' ');
        }
        putchar('|');
        for (j = 0; j < n; j++) {
            print_value(*(a + i*n + j), length, *(lena + i*n + j), FIXED);
			// Mark element, which will be removed from the basis
            if ((map & PRINTELEMENT) && (r == i) && (s == j))
                putchar('*');
            else
                putchar(' ');
            putchar('|');
        }
		// Print equation signs, if it's objective printing
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
		// Print Vector B in last cells
        if (!(map & PRINTBASE)) {
            print_value(*(x + base[i]), length, *(lena + m*n + i), FIXED);
            printf(" |");
        }
        putchar('\n');
    }
	
	// Printing footer
	print_border(border, n, map);
    if (map & PRINTBASE) {
        printf("| F= |");
        print_value(F, length, f_length, FIXED);
        putchar(' ');
    }

    putchar('|');
    for (i = 0; i < n; i++) {
        print_value(*(d + i), length, *(lena + m*n + n + i), FIXED);
        printf(" |");
    }
    if (map & PRINTTARGET)
        printf(" -> |%-*s|", length+3, (signs[DIMENSION-1] & MAX ? " max " : " min "));
    putchar('\n');
    print_border(border, n, map);
	
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
	printf("%s+\n", border);		// Another one border for <X> or <B> cols.
}


void print_xvector (double *x, int n){
    char lena[DIMENSION];
    int i, length;

    length = 0;
    for (i=0; i<n; i++)
        if ((lena[i] = get_length(*(x+i), NULL)) > length)
            length = lena[i];

    if (length < DBL_DIG/3)
        length = DBL_DIG/3;
    else if (length < DBL_DIG)
        length += 2;

    for (i = 0; i < n; i++) {
        printf(" x[%u] = ", i+1);
        print_value(*(x + i), length, lena[i], 0);
        putchar(';');
    }
}


// Return (lenght of string with integer of val) + 2. Ex: 123.3456 -> length('123') = 3
int get_length (double val, char *map){
    double integer;
    int length;
    char s[2 * DBL_DIG];

	modf(fabs(val), &integer);
	if (map != NULL && !is_integer(val)) {
		*map |= FRACTIONALS;
		//printf("FRACT: %20.18lf\n", val);
	}
    sprintf(s, "%.f", integer);
    length = strlen(s);
    if (val < 0)
        length++;
	//printf("VAL: %20.18lf = %i\n", val, length);
    return length;
}


void print_value (double num, int a, int b, char map) {
    char s[2 * DBL_DIG];
    int i;
	
    sprintf(s, "%*.*f", a, a-b, num);	// Convert double val to string
    i = strlen(s);
	if (strchr(s, '.') != NULL) {		// Value is fractional, so we need to make it nicer
		while (s[i-1] == '0')			// Remove trailing zeros from end of fractpart.
			i--;
		if(s[i-1] == '.')				// if whole fractpart consist of zeros remove leading '.', so there will be only integer part
			i--;
		s[i] = '\0';
	}
    if (map & FIXED)					// Element will be printed in table's cell, so we add more precision (two additional digits will be printed)
        a += 2;
    else
        a = strlen(s);
    if (fabs(atof(s)) == 0)				// If rest of val string is only "-0" change it to "0"
        printf("%*i", a, (int) atof(s));
    else
        printf("%*s", a, s);
}


int get_objective(symplex *obj, char map) {
    unsigned int i, j;

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
            for (j=0; j < obj->n; j++)					// Filling <A> row
				get_value(obj->a+i*obj->n+j);
			obj->signs[i] = getSign();					// Getting sign type {>=, <=, =}
			get_value(obj->b+i);						// Filling <B> value
        }
    }

    if (map & C & obj->mem) {
        printf("\nEnter the target function (C)\n");
        for (i=0; i < obj->n; i++)
			get_value(obj->c+i);
		clear_stdin();									// Remove all junk from input buffer
        printf("\n1 - maximize\n2 - minimize\nEnter the objective: ");
        scanf("%u", &j);
        obj->signs[DIMENSION-1] = (j == 1 ? MAX : MIN );
    }

    return SUCCESS;
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
						get_value(obj->a+(i-1)*obj->n+(j-1));
                        clear_stdin();
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
					get_value(obj->b+(i-1));
					clear_stdin();
                } else
                    printf("\nIndex out of bounds.\n");
                break;
            case 'c':
                printf("\nEnter col: ");
                scanf("%u", &j);
                if ((j > 0) && (j <= obj->n)) {
                    printf("\nEnter new value of C[%u] = ", j);
                    get_value(obj->c+(j-1));
					clear_stdin();
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

    return SUCCESS;
}


int get_value(double *F) {
	while (scanf("%lf", F) != 1)
		getchar();

	return SUCCESS;
}


int getSign() {
	char s[2];
	
	while (strchr(s,'=') == NULL) {
		scanf("%[><=]", s);
		getchar();
	}

	if (strstr(s, ">=") != NULL)
        return BIGGER;
	else if (strstr(s, "<=") != NULL)
		return LOWER;
	else if (strstr(s, "=") != NULL)
		return EQUAL;
	else
		return 0;
}


int clear_stdin() {
    while (getchar()!='\n')
		;
    return 1;
}