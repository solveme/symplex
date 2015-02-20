#include <stdio.h>
#include <string.h>
#include <float.h>
#include "defines.h"

symplex obj;

int main(void) {
    enum { OBJ = 01, SYM = 02, EDT = 040};
    char s, state;
    int t;
    double F;
    symplex sym;

    obj.mem = sym.mem = state = 0;
    
    printf("Welcome to symplex method program.\n");
    printf(" 1 - Define new objective\n 2 - Edit current objective\n 3 - Print objective\n 4 - Find a solution\n q - Quit\n");

    do {
        printf("\nCommand: ");
        switch (s = getchar()) {
            case '1':
                if (state & OBJ)
                    clean_mem(&obj);
                if (state & SYM)
                    clean_mem(&sym);
                state = 0;
                get_objective(&obj, DIM);
                if (get_mem(&obj, A | B | C | X) == SUCCESS) {
                    get_objective(&obj, A | B | C);
                    state |= OBJ;
                }

                break;
            case '2':
                if (state & OBJ) {
                    getchar();              // This will remove <Enter> character from input buffer, so next call of getchar() will wait until another button is pressed
                    edit_objective(&obj);
                    if (state & SYM)
                        clean_mem(&sym);    // Clear all current computations
                    state = OBJ | EDT;      // Set flag, that there is no need to release <Enter> character at the end of _while_ loop
                }
                break;
            case '3':
                if (state & OBJ) {
                    printf("\nObjective\n");
                    print_symplex(&obj, PRINTTARGET | BINX);
                    putchar('\n');
                } else
                    printf("\nObjective is not defined\n");
                break;
            case '4':
                if (state & OBJ) {
                    printf("\nSearch for solutions\n");
                    if ( !(state & SYM) ) {
                        sym.n = obj.n;
                        sym.m = obj.m;
                        if (get_mem(&sym, A | B | C | D | X ) == SUCCESS)
                            state |= SYM;
                        else {
                            printf("\nError: Can't get enough memory");
                            clean_mem(&sym);
                            break;
                        }
                    }
                    if (state & SYM) {
                        sym.n = obj.n;
                        sym.m = obj.m;
                        copy_symplex(&obj, &sym, 0, A | B | C | SIGNS);
                        printf("\nReduction to the standard form\n");
                        make_standard(&sym);
                        print_symplex(&sym, PRINTTARGET | BINX);
                        printf("\n\nReduction to the canonical form\n");
                        make_canonic(&sym);
                        print_symplex(&sym, PRINTTARGET | BINX);
                        printf("\n\nReceiving a basic plan\n");
                        if (get_first_plan(&sym) != SUCCESS) {
                            printf("\nNo basic plan\n");
                        } else {
                            printf("\nBasic plan:");
                            print_xvector(sym.x, sym.n);
                            putchar('\n');
                            putchar('\n');
                            get_deltas(&sym);
                            if (iteration(&sym) == SUCCESS) {
                                copy_symplex(&sym, &obj, obj.n, X);
                                printf("\nSolved:");
                                print_xvector(obj.x, obj.n);
                                printf(" F = ");
                                F = find_value(obj.x);
                                t = get_length(F, NULL);
                                print_value(F, (t < DBL_DIG/3 ? DBL_DIG/3 : t + 2), t, 0);
                                putchar('\n');
                            } else
                                printf("\nThere is no solutions\n");

                        }
                    }
                } else
                    printf("\nThere is no data to process\n");
                break;
            default:
                if (s != 'q')
                    printf("Invalid command %c\n", s);
        }
        state & EDT ? state &= ~EDT : getchar();
    } while (s != 'q');

    if (state & OBJ)
        clean_mem(&obj);
    if (state & SYM)
        clean_mem(&sym);
    return SUCCESS;
}


double find_value(double *x) {
    int i;
    double S = 0;

    for (i = 0; i < obj.n; i++)
        S += *(x + i) * *(obj.c + i);

    return S;
}
