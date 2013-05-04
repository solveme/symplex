#include <stdio.h>
#include <string.h>
#include <float.h>
#include "defines.h"

int main(void) {
    enum { OBJ = 01, SYM = 02, PRC = 04, BAS = 010, SCS = 020, EDT = 040};
    char s, state;
    int t;
    double F;
    symplex obj, sym;

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
                get_objective(&obj, 010);
                if (get_mem(&obj, A | B | C | X) == 0) {
                    get_objective(&obj, 07);
                    state |= OBJ;
                }

                break;
            case '2':
				if (state & OBJ) {
                    getchar();
                    edit_objective(&obj);
                    if (state & SYM)
                        clean_mem(&sym);
                    state = OBJ | EDT;
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
                    if ((state & SYM) == 0) {
                        sym.n = obj.n;
                        sym.m = obj.m;
                        if (get_mem(&sym, A | B | C | D | X ) == 0)
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
	                    if (get_first_plan(&sym)) {
	                        printf("\nNo basic plan\n");
	                    } else {
                            state |= BAS;
	                        printf("\nBasic plan:");
	                        print_xvector(sym.x, sym.n);
                            putchar('\n');
                            putchar('\n');
                            get_deltas(&sym);
	                        if (iteration(&sym) == 0) {
                                copy_symplex(&sym, &obj, obj.n, X);
	                            printf("\nSolved:");
                                print_xvector(obj.x, obj.n);
                                printf(" F = ");
                                F = find_value(obj.x, obj.c, obj.n);
                                t = get_length(F);
                                print_number(F, (t < DBL_DIG/3 ? DBL_DIG/3 : t + 2), t, 0);
                                putchar('\n');
                                state |= SCS;
	                        } else
                                printf("\nNo solutions\n");

	                    }
                        state |= PRC;
                    }
                } else
                	printf("\nNo data to process\n");
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
    return 0;
}
