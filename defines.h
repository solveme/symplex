#define DIMENSION 16

typedef struct ad {
    double *a;
    double *d;
} short_symplex;

typedef struct abcdx {
    char mem;
    unsigned int m;
    unsigned int n;
    double *a;
    double *b;
    double *c;
    double *d;
    double *x;
    char signs[DIMENSION];
    int base[DIMENSION];
} symplex;

enum { A = 01, B = 02, C = 04, D = 010, X = 020, BASE = 040, SIGNS = 0100 };
enum { BIGGER = 01, LOWER = 02, EQUAL = 04, MAX = 010, MIN = 020 };
enum { PRINTBASE = 01, PRINTTARGET = 02, PRINTELEMENT = 04, PRINTLAST = 010, BINX = 020};

int get_mem (symplex *, char);
void clean_mem (symplex *);
void copy_symplex (symplex *, symplex *, int, char);
void get_deltas (symplex *);
int get_first_plan (symplex *);
void make_standard (symplex *);
int make_canonic (symplex *);
double find_value (double *, double *, int);

void print_symplex (symplex *, char);
void print_main (double *, double *, double *, int, int, int [], char [], int, int, double, char);
void print_xvector (double *, int);
int get_length (double);
void print_number (double, int, int, char);
int get_objective (symplex *, char);
int edit_objective(symplex *);

int iteration (symplex *);

