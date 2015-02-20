// add this for some tests of github
#define DIMENSION 16			// Number of constraints (max m in symplex structure)

#define SUCCESS 0
#define NOT_ENOUGH_MEMORY 01000
#define DIM 0200				// Flag for <get_objective function>: we should get dimensions of objective before we do any allocation

typedef struct ad {
    double *a;
    double *d;
} short_symplex;

typedef struct abcdx {
    char mem;					// bitmap that show for which vectors (a,b,c ...) memeory was allocated
    unsigned int m;				// number of costraints
    unsigned int n;				// number of variables
    double *a;					// matrix for constraints coefficients
    double *b;					// B vector
    double *c;					// C vector (coefficients of target function)
    double *d;					// Vector of deltas
    double *x;					// X vector (variable values)
    char signs[DIMENSION];      // Stores signs of equations (e.g. =, <, >=), signs[DIMENSION-1] (last cell of this Vector) used for storing type of objective (MIN or MAX)
    int base[DIMENSION];		// Stores which of A vectors are basic
} symplex;

enum { A = 01, B = 02, C = 04, D = 010, X = 020, BASE = 040, SIGNS = 0100};		// Flags, which represent appropriate vectors an matrix
enum { BIGGER = 01, LOWER = 02, EQUAL = 04, MAX = 010, MIN = 020 };				// Comparison signs
enum {	ITERATION_OPTIMAL = 0,			// Solution was found
		ITERATION_FAIL = 01,			// There is no solution for this objective
		ITERATION_NEXT = 02				// We should do another loop
};
enum { 	PRINTBASE = 01,					// Print basis col in table
		PRINTTARGET = 02,				// Print target functions coefficients in table
		PRINTELEMENT = 04,				// Mark element, which will be removed from basis
		PRINTLAST = 010,
		BINX = 020						// Instead of X print B vector in X's cells
	};

// Memory functions (functions.c)
int get_mem (symplex *, char);
void clean_mem (symplex *);
void copy_symplex (symplex *, symplex *, int, char);

// Computing functions (functions.c)
void get_deltas (symplex *);
int get_first_plan (symplex *);
void make_standard (symplex *);
int make_canonic (symplex *);
double find_value (double *);
int iteration (symplex *);
int is_integer (double);

// IO functions (interaction.c)
void print_symplex (symplex *, char);
void print_main (double *, double *, double *, int, int, const int [], char [], int, int, double, char);
void print_xvector (double *, int);
int get_length (double, char *);
void print_value (double, int, int, char);
int get_objective (symplex *, char);
int edit_objective(symplex *);



