//reads in a file to a matrix given a char string
void readfile(char * NAME, struct matrix_t * matrix);

//this is a utility function for reading in the matrix file type.
//this is so we can allocate the proper storage size before reading it
// into the matrix.
int getmatrixsize(char * NAME, int * m, int * n);

//prints a given array to std::out
void printmatrix(struct matrix_t * matrix);

//prints array of doubles to standard out.
void printarray(double * array, int size);

//prints array of ints to standard out.
void printarrayi(int * array, int size);

//returns the trace of a given square matrix
double trace(struct matrix_t * matrix);

//returns the euclidian norm of an array
double euclidian_norm_array(double * array, int size);

//returns the euclidian norm of a matrix
double euclidian_norm(struct matrix_t * matrix);

//copys a matrix "from" into the matrix "to"
void copymatrix(struct matrix_t * from, struct matrix_t * to);

//subtracts two arrays elmentwise and stores in "sol":  sol = A-B
void arraysubtract(double * A, double * B, double * sol, int size);

//adds two arrays elmentwise and stores in "sol"; sol = A+B
void arrayadd(double * A, double * B, double * sol, int size);

// subtracts two matrices element wise and stores in sol sol = A-B
void matrixsubtract(struct matrix_t * A, struct matrix_t * B, struct matrix_t * sol);

//adds two matrices element wise and stores in sol sol = A+B
void matrixadd(struct matrix_t * A, struct matrix_t * B, struct matrix_t * sol);

//multiplies two matrices and stores in sol. Has checks to make sure A and B
//as well as the sol matrix are the correct size. sol = A*B
void matrixmultiply(struct matrix_t * A, struct matrix_t * B, struct matrix_t * sol);

//multiplies a matrix A by a vector B and stores in sol. Has checks to make sure A and B
//as well as the sol matrix are the correct size. sol = Ab
void matrixvecmultiply(struct matrix_t * A, double * B, double * sol);

//Given a system Ax=B, where a is a matrix and B is an array, this computes the
//gaussian elimination of A and stores it back in A. This will result in an upper
//triangular matrix. If this method fails, the return will be true for A is a
//singular matrix. If false is returned, A is succesfully in GE form.
bool gaussian_elimination(struct matrix_t * A, double * B);

//Given a matrix U in GE form, this backwards subsitutes solving Ux = y.
//where x and y are vectors. This will fail if U is singular.
void backwards_subsitution(struct matrix_t * U, double * y, double * x);

//Given a matrix A, this computes the LU decomposition and stores it back in A.
//Where L is the lower diagonal matrix and 1s on the diagonal and U is the
//upper diagonal with the actual values store in the diagonal. There is also a
//permutation vector s one must provide to keep track of how we permute A to
//get it into this form. This method also returns singular as a bool, failing if
//the matrix supplied is singular and returning true.
bool LU_decomposition(struct matrix_t * A, int * s);

//This does LU backsubsitution given a matrix A in LU form like in the provided
//method LU_decomposition(). This solves the system Ax=b where s it the
//permutation vector provided in LU_decomposition, b is the rhs vector, and x
//is a vector to hold the solutions.
void LU_backsubsitution(struct matrix_t * A, double *b, double * x, int * s);
