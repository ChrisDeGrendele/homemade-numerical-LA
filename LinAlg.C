#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
//#include <iostream>

//The primary struct we will be using to store matrices. A matrix can
//be indexed like matrix.data[column + matrix.n * row]
typedef struct matrix_t {
    int m;
    int n;
    double * data;
} matrix_t;



void arraysubtract(double * A, double * B, double * sol, int size){
    for (int i = 0; i < size; i++){
        sol[i] = A[i] - B[i];
    }
}

void arrayadd(double * A, double * B, double * sol, int size){
    for (int i = 0; i < size; i++){
        sol[i] = A[i] + B[i];
    }
}


void matrixsubtract(struct matrix_t * A, struct matrix_t * B, struct matrix_t * sol){

    assert(A->m == B->m);
    assert(B->m == sol->m);
    assert(A->n == B->n);
    assert(B->n == sol->n);


    for (int i = 0; i< A->m; i++){
        for (int j = 0; j < A->n; j++){
            sol->data[j+A->n*i] = A->data[j+A->n*i] - B->data[j+A->n*i];
        }
    }
}



void matrixadd(struct matrix_t * A, struct matrix_t * B, struct matrix_t * sol){

    assert(A->m == B->m);
    assert(B->m == sol->m);
    assert(A->n == B->n);
    assert(B->n == sol->n);


    for (int i = 0; i< A->m; i++){
        for (int j = 0; j < A->n; j++){
            sol->data[j+A->n*i] = A->data[j+A->n*i] + B->data[j+A->n*i];
        }
    }
}


void matrixvecmultiply(struct matrix_t * A, double * B, double * sol){


    for (int i = 0; i < A->m; i++){// loop through row of A
        double dotproduct = 0.0;
        for (int j = 0; j < A->n; j++){//loop through columns of B
            dotproduct += A->data[j+A->n*i] * B[j];
        }
        sol[i] = dotproduct;
    }

}


void matrixmultiply(struct matrix_t * A, struct matrix_t * B, struct matrix_t * sol){

    assert(A->n == B->m);
    assert(sol->m == A->m);
    assert(sol->n == B->n);

    for (int i = 0; i < A->m; i++){// loop through row of A

        for (int k = 0; k < B->n; k++){//loop through columns of B

            double dotproduct = 0.0;

            for (int j = 0; j < B->m; j++){ //loop through columns of B

                dotproduct += A->data[j+A->n*i] * B->data[k+B->n*j];

            }

            sol->data[k + sol->n*i] = dotproduct;

        }
    }

}

void copymatrix(struct matrix_t * from, struct matrix_t * to){
    assert(from->m == to->m);
    assert(from->n == to->n);

    for (int j = 0; j < from->m; j++){
        for (int i = 0; i < from->n; i++){
            to->data[i + to->n*j] = from->data[i + from->n*j];
        }
    }
}

void printarray(double * array, int size){
    for (int i = 0; i < size; i++){
        printf("%lf\t", array[i]);
    }
    printf("\n");
}

void printarrayi(int * array, int size){
    for (int i = 0; i < size; i++){
        printf("%d\t", array[i]);
    }
    printf("\n");
}


void printmatrix(struct matrix_t * matrix){
    printf("%d x %d Matrix: \n ", matrix->m, matrix->n);
    for (int j = 0; j < matrix->m; j++){
        for (int i = 0; i < matrix->n; i++){
            printf("%lf\t", matrix->data[i + matrix->n*j]);
        }
        printf("\n");
    }
}


double trace(struct matrix_t * matrix){
    assert(matrix->m == matrix->n);
    double trace = 0;
    for (int j = 0; j < matrix->m; j++){
        for (int i = 0; i < matrix->n; i++){
            if (i == j){
                trace += matrix->data[i + matrix->n*j];
            }
        }
    }
    return trace;
}

double euclidian_norm_array(double * array, int size){
    double norm = 0;
    for (int i = 0; i < size; i++){
            norm += pow(array[i],2);
    }
    return pow(norm,0.5);
}


double euclidian_norm_matrix(struct matrix_t * matrix){
    double norm = 0;
    for (int j = 0; j < matrix->m; j++){
        for (int i = 0; i < matrix->n; i++){
            if (i == j){
                norm += pow(matrix->data[i + matrix->n*j],2);
            }
        }
    }
    return pow(norm,0.5);
}

void getmatrixsize(char * NAME, int * m, int * n){
    // Opens a file and puts the size into pointer m and n
    FILE *fptr;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int linenum = 0;
    //array size

    if ((fptr = fopen(NAME,"r")) == NULL){
        printf("Error! opening file");
        exit(1);
    }


    while ((read = getline(&line, &len, fptr)) != -1) {
        if (linenum==0){
            sscanf(&line[0], "%d %d", m, n);
        }
        linenum++;
    }

    fclose(fptr);


}

void readfile(char * NAME, struct matrix_t * matrix){
    // Function for reading input files and returning data.
    FILE *fptr;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    int linenum = 0;
    //array size
    int m;
    int n;


    fptr = fopen(NAME,"r");
    linenum=0;
    fscanf(fptr, "%d %d", &m,&n);

    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            //fscanf(fptr, "%lf", &matrix[i][j]);
            fscanf(fptr, "%lf", (matrix->data + j*n + i));
        }
    }

}

bool gaussian_elimination(struct matrix_t * A, double * B){

    //Fails if singular is true.
    bool singular = false;


    for (int j = 0; j < A->m; j++){

        //find max |Aₖⱼ|.
        double max = -1;
        double Akj;
        int max_loc;
        for (int k = j; k < A->m; k++){
            Akj = A->data[k + A->n*j];
            if (abs(Akj) > max){
                max = Akj;
                max_loc = k;
            }
        }

        //if we found a higher max, we need to pivot.
        double temp;
        if (max_loc != j){
            for (int k = j; k < A->n; k++){
                //Swap A(j,k) with A(max_loc, k) note: A(row, column)
                temp = A->data[k + A->n*j];
                A->data[k + A->n*j] = A->data[k + A->n*max_loc];
                A->data[k + A->n*max_loc] = temp;
            }
            temp = B[j];
            B[j] = B[max_loc];
            B[max_loc] = temp;
        }

        //if singular, just break. GE fails.
        if (A->data[j + A->n*j] == 0){
            singular = true;
            break;
        }

        //Now do gaussian elimination.
        for (int i = j+1; i < A->m; i++){// loop down the column under diagonal

            //ri = ri -rj*Aij/Ajj
            double Aij = A->data[j+A->n*i];
            double Ajj = A->data[j+A->n*j];

            for (int k = 0; k < A->n; k++){//loop over whole row to apply scaling
                A->data[k + A->n*i] -= A->data[k + A->n*j]*Aij/Ajj;
            }
            B[i] -= B[j]*Aij/Ajj;

        }
    }

    return singular;

}

void backwards_subsitution(struct matrix_t * U, double * y, double * x){

    for (int i = U->m -1; i >= 0; i--){

        double sum = 0.0;
        for (int k = i+1; k < U->m; k++){
            sum += U->data[k + U->n*i]*x[k];
        }

        //If singular, this doesn't work!
        assert(U->data[i+U->n*i] != 0);

        x[i] = (y[i]-sum)/U->data[i + U->n*i];
    }
}


int LU_decomposition(struct matrix_t * A, int * s){

    bool singular = false;

    //Initialize permutation vector
    for (int j = 0; j < A->m; j++){
        s[j] = j;
    }


    for (int j = 0; j < A->m; j++){

        //find max |Aₖⱼ|.
        double max = -1;
        double Akj;
        int max_loc;
        for (int k = j; k < A->m; k++){
            Akj = A->data[k + A->n*j];
            if (abs(Akj) > max){
                max = Akj;
                max_loc = k;
            }
        }

        //if we found a higher max, we need to pivot.
        double temp;
        if (max_loc != j){
            for (int k = j; k < A->n; k++){
                //Swap A(j,k) with A(max_loc, k) note: A(row, column)
                temp = A->data[k + A->n*j];
                A->data[k + A->n*j] = A->data[k + A->n*max_loc];
                A->data[k + A->n*max_loc] = temp;
            }
            temp = s[j];
            s[j] = s[max_loc];
            s[max_loc] = temp;
        }

        if (A->data[j + A->n*j] == 0){
            singular = true;
            break;
        }


        for (int i = j+1; i<A->m; i++){


            //aij = aij/ajj
            A->data[j + A->n*i] = A->data[j + A->n*i] / A->data[j + A->n*j];

            for (int k = j+1; k < A->m; k++ ){
                //aik = aik -aij ajk
                A->data[k + A->n*i] -= A->data[j + A->n*i]* A->data[k + A->n*j];
            }
        }

    }

    return singular;
}


void LU_backsubsitution(struct matrix_t * A, double *b, double * x, int * s){
    //solve Ax = b with permutation vector s.

    //A should be in LU form. Where U is in the upper traingular and L is in the lower.


    //compute L and U. Note this is not necessary and it is slower, but it was
    //nice to have for debugging purposes so I'm leaving it for now.
    double * matrix_ptrL = (double *)malloc(A->m*A->n*sizeof(double));
    static struct matrix_t matrix_L;
    matrix_L.m = A->m;
    matrix_L.n = A->n;
    matrix_L.data = matrix_ptrL;

    double * matrix_ptrU = (double *)malloc(A->m*A->n*sizeof(double));
    static struct matrix_t matrix_U;
    matrix_U.m = A->m;
    matrix_U.n = A->n;
    matrix_U.data = matrix_ptrU;

    for (int j = 0; j < A->m; j++){ // row
        for (int i = 0; i < A->n; i++){ //column
            if (i == j){//diagonal

                matrix_L.data[i + A->n*j] = 1.0;
                matrix_U.data[i + A->n*j] = A->data[i + A->n*j];
            }
            else if (i > j){//upperdiagonal
                matrix_L.data[i + A->n*j] = 0.0;
                matrix_U.data[i + A->n*j] = A->data[i + A->n*j];
            }
            else{
                matrix_L.data[i + A->n*j] = A->data[i + A->n*j];
                matrix_U.data[i + A->n*j] = 0.0;

            }

        }
    }

    //initialize y
    double y[A->m];
    for (int j = 0; j < A->m; j++){
        y[j] = b[s[j]];
    }

    //Step 1 Forward Subsitution
    // Ly = b

    //a_ij = ith row and the jth column
    //data[column + n*row]

    for (int j = 0; j < A->m-1; j++){
        for (int i = j+1; i < A->m; i++){
            y[i] -= y[j]*matrix_L.data[j + A->n*i];
        }
    }


    //Step 2 backwards Subsitution
    // Ux = y
    for (int i = A->m-1; i >= 0; i--){

        //If singular, this doesn't work!
        assert(A->data[i+A->n*i] != 0);
        double sum = 0.0;

        for (int k = i+1; k < A->m; k++){
            sum += matrix_U.data[k+A->n*i]*x[k];
        }

        x[i] = (y[i]-sum)/matrix_U.data[i+A->n*i];

    }

}
