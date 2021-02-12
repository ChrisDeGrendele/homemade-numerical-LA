#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

//#include <iostream>
#include "LinAlg.H"

//The primary struct we will be using to store matrices. A matrix can
//be indexed like matrix.data[column + matrix.n * row]
typedef struct matrix_t {
    int m;
    int n;
    double * data;
} matrix_t;


int main()
{
    char filenameA[] = "Amat.dat";
    char filenameB[] = "Bmat.dat";

    //BUILD AND ALLOCATE MATRIX A
    int m_a;
    int n_a;
    getmatrixsize(&filenameA[0], &m_a, &n_a);
    //build matrix type
    double * matrix_ptrA = (double *)malloc(m_a*n_a*sizeof(double));
    static struct matrix_t matrix_A;
    matrix_A.m = m_a;
    matrix_A.n = n_a;
    matrix_A.data = matrix_ptrA;
    readfile(&filenameA[0], &matrix_A);

    //A"save"
    double * matrix_ptrAs = (double *)malloc(m_a*n_a*sizeof(double));
    static struct matrix_t matrix_As;
    matrix_As.m = m_a;
    matrix_As.n = n_a;
    matrix_As.data = matrix_ptrAs;
    readfile(&filenameA[0], &matrix_As);

    //BUILD AND ALLOCATE MATRIX B
    int m_b;
    int n_b;
    getmatrixsize(&filenameB[0], &m_b, &n_b);
    //build matrix type
    double * matrix_ptrB = (double *)malloc(m_b*n_b*sizeof(double));
    static struct matrix_t matrix_B;
    matrix_B.m = m_b;
    matrix_B.n = n_b;
    matrix_B.data = matrix_ptrB;
    readfile(&filenameB[0], &matrix_B);

    double * matrix_ptrBs = (double *)malloc(m_b*n_b*sizeof(double));
    static struct matrix_t matrix_Bs;
    matrix_Bs.m = m_b;
    matrix_Bs.n = n_b;
    matrix_Bs.data = matrix_ptrBs;
    readfile(&filenameB[0], &matrix_Bs);

    //QUESTION 2
    printf("---------------------PROBLEM 2-----------------------\n");
    printf("-----------------------------------------------------\n");
    printmatrix(&matrix_A);
    printf("Trace: %lf \n", trace(&matrix_A));

    for (int i=0; i < matrix_A.n; i++){ // loop over columns
        double column_vector[matrix_A.m];
        for (int j=0; j < matrix_A.m; j++){ //loop over rows
            column_vector[j] = matrix_A.data[i + matrix_A.n*j];
        }
        printf("Norm of %d th column: %lf \n", i+1,
                euclidian_norm_array(&column_vector[0],matrix_A.m));
    }



    //QUESTION 3
    printf("---------------------PROBLEM 3-----------------------\n");
    printf("-----------------------------------------------------\n");


    printf("Before Gaussian Elimination: \n");
    printf("A: \n");
    printmatrix(&matrix_A);
    printf("B: \n");
    printmatrix(&matrix_B);


    //BUILD SOLUTION MATRIX X
    double * matrix_ptrX = (double *)malloc(m_b*n_b*sizeof(double));
    static struct matrix_t matrix_X;
    matrix_X.m = m_b;
    matrix_X.n = n_b;
    matrix_X.data = matrix_ptrX;



    //make column vectors from B

    for (int i=0; i < matrix_B.n; i++){
        double column_vector[matrix_B.m];
        for (int j=0; j < matrix_B.m; j++){
            column_vector[j] = matrix_B.data[i + matrix_B.n*j];
        }
        bool singular;

        double solution[matrix_B.m];

        //printf("B: \n");
        //printarray(&column_vector[0], matrix_B.m);

        singular = gaussian_elimination(&matrix_A, &column_vector[0]);

        printf("After Gaussian Elimination \n");
        printf("A: \n");
        printmatrix(&matrix_A);
        printf("Vector b, column number: %d \n", i+1);
        printarray(&column_vector[0], matrix_B.m);
        printf("\n");


        backwards_subsitution(&matrix_A, &column_vector[0], &solution[0]);

        //reset A once we get solution vector x
        copymatrix(&matrix_As, &matrix_A);


        //store solution vector in x
        for (int j = 0; j < matrix_B.m; j++){
            matrix_X.data[i+matrix_X.n*j] = solution[j];
        }

    }

    printf("MATRIX X (SOLUTION):\n");
    printmatrix(&matrix_X);


    // E = AX-B

    //BUILD SOLUTION MATRIX E
    double * matrix_ptrE = (double *)malloc(m_b*n_b*sizeof(double));
    static struct matrix_t matrix_E;
    matrix_E.m = m_a;
    matrix_E.n = n_b;
    matrix_E.data = matrix_ptrE;

    printf("ERROR MATRIX AX-B\n");
    matrixmultiply(&matrix_A, &matrix_X, &matrix_E);
    matrixsubtract(&matrix_E, &matrix_B, &matrix_E);
    printmatrix(&matrix_E);

    for (int i=0; i < matrix_E.n; i++){ // loop over columns
        double column_vector[matrix_E.m];
        for (int j=0; j < matrix_E.m; j++){ //loop over rows
            column_vector[j] = matrix_E.data[i + matrix_E.n*j];
        }
        printf("Norm Error of %d th column: %lf \n", i+1,
                euclidian_norm_array(&column_vector[0],matrix_A.m));
    }


    //reset matrices
    copymatrix(&matrix_Bs, &matrix_B);
    copymatrix(&matrix_As, &matrix_A);



    //QUESTION 4
    printf("---------------------PROBLEM 4-----------------------\n");
    printf("-----------------------------------------------------\n");

    printf("Matrix A before LU decompositon:\n");
    printmatrix(&matrix_A);

    //LU Decomposition
    int permutation[matrix_A.m];
    bool singular;
    singular = LU_decomposition(&matrix_A, &permutation[0]);
    assert(singular == false);

    printf("Matrix A after LU decompositon:\n");
    printmatrix(&matrix_A);
    printf("\n");

    //calculate L and U to print
    double * matrix_ptrL = (double *)malloc(matrix_A.m*matrix_A.n*sizeof(double));
    static struct matrix_t matrix_L;
    matrix_L.m = matrix_A.m;
    matrix_L.n = matrix_A.n;
    matrix_L.data = matrix_ptrL;

    double * matrix_ptrU = (double *)malloc(matrix_A.m*matrix_A.n*sizeof(double));
    static struct matrix_t matrix_U;
    matrix_U.m = matrix_A.m;
    matrix_U.n = matrix_A.n;
    matrix_U.data = matrix_ptrU;

    for (int j = 0; j < matrix_A.m; j++){ // row
        for (int i = 0; i < matrix_A.n; i++){ //column
            if (i == j){//diagonal

                matrix_L.data[i + matrix_A.n*j] = 1.0;
                matrix_U.data[i + matrix_A.n*j] = matrix_A.data[i + matrix_A.n*j];
            }
            else if (i > j){//upperdiagonal
                matrix_L.data[i + matrix_A.n*j] = 0.0;
                matrix_U.data[i + matrix_A.n*j] = matrix_A.data[i + matrix_A.n*j];
            }
            else{
                matrix_L.data[i + matrix_A.n*j] = matrix_A.data[i + matrix_A.n*j];
                matrix_U.data[i + matrix_A.n*j] = 0.0;

            }

        }
    }
    printf("Matrix L:\n");
    printmatrix(&matrix_L);
    printf("\n");
    printf("Matrix U:\n");
    printmatrix(&matrix_U);
    printf("\n");



    //make column vectors from B
    for (int i=0; i < matrix_B.n; i++){
        double column_vector[matrix_B.m];
        for (int j=0; j < matrix_B.m; j++){
            column_vector[j] = matrix_B.data[i + matrix_B.n*j];
        }


        double x[matrix_A.m];

        LU_backsubsitution(&matrix_A, &column_vector[0], &x[0], &permutation[0]);

        printf("Peforming backsubsitution with b: \n");
        printarray(&column_vector[0], matrix_A.m);

        printf("Solution x Obtained:\n");
        printarray(&x[0], matrix_A.m);

        // E = Ax-b
        double E[matrix_B.m];
        matrixvecmultiply(&matrix_As, &x[0], &E[0]);
        arraysubtract(&E[0], &column_vector[0], &E[0], matrix_B.m);
        double error = euclidian_norm_array(&E[0], matrix_B.m);
        printf("Which has the following error (|Ax-b|) %lf:\n", error);
        printf("\n");
    }



    //QUESTION 5
    printf("---------------------PROBLEM 5-----------------------\n");
    printf("-----------------------------------------------------\n");

    //build matrix type
    double * matrix_ptr = (double *)malloc(3*3*sizeof(double));
    static struct matrix_t matrix;
    matrix.m = 3;
    matrix.n = 3;
    matrix.data = matrix_ptr;

    //matrix.data[j + matrix.n*i]

    matrix.data[0 + matrix.n*0] = 1.0;
    matrix.data[1 + matrix.n*0] = 2.0;
    matrix.data[2 + matrix.n*0] = 1.0;

    matrix.data[0 + matrix.n*1] = -3.0;
    matrix.data[1 + matrix.n*1] = 2.0;
    matrix.data[2 + matrix.n*1] = 1.0;

    matrix.data[0 + matrix.n*2] = 3.14159265359;
    matrix.data[1 + matrix.n*2] = 2.71828182845;
    matrix.data[2 + matrix.n*2] = 1.0;

    double vec[3] = {3.0,5.0,-1.41421356237};
    double solu[3];

    singular = gaussian_elimination(&matrix, &vec[0]);
    backwards_subsitution(&matrix, &vec[0], &solu[0]);

    printf("The solution found:\n");
    printarray(&solu[0], 3);



    return 0;

}
