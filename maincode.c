#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cilk/cilk.h>
#include "matrix_operations.h"

typedef struct {
    double *values;
    int *indptr;
    int *indices;
} CSRMatrix;

typedef struct {
    double *values;
    int *indptr;
    int *indices;
} CSCMatrix;

CSCMatrix initialize_csc_matrix() {
    CSCMatrix matrix;
    matrix.values = NULL;
    matrix.indices = NULL;
    matrix.indptr = NULL;
    return matrix;
}
CSRMatrix initialize_csr_matrix() {
    CSRMatrix matrix;
    matrix.values = NULL;
    matrix.indices = NULL;
    matrix.indptr = NULL;
    return matrix;
}

static const CSCMatrix ERROR_CSC_MATRIX = {NULL, NULL, NULL};
static const CSRMatrix ERROR_CSR_MATRIX = {NULL, NULL, NULL};
void printArray(int arr[], int size) {
    for (int i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}
void print_csr_matrix(CSRMatrix matrix, int c) {
    printf("Result in CSR form:\n");

    // Print values
    printf("Values: ");
    for (int i = 0; i < matrix.indptr[c]; i++) {
        printf("%f ", matrix.values[i]);
    }
    printf("\n");

    // Print indices
    printf("Indices: ");
    for (int i = 0; i  < matrix.indptr[c]; i++) {
        printf("%d ", matrix.indices[i]);
    }
    printf("\n");

    // Print indptr
    printf("Indptr: ");
    for (int i = 0; i <= c; i++) {
        printf("%d ", matrix.indptr[i]);
    }
    printf("\n");
}
void print_csc_matrix(CSCMatrix matrix, int c) {
    printf("Result in CSC form:\n");

    // Print values
    printf("Values: ");
    for (int i = 0; i < matrix.indptr[c]; i++) {
        printf("%f ", matrix.values[i]);
    }
    printf("\n");

    // Print indices
    printf("Indices: ");
    for (int i = 0; i < matrix.indptr[c]; i++) {
        printf("%d ", matrix.indices[i]);
    }
    printf("\n");

    // Print indptr
    printf("Indptr: ");
    for (int i = 0; i <= c; i++) {
        printf("%d ", matrix.indptr[i]);
    }
    printf("\n");
}

 
CSRMatrix sparse_matrix_multiplication(CSRMatrix matrix, int *vector, int c, int n) {

    int s;
    double el;
    CSCMatrix csc = initialize_csc_matrix();
    int val, start_idx, end_idx;
    int max_nnz = matrix.indptr[n];
    csc.values = (double *)malloc(max_nnz * sizeof(double));
    csc.indices = (int *)malloc(max_nnz * sizeof(int));
    csc.indptr = (int *)malloc((c + 1) * sizeof(int));

    if (csc.values == NULL || csc.indices == NULL || csc.indptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(csc.values);
        free(csc.indices);
        free(csc.indptr);
        return ERROR_CSR_MATRIX;
    }
    int num_threads = c;
    double** local_csc_values = (double **)malloc(num_threads * sizeof(double*));
    int** local_csc_indices = (int **)malloc(num_threads * sizeof(int*));
    int** local_csc_indptr = (int **)malloc(num_threads * sizeof(int*));

    for (int j = 0; j < num_threads; j++) {
        local_csc_values[j] = (double *)calloc(max_nnz, sizeof(double));
        local_csc_indices[j] = (int *)calloc(max_nnz, sizeof(int));
        local_csc_indptr[j] = (int *)calloc(1, sizeof(int));
    }
    csc.indptr[0] = 0;
    int new_col = csc.indptr[0];

    int* local_nnz = (int*)malloc(num_threads * sizeof(int));
     
    
     
    cilk_for (int col = 0; col <c; col++) {
        
        local_nnz[col] = 0;
	
        for (int i = 0; i < n; i++) {
            val = 0;
            start_idx = matrix.indptr[i];
            end_idx = matrix.indptr[i + 1];

            for (int j = start_idx; j < end_idx; j++) {
                s = matrix.indices[j];
                el = matrix.values[j];

                if (vector[s] == col) {
                    val += el;
                }
            }

            if (val != 0) {
		 
                local_csc_values[col][local_nnz[col]] = val;
                local_csc_indices[col][local_nnz[col]] = i;
                local_nnz[col]++;
		
            }
        }
	
	local_csc_indptr[col][0] = local_nnz[col];
    }

	

    
	//printArray(local_nnz,c);
	 

         //COMBINE thelv ordered idanika
        //cilk_sync;
       
	
        for (int t = 0; t < c; t++) {
             int new_col = csc.indptr[t];
            int nnz = local_csc_indptr[t][0];


            for (int j = 0; j < nnz; j++) {
                csc.values[new_col + j] = local_csc_values[t][j];
                csc.indices[new_col + j] = local_csc_indices[t][j];
            }

            new_col += nnz;
        

        csc.indptr[t + 1] = new_col;
    }

    //print_csc_matrix(csc,c);
    free(local_nnz);

    csc.values = realloc(csc.values, csc.indptr[c] * sizeof(double));
    csc.indices = realloc(csc.indices, csc.indptr[c] * sizeof(int));

    // Free the arrays of pointers
    for (int j = 0; j < num_threads; j++) {
        free(local_csc_values[j]);
        free(local_csc_indices[j]);
        free(local_csc_indptr[j]);
    }
    free(local_csc_values);
    free(local_csc_indices);
    free(local_csc_indptr);

   //2os ypologismos---------------------------
    num_threads=c;
    s=0;
    el=0.0;
    CSRMatrix csr=initialize_csr_matrix();
    int nnz_res=0;
    //max_nnz = nnz;  //to megisto nnz pou mporei na exei einai osa nnz eixe o csr prin.
    max_nnz=csc.indptr[c];
    csr.values=(double *)malloc(max_nnz * sizeof(double));
    csr.indices=(int *)malloc(max_nnz * sizeof(int));
    csr.indptr=(int *)malloc((c + 1) * sizeof(int));
    if (csr.values == NULL || csr.indices == NULL || csr.indptr == NULL) {
        fprintf(stderr, "Memory allocation failed\n");

        free(csr.values);
        free(csr.indices);
        free(csr.indptr);
        return ERROR_CSR_MATRIX;
    }

    double** local_csr_values = (double **)malloc(num_threads * sizeof(double*));
    int** local_csr_indices = (int **)malloc(num_threads * sizeof(int*));
    int** local_csr_indptr = (int **)malloc(num_threads * sizeof(int*));

    for (int j = 0; j < num_threads; j++) {
        local_csr_values[j] = (double *)calloc(max_nnz, sizeof(double));
        local_csr_indices[j] = (int *)calloc(max_nnz, sizeof(int));
        local_csr_indptr[j] = (int *)calloc(1, sizeof(int));
    }
    csr.indptr[0] = 0;
    int new_row = csr.indptr[0];

    local_nnz = (int*)malloc(num_threads * sizeof(int)); 
    cilk_for (int row = 0; row <c; row++) {
        
        local_nnz[row] = 0;
	
        for (int i = 0; i < c; i++) {
            val = 0;
            start_idx = matrix.indptr[i];
            end_idx = matrix.indptr[i + 1];

            for (int j = start_idx; j < end_idx; j++) {
                s = matrix.indices[j];
                el = matrix.values[j];

                if (vector[s] == row) {
                    val += el;
                }
            }

            if (val != 0) {
		 
                local_csr_values[row][local_nnz[row]] = val;
                local_csr_indices[row][local_nnz[row]] = i;
                local_nnz[row]++;
		
            }
        }
	
	local_csr_indptr[row][0] = local_nnz[row];
    }
	//combine
    for (int t = 0; t < c; t++) {
             int new_row = csr.indptr[t];
            int nnz = local_csr_indptr[t][0];


            for (int j = 0; j < nnz; j++) {
                csr.values[new_row + j] = local_csr_values[t][j];
                csr.indices[new_row + j] = local_csr_indices[t][j];
            }

            new_row += nnz;
        

        csr.indptr[t + 1] = new_row;
    }

   free(local_nnz);

    csr.values = realloc(csr.values, csr.indptr[c] * sizeof(double));
    csr.indices = realloc(csr.indices, csr.indptr[c] * sizeof(int));
 free(csc.values);
free(csc.indptr);
free(csc.indices);

//   print_csr_matrix(csr, c);

    for (int j = 0; j < num_threads; j++) {
    free(local_csr_values[j]);
    free(local_csr_indices[j]);
    free(local_csr_indptr[j]);
    }

    // Free the arrays of pointers
    free(local_csr_values);
    free(local_csr_indices);
    free(local_csr_indptr);
    return csr;
    
     
    



}
int main() {


    clock_t start_time = clock();
    int nrows, ncols, nnz;
    // CSR arrays
    int *row_ptr, *col_ind;
    double *values;

    //PARALLEL
    if (read_matrix_market_to_csr("mbeaflw.mtx", &nrows, &ncols, &nnz, &row_ptr, &col_ind, &values) != 0) {
        return 1;  // Error reading mat
    }

    // Now you can use the sparse matrix in CSR format for your functions
    CSRMatrix csr_matrix;
    CSRMatrix csr_result;
    csr_matrix.values = values;
    csr_matrix.indices = col_ind;
    csr_matrix.indptr = row_ptr;

    // Create a vector with random values between 0 and ncols - 1
    int FIXED_SEED=5;
    srand(FIXED_SEED); // Seed for random number generation
    int *vector = (int *)malloc(nrows * sizeof(int));
    for (int i = 0; i < nrows; i++) {
        vector[i] = rand()%ncols;
    }

    // Call the function from matrix_ops.c
    csr_result = sparse_matrix_multiplication(csr_matrix, vector, ncols, nrows);
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Elapsed time: %f seconds\n", elapsed_time);

    // Print the result
    //print_csr_matrix(csr_result, ncols);

    // Free allocated memory for CSR format
    free(row_ptr);
    free(col_ind);
    free(values);

    // Free the vector
    free(vector);

    // Free the allocated memory for the CSC matrix
    free(csr_result.values);
    free(csr_result.indices);
    free(csr_result.indptr);

    return 0;

}
