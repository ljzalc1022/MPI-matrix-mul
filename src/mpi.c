#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "matrix.h"

const int max_size = 512;
const int max_elements = 512 * 512;

int main(int argc, char *argv[]) {
    /** Matrix Properties
     * [0] = Rows of Matrix A
     * [1] = Cols of Matrix A
     * [2] = Rows of Matrix B
     * [3] = Cols of Matrix B
     **/
    int matrix_properties[4];
    double *local_A;
    double *local_B;
    double *local_C;
    int block_size;

    // used only at master
    int *row_counts;
    int *displs;
    
    int num_worker, rank;
    double start_time, end_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    /** the master initializes the data **/
    if (rank == 0) {

        if(argc != 3){
            printf("ERROR: Please specify only 2 files.\n");
            exit(EXIT_FAILURE);
        }
            
        matrix_struct *m_1 = get_matrix_struct(argv[1]);
        matrix_struct *m_2 = get_matrix_struct(argv[2]);

        if(m_1->cols != m_2->rows){
            printf("ERROR: The number of columns of matrix A must match the number of rows of matrix B.\n");
            exit(EXIT_FAILURE);
        }
        
        // fill the property-array for workers
        matrix_properties[0] = m_1->rows;
        matrix_properties[1] = m_1->cols;
        matrix_properties[2] = m_2->rows;
        matrix_properties[3] = m_2->cols;
        
        start_time = MPI_Wtime();

        MPI_Bcast(matrix_properties, 4, MPI_INT, 0, MPI_COMM_WORLD);

        int A_size = m_1->rows * m_1->cols;
        local_A = malloc(A_size * sizeof(double));
        for (int i = 0; i < m_1->rows; i++) {
            memcpy(local_A + i * m_1->cols, m_1->mat_data[i], m_1->cols * sizeof(double));
        }
        block_size = (m_1->rows + num_worker - 1) / num_worker;
        row_counts = malloc(num_worker * sizeof(int));
        displs = malloc(num_worker * sizeof(int));
        for (int i = 0; i < num_worker; i++) {
            displs[i] = i * block_size * m_1->cols;
            row_counts[i] = i != num_worker - 1 ? block_size * m_1->cols : A_size - displs[i];
        }
        MPI_Scatterv(local_A, row_counts, displs, MPI_DOUBLE, 
                     MPI_IN_PLACE, 0, 0, 0, MPI_COMM_WORLD);

        int B_size = m_2->rows * m_2->cols;
        local_B = malloc(B_size * sizeof(double));        
        for (int i = 0; i < m_2->rows; i++) {
            memcpy(local_B + i * m_2->cols, m_2->mat_data[i], m_2->cols * sizeof(double));
        }
        MPI_Bcast(local_B, B_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        free_matrix(m_1);
        free_matrix(m_2);
    }
    else {
        // to receive data from the master
        MPI_Bcast(matrix_properties, 4, MPI_INT, 0, MPI_COMM_WORLD);

        block_size = (matrix_properties[0] + num_worker - 1) / num_worker;
        if (rank == num_worker - 1) {
            block_size = matrix_properties[0] - (num_worker - 1) * block_size;
        }
        local_A = malloc(block_size * matrix_properties[1] * sizeof(double));
        MPI_Scatterv(NULL, NULL, NULL, 0, 
                     local_A, block_size * matrix_properties[1], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
        int B_size = matrix_properties[2] * matrix_properties[3];
        local_B = malloc(B_size * sizeof(double));
        MPI_Bcast(local_B, B_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        // allocate more space for gathering at master
        local_C = malloc(matrix_properties[0] * matrix_properties[3] * sizeof(double));
    }
    else {
        local_C = malloc(block_size * matrix_properties[3] * sizeof(double));
    }
    for (int i = 0; i < block_size; i++) {
        for (int j = 0; j < matrix_properties[3]; j++) {
            local_C[i * matrix_properties[3] + j] = 0;
            for (int k = 0; k < matrix_properties[1]; k++) {
                local_C[i * matrix_properties[3] + j] += 
                    local_A[i * matrix_properties[1] + k] * local_B[k * matrix_properties[3] + j];
            }
        }
    }
    if (rank == 0) {
        MPI_Gatherv(MPI_IN_PLACE, 0, 0, 
                    local_C, row_counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Gatherv(local_C, block_size * matrix_properties[3], MPI_DOUBLE, 
                    NULL, NULL, NULL, 0, 0, MPI_COMM_WORLD);
    }

    /** The master presents the results on the console */
    if (rank == 0){

        end_time = MPI_Wtime();
        printf("MPI Execution time: %f seconds\n", end_time - start_time);

        // matrix_struct *result_matrix = malloc(sizeof(matrix_struct));
        // result_matrix->rows = matrix_properties[0];
        // result_matrix->cols = matrix_properties[3];
        // result_matrix->mat_data = malloc(result_matrix->rows * sizeof(double));
        // for (int i = 0; i < result_matrix->rows; i++) {
        //     result_matrix->mat_data[i] = local_C + i * result_matrix->cols;
        // }
        // print_matrix(result_matrix);
        // // can't use free_matrix there because mat_data[i] points to local_C
        // free(result_matrix->mat_data); 
        // free(result_matrix);
    }

    if (rank == 0) {
        free(row_counts);
        free(displs);
    }
    free(local_A);
    free(local_B);
    free(local_C);
    
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
