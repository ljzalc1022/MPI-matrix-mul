#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "matrix.h"

int main(int argc, char *argv[]) {
    /** Matrix Properties
     * [0] = Rows of Matrix A
     * [1] = Cols of Matrix A
     * [2] = Rows of Matrix B
     * [3] = Cols of Matrix B
     **/
    int matrix_properties[4];
    
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

        // TODO

        free_matrix(m_1);
        free_matrix(m_2);
    }

    // TODO

    /** The master presents the results on the console */
    if (rank == 0){

        end_time = MPI_Wtime();
        printf("MPI Execution time: %f seconds\n", end_time - start_time);

    }
    
    
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
