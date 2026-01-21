#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <mpi.h>
#include <unistd.h>

double my_elapsed;
struct timespec my_start, my_finish;
double local_time_start, local_finish, local_elapsed;
double elapsed;
double local_total_start;

void generate_main_array(int **main_array, int n, long long num_zeros){

    clock_gettime(CLOCK_REALTIME, &my_start );

    long long num_nonzeros = (n * n) - num_zeros;
    int i, j;

    for(int k = 0; k < num_nonzeros; k++){
        i = rand() % n;
        j = rand() % n;
        while(main_array[i][j] != 0){
            
            if(j < n - 1){
                j++;
            }

            else if(i < n - 1){
                j = 0;
                i ++;
            }

            else{
                j = 0;
                i = 0;
            }

        }
        main_array[i][j] = rand() + 1;
        main_array[i][j]*= (rand() % 2 ? 1 : -1);
            
    }    

    clock_gettime( CLOCK_REALTIME, &my_finish );    //finish time of generating main array

    my_elapsed = (my_finish.tv_sec - my_start.tv_sec) + (my_finish.tv_nsec - my_start.tv_nsec) / 1e9;

    printf("Creation of main array time = %lf seconds\n", my_elapsed);

}

void generate_vector(int n, double *dian){
    for(int i = 0; i < n; i++){
        dian[i] = rand();
        dian[i]*= (rand() % 2 ? 1 : -1);
    }
}

void generate_csr_arrays(int n, int **main_array, int *v, int *row_index, int *col_index, int num_nonzeros){

    
    int *row_nonz;

    row_nonz = calloc(n, sizeof(int));

    row_index[0] = 0;

    clock_gettime(CLOCK_REALTIME, &my_start );             //start  time of generating csr


    for(int i = 0; i < n; i++){
        int counter = 0;                        //how many zeros per row
        for(int j = 0; j < n; j++){
            if(main_array[i][j] != 0){
                counter ++;
            }
        }
        row_nonz[i] = counter;
    }

    for(int i = 0; i < n; i++){
        row_index[i+1] = row_index[i] + row_nonz[i];        //fill row_index array
    }




    for(int i = 0; i < n; i++){
        int pos = row_index[i];
        for(int j = 0; j < n; j++){
            if(main_array[i][j] != 0){
                v[pos] = main_array[i][j];
                col_index[pos] = j;
                pos++;
            }
        }
    }


    clock_gettime( CLOCK_REALTIME, &my_finish );    //finish time of generating csr

    my_elapsed = (my_finish.tv_sec - my_start.tv_sec) + (my_finish.tv_nsec - my_start.tv_nsec) / 1e9;

    printf("Creation of csr time = %lf seconds\n", my_elapsed);

    free(row_nonz);

}

double* dense_mult(int n, int *local_array, double *dian, int my_rank, int * row_count){

    double *result = calloc(row_count[my_rank], sizeof(double));
    if (!result) {
        fprintf(stderr, "calloc failed\n"); 
        exit(EXIT_FAILURE);
    }

    //densen mult
    for(int i = 0; i < row_count[my_rank]; i ++){
        for(int j = 0; j < n; j++){
            result[i] += local_array[i*n+j] * dian[j]; 
        }
    }
    return result;

}

double* csr_mult(int *local_v, int *local_col_index, int *local_row_index, int local_rows, double *dian){
    
    double *result = calloc(local_rows, sizeof(double));
    if (!result) {
        fprintf(stderr, "calloc failed\n"); 
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < local_rows; i++) {
        for (int k = local_row_index[i]; k < local_row_index[i+1]; k++){
            result[i] += local_v[k] * dian[local_col_index[k]]; 
        }
    }

    return result;
}


int compare_results(int n, double *dian_a, double * dian_b){

    for(int i = 0; i < n; i++){
        if (fabs(dian_a[i] - dian_b[i]) > 1e-6){
            return 0;
        }
    }
    return 1;
}

double* serial_mult(int n, int **main_array, double *dian) {

    // Δέσμευση μνήμης για το αποτέλεσμα και αρχικοποίηση σε μηδέν με την calloc
    double *result = calloc(n, sizeof(double));
    if (!result) {
        fprintf(stderr, "calloc failed\n"); 
        exit(EXIT_FAILURE);
    }

    // Σειριακή εκτέλεση των βρόχων (loops)
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){ 
            // Υπολογισμός του γινομένου για κάθε γραμμή i [cite: 683, 786, 868]
            result[i] += (double)main_array[i][j] * dian[j];
        }
    }

    return result;
}


int main(int argc, char *argv[]){

    int n, pct_zeros, mult_count, rem, local_count, *row_count, *starting_row;
    int **main_array = NULL, *v = NULL, *col_index = NULL, *row_index = NULL, *local_array, local_nz, *local_v, *local_col_index, *local_row_index = NULL;
    double *dian, *local_dian;
    int *main_data = NULL   ;
    long long num_zeros, num_nonzeros;
    double *dian_a = NULL, *dian_b = NULL, *dian_c = NULL, *temp = NULL;
    int comm_sz;                //number of proccesses
    int my_rank;                //rank of current proccess
    int *nz_displ, *nz_count;


    n = strtol(argv[1], NULL, 10);                      //array dimensions

    pct_zeros = strtol(argv[2], NULL, 10);              //percentage of zeros
    
    mult_count = strtol(argv[3], NULL, 10);               //multiplication count

    num_zeros = (long long) llround((double)n * n * pct_zeros / 100.0);                //find number of zeros needed
    
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
                 

    if(my_rank == 0){                                //if process 0
        main_data = calloc(n * n, sizeof(int));                    //flat version of main_array
        if (!main_data) { perror("malloc main_data"); return 1; }        //main_array stores all elements contiguously 
        main_array = malloc(n * sizeof(int*));

        for (int i = 0; i < n; i++){
            main_array[i] = &main_data[i*n];
        }

        generate_main_array(main_array, n, num_zeros);                

        dian = calloc(n, sizeof(double));               //allocate and initialize dian

        generate_vector(n, dian);

        num_nonzeros = (n * n) - num_zeros;

        v = malloc(num_nonzeros * sizeof(int));
        col_index = malloc(num_nonzeros * sizeof(int));
        row_index = malloc((n + 1) * sizeof(int));

        generate_csr_arrays(n, main_array, v, row_index, col_index, num_nonzeros);

         dian_a = calloc(n, sizeof(double));
        if (!dian_a) {
            fprintf(stderr,"malloc dian_b failed\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }


        clock_gettime( CLOCK_REALTIME, &my_start );     //start time of Dense mult
    }


    if (my_rank != 0) {
        dian = malloc(n * sizeof(double));
        if (!dian) { fprintf(stderr, "malloc failed for dian on rank %d\n", my_rank); MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);                            //wait for every proccess
    local_time_start = MPI_Wtime();

    MPI_Bcast(dian, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);          //send dian to every process


    row_count = malloc(comm_sz * sizeof(int));          //number of rows per row
    starting_row = malloc(comm_sz * sizeof(int));       //starting row for every process
    if (!row_count || !starting_row) {
        fprintf(stderr,"malloc failed\n"); MPI_Abort(MPI_COMM_WORLD,1); 
    }

    local_count = n / comm_sz;
    rem = n % comm_sz;


    for(int i = 0; i < comm_sz; i ++){
        row_count[i] = local_count + (i < rem ? 1 : 0);         //calculate number of rows of main_array for every process
    }

    starting_row[0] = 0;
    for(int i = 1; i < comm_sz; i ++){
        starting_row[i] = starting_row[i-1] + row_count[i-1];       //calculate first row of main_array for every process
    }

    int *sendcounts = malloc(comm_sz * sizeof(int));
    int *displs     = malloc(comm_sz * sizeof(int));
    
    for (int i = 0; i < comm_sz; i++) {
        sendcounts[i] = row_count[i] * n;   //number of ints that will be send to each process
    }

    displs[0] = 0;
    for (int p = 1; p < comm_sz; ++p) {
        displs[p] = displs[p-1] + sendcounts[p-1];
    }

    int local_rows = row_count[my_rank];
    int recvcount = local_rows * n;             //number of ints this process will receive

    local_array = malloc(recvcount * sizeof(int)); 


    MPI_Scatterv(
    (my_rank == 0 ? main_data : NULL),
    sendcounts,
    displs,
    MPI_INT,
    local_array,
    recvcount,
    MPI_INT,
    0,
    MPI_COMM_WORLD
    );

    local_finish = MPI_Wtime();
    local_elapsed = local_finish - local_time_start;
    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);        //find max time
    
    if(my_rank == 0)
        printf("MPI_bcast and scatter time = %lf\n", elapsed);                  //print time of sending data to other processes


    local_dian = malloc(recvcount * sizeof(double));
    if(!local_dian){ 
            fprintf(stderr, "malloc failed for local_dian on rank %d\n", my_rank); MPI_Abort(MPI_COMM_WORLD, 1);
    }
    temp = malloc(recvcount * sizeof(double));


    MPI_Barrier(MPI_COMM_WORLD);                            //wait for every proccess
    local_time_start = MPI_Wtime();
    local_dian = dense_mult(n, local_array, dian, my_rank, row_count);


    if(my_rank != 0) {
        dian_a = calloc(n, sizeof(double));
        if(!dian_a){ 
            fprintf(stderr, "malloc failed for dian_a on rank %d\n", my_rank); MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Gatherv(
        local_dian,                      
        local_rows,                   
        MPI_DOUBLE,
        dian_a,                       // receive buffer 
        row_count,                    // counts 
        starting_row,                 // displacements 
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );



    if (my_rank != 0) {
        row_index = malloc((n + 1) * sizeof(int));      //allocate row_index for every process != 0x    
    }

    MPI_Bcast(dian_a, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);                    
    
    for(int i = 1; i < mult_count; i++){

        temp = dense_mult(n, local_array, dian_a, my_rank, row_count);          //complete mult mult_count-1 times
        free(local_dian);
        local_dian = temp;
        MPI_Gatherv(
        local_dian,                    
        local_rows,                    
        MPI_DOUBLE,
        dian_a,                       // receive buffer 
        row_count,                    // counts 
        starting_row,                 // displacements 
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD);
        MPI_Bcast(dian_a, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);                    
        

    }
    local_finish = MPI_Wtime();
    local_elapsed = local_finish - local_time_start;
    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);        //find max time
    
    if(my_rank == 0)
        printf("Dense mult time time = %lf\n", elapsed);                  //print time of dense mult

    MPI_Barrier(MPI_COMM_WORLD);                            //wait for every proccess
    local_total_start = MPI_Wtime();                //start of total csr time
    MPI_Barrier(MPI_COMM_WORLD);

    //broadcast row_index
    MPI_Bcast(row_index, n + 1, MPI_INT, 0, MPI_COMM_WORLD);  
    nz_count = malloc(comm_sz * sizeof(int));                                //number of nonzeros per process 
    for(int i = 0; i < comm_sz; i++){
        nz_count[i] = row_index[starting_row[i] + row_count[i]] - row_index[starting_row[i]];    
    }
    local_nz = nz_count[my_rank];


    local_v = malloc(nz_count[my_rank] * sizeof(int));
    local_col_index = malloc(nz_count[my_rank] * sizeof(int));
    local_row_index = malloc((row_count[my_rank] + 1) * sizeof(int));
    if(!local_v || !local_col_index || !local_row_index){ 
            fprintf(stderr, "malloc failed on rank %d\n", my_rank); MPI_Abort(MPI_COMM_WORLD, 1);
    }
    else{
        local_v = NULL;
        local_col_index = NULL;
    }

    nz_displ = malloc(comm_sz * sizeof(int));
    nz_displ[0] = 0;                                    //number of nonzeros per process
    for (int i = 1; i < comm_sz; i++) {
        nz_displ[i] = nz_displ[i-1] + nz_count[i-1];
    }


    //scatter csr to every process
 
    //scatter v
    MPI_Scatterv(my_rank==0? v:NULL, nz_count, nz_displ, MPI_INT, local_v, local_nz, MPI_INT, 0, MPI_COMM_WORLD);
    
    //scatter col_index
    MPI_Scatterv(my_rank==0? col_index:NULL, nz_count, nz_displ, MPI_INT, local_col_index, local_nz, MPI_INT, 0, MPI_COMM_WORLD);
    
    for (int i = 0; i <= row_count[my_rank]; ++i) {
        local_row_index[i] = row_index[starting_row[my_rank] + i] - row_index[starting_row[my_rank]];
    }


    //csr mult

    MPI_Barrier(MPI_COMM_WORLD);                            //wait for every proccess
    local_time_start = MPI_Wtime();
    local_dian = csr_mult(local_v, local_col_index, local_row_index, row_count[my_rank], dian);
    
    if (my_rank == 0) {
        dian_b = malloc(n * sizeof(double));
    }


    if(my_rank != 0) {
        dian_b = malloc(n * sizeof(double));
        if(!dian_b){ 
            fprintf(stderr, "malloc failed for dian_b on rank %d\n", my_rank); MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Gatherv(
    local_dian,                      
    local_rows,                   
    MPI_DOUBLE,
    dian_b,                       // receive buffer 
    row_count,                    // counts 
    starting_row,                 // displacements 
    MPI_DOUBLE,
    0,
    MPI_COMM_WORLD
    );



    MPI_Bcast(dian_b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);                    
    
    for(int i = 1; i < mult_count; i++){

        temp = csr_mult(local_v, local_col_index, local_row_index, row_count[my_rank], dian_b);         //complete mult mult_count-1 times
        free(local_dian);
        local_dian = temp;
        MPI_Gatherv(
        local_dian,                    
        local_rows,                    
        MPI_DOUBLE,
        dian_b,                       // receive buffer 
        row_count,                    // counts 
        starting_row,                 // displacements 
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD);
        MPI_Bcast(dian_b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);                    
    }
        local_finish = MPI_Wtime();
    local_elapsed = local_finish - local_time_start;
    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);        //find max time
    
    if(my_rank == 0)
        printf("CSR mult time = %lf\n", elapsed);                  //print time of csr mult

    MPI_Barrier(MPI_COMM_WORLD);
    double total_finish = MPI_Wtime();
    double total_local = total_finish - local_total_start;
    MPI_Reduce(&total_local, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (my_rank==0) printf("Total csr time = %lf\n", elapsed + my_elapsed);         //my_elapsed = csr creation time
    
    if(my_rank == 0){


    //compare results of both algorithms
        if(compare_results(n, dian_a, dian_b)){
            printf("Success: both algorithms gave the same result!\n");
        }

        else{
            printf("Failure: Algorithms gave different results\n");
        }

        //serial mult for comparison
        dian_c = malloc(n * sizeof(int));
        dian_c = serial_mult (n, main_array, dian);
        for(int i = 1; i < n; i++){
            temp = serial_mult(n, main_array, dian_c);
            free(dian_c);
            dian_c = temp;
        }

        //compare results of dense mult with serial mult
        if(compare_results(n, dian_a, dian_b)){
            printf("Success: both parallel dense and serial algorithms gave the same result!\n");
        }

        else{
            printf("Failure: Algorithms gave different results\n");
        }
    }




    free(dian);

    free(dian_a);
    free(dian_b);
    free(dian_c);

    free(row_index);
    if(my_rank == 0){
        free(col_index);
        free(v);
        free(main_data);
    }

    free(local_array);
    free(local_col_index);
    free(local_dian);
    free(local_row_index);
    free(row_count);
    free(nz_count);
    free(displs);
    free(nz_displ);    

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    

    return 0;
}

    
 
