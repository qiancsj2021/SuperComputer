#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define real_integral M_PI * (M_E - 2.0) / 2.0 // analytical value of the integral

double f(double x, double y, double z){
    return z * exp(x * x + y * y);
}

double generate_point(double left_range, double right_range){
    return ((double)rand() / RAND_MAX) * (right_range - left_range) + left_range; //generate point between [left_range, right_range]
}

int main(int argc, char *argv[]){

    int n_points = atoi(argv[1]); // Each process adds n_points random points at a time
    double eps = atof(argv[2]);
    int shift = atoi(argv[3]);  // = 0; = 10; = 100; = 1000; = 10000;

    double start, end;
    int my_id, num_procs;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  // Get the total number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);  // Get the rank of the current process
    start = MPI_Wtime();
    
    double sum_local = 0.0, sum = 0.0, integral = 0.0, current_gap = fabs(real_integral);

    int i, seed = my_id + shift;
    int n_times = 0;
    while (current_gap > eps){
        n_times += 1;
        srand(seed);
        seed += num_procs;

        for (i = 1; i <= n_points; ++i){
            double x = generate_point(-1.0, 1.0);
            double y = generate_point(-1.0, 1.0);
            double z = generate_point(0.0, 1.0);

            if (x * x + y * y + z * z <= 1){
                sum_local += f(x, y, z);
            }
        }
        
        MPI_Allreduce(&sum_local, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        integral = 4.0 * sum / (n_points * num_procs * n_times);
        current_gap = fabs(integral - real_integral);
    }
    end = MPI_Wtime();

    if(my_id == 0){
        printf("Finally, integral is approximated as %.8f \n", integral);
        printf("Random points: %d, the current gap: %.8f, runtime: %.4f. \n\n", n_points * num_procs * n_times, current_gap, end - start);
    }

    MPI_Finalize();
    return 0;
}
