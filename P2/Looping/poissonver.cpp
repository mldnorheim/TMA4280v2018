/**
 * Rewriting in C++ the provided C program (by Einar M. Rønquist) to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions and fast sine transforms.
 */

#include <cstdlib>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
//#include <cmath>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <string>
#include <numeric>
#include <vector>
#include <iostream>
#include <fstream>
#include <functional>

#define PI 3.14159265358979323846
#define true 1
#define false 0

using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////
/* THIS IS THE VERIFICATION TEST FOR POISSON.CPP --- RUN WITH AT LEAST 2 AND 4 PROCESSORS */
////////////////////////////////////////////////////////////////////////////////////////////


//typedef double real;
//typedef int bool;

// Function prototypes
double *mk_1D_array(size_t n, bool zero);
double **mk_2D_array(size_t n1, size_t n2, bool zero);
//void transpose(real **bt, real **b, size_t m);
double rhs(double x, double y);

extern "C" {
// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(double *v, int *n, double *w, int *nn);
void fstinv_(double *v, int *n, double *w, int *nn);
}

int main(int argc, char **argv)
{

    ofstream out_data1("PoissonErrors.txt");
    ofstream out_data2("PoissonTimings.txt");
    ofstream out_data3("Speedups.txt");
    ofstream out_data4("ParEffs.txt");


    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */

    int n;     // system size is n+1 including the boundary (n is declared in the makefile or as a part of the verification loop)
    int myrank, P, t;

    //cout << "hello, I want threads" << endl;

    /* -----------------------------------------------------------------------------------------------------*/
    /* MPI begins ------------------------------------------------------------------------------------------*/

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double time_start;

    size_t verstart = 3;
    size_t verend = 14;

    for(size_t ver=verstart; ver<verend; ver++)
    {


        //n = atoi(argv[1]); // makefile
        n = pow(2,ver);     // verification
        int m = n-1;
        double h = 1.0 / n; // step size in each direction

        if (n<P){
            cout << "P > n-1, some of the processors get zero rows!" << endl;
            return 0;
        }


        /*    Divide the domain into blocks for each processor (done by processor 0) */
        if (myrank == 0) {
            cout << "ver=" << ver << endl;
            cout << "n=" << n << endl;
            cout << "P=" << P << endl;
            time_start = MPI_Wtime();
            if (P > m)
                printf("Not possible to use all processors on this gridded domain. Please decrease P or increase n.");
        }

        int *blocksize = (int*) calloc(P,sizeof(int));   // length of blocks
        int *domsum = (int*) calloc(P,sizeof(int));   // partial sums of these lengths

        int *counts,*displs;
        counts = (int*)calloc(P,sizeof(int)); // MPI counts-vector, processor by processor, equal for send and recv
        displs = (int*)calloc(P,sizeof(int)); // MPI displacements-vector, equal for send and recv, zeros by default

        blocksize[0] = ceil(m/P);
        domsum[0] = 0;

        int sum = 0;
        int Sum = 0;

        displs[0] = 0;

        #pragma omp parallel for
        for (size_t i=1; i<P; i++){
            sum += blocksize[i-1];
            blocksize[i] = ceil((m-sum)/(P-i));
            domsum[i] = sum;
        }
        counts[0] = blocksize[myrank]*blocksize[0];
        #pragma omp parallel for
        for(size_t l=1; l<P; l++){
            counts[l] = blocksize[myrank]*blocksize[l];
            Sum += counts[l-1];
            displs[l] = Sum;
        }

       // if (myrank==0){
       //     int bstot = accumulate(blocksize.begin(), blocksize.end(), 0);
       //     if ( bstot!= m)
       //         cout << "The elements of blocksize do not add up to m!" << endl;
       // }



    //	  MPI_Barrier(MPI_COMM_WORLD);
    //    MPI_Bcast(&blocksize, 1, MPI_INT, 0, MPI_COMM_WORLD); // sharing blocksize with the other processors
    //    MPI_Bcast(&domsum, 1, MPI_INT, 0, MPI_COMM_WORLD);    // sharing domsum with the other processors



        /* ------------------------------------------------------------------------------
         * Grid points are generated with constant mesh size on both x- and y-axis.
         */
        double *grid = mk_1D_array(m, false);
        for (size_t i = 0; i < m; i++){
            grid[i] = i * h; // As the boundary grid-points aren't used in the code anyway, they won't be calculated
        }
        //     We let processor P-1 check the partition:
        if (myrank == P-1) {
            if (domsum[myrank]+blocksize[myrank] < m)
                cout << "Not all of the domain is covered by the blocks! Check domsum and blocksize." << endl;
            else if (domsum[myrank]+blocksize[myrank] > m)
                cout << "The blocks exceed the domain! Check domsum and blocksize." << endl;
        }

        /* Exact solution for the given RHS */
        double verexact [blocksize[myrank]][m];
        #pragma omp parallel
        {
            t = omp_get_num_threads();
            #pragma omp for collapse(2)
            for (size_t i = 0; i < blocksize[myrank]; i++){
                for (size_t j=0; j<m; j++){
                    verexact[i][j] = sin(PI*grid[domsum[myrank]+i]) *sin(2*PI*grid[j]);
                }
            }
        }
	//cout << "t=" << t << endl;

        /* ------------------------------------------------------------------------------
         * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
         * defined Chapter 9. page 93 of the Lecture Notes.
         * Note that the indexing starts from zero here, thus i+1.
         */
        double *diag = mk_1D_array(m, false);
        #pragma omp parallel for
        for (size_t i = 0; i < m; i++) {
            diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
        }


        /* ------------------------------------------------------------------------------
         * Allocate the matrices b and bt which will be used for storing value of
         * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
         */
        double **b = mk_2D_array(blocksize[myrank], m, false);
        double **bt = mk_2D_array(blocksize[myrank], m, false);
      //  double **Error = mk_2D_array(m, m, false);


        /* ------------------------------------------------------------------------------
         * Initialize the right hand side data for a given rhs function.
         * Note that the right hand-side is set at nodes corresponding to degrees
         * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
         *
         */
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < blocksize[myrank]; i++) {
            for (size_t j = 0; j < m; j++) {
                b[i][j] = h * h * rhs(domsum[myrank]+grid[i], grid[j]); // grid boundary excluded
            }
        }


        /* ------------------------------------------------------------------------------
         * This vector will holds coefficients of the Discrete Sine Transform (DST)
         * but also of the Fast Fourier Transform used in the FORTRAN code.
         * The storage size is set to nn = 4 * n, look at Chapter 9. pages 98-100:
         * - Fourier coefficients are complex so storage is used for the real part
         *   and the imaginary part.
         * - Fourier coefficients are defined for j = [[ - (n-1), + (n-1) ]] while
         *   DST coefficients are defined for j [[ 0, n-1 ]].
         * As explained in the Lecture notes coefficients for positive j are stored
         * first.
         * The array is allocated once and passed as arguments to avoid doings
         * reallocations at each function call.
         */
        int nn = 4 * n;

        //double *blockvectorb, *bvpretransp, *blockvectorbt;
        //blockvectorb = (double*)malloc(m*blocksize[myrank]*sizeof(double));
        //bvpretransp = (double*)malloc(m*blocksize[myrank]*sizeof(double));  // these vectors are used for MPI_Alltoallv
        //blockvectorbt = (double*)malloc(m*blocksize[myrank]*sizeof(double));
        double *blockvectorb = mk_1D_array(blocksize[myrank]*m,false);
        double *bvpretransp = mk_1D_array(blocksize[myrank]*m,false);
        double *blockvectorbt = mk_1D_array(blocksize[myrank]*m,false);

        double u_max = 0.0;
        double errmax = 0.0;

        #pragma omp parallel
        {

            double *z = mk_1D_array(nn, false);

            /* ------------------------------------------------------------------------------
             * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
             * Instead of using two matrix-matrix products the Discrete Sine Transform
             * (DST) is used.
             * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
             * The array zz is used as storage for DST coefficients and internally for
             * FFT coefficients in fst_ and fstinv_.
             * In functions fst_ and fst_inv_ coefficients are written back to the input
             * array (first argument) so that the initial values are overwritten.
             */

            // Using fst on b:
            for (size_t i = 0; i < blocksize[myrank]; i++) {
                fst_(b[i], &n, z, &nn);
            }

            // Shuffling matrix b (stored row-by-row) into a vector of blocks (subrow-by-subrow):
            #pragma omp for collapse(2)
            for(size_t i = 0; i < blocksize[myrank]; i++){
                for(size_t k=0; k<P; k++){
                    for(size_t j=domsum[k]; j<(domsum[k]+blocksize[k]); j++){
                        blockvectorb[ displs[k]+i*blocksize[k]+j ] = b[i][j];
                    }
                }
            }

            // Exchanging data (before transposing):
            #pragma omp master
            {
                MPI_Alltoallv(blockvectorb, counts, displs, MPI_DOUBLE, bvpretransp, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
            }

            // Transposing, block-wise:
            for(size_t k=0; k<P; k++){
                int M = counts[k]/blocksize[k];          // #rows of original block
                #pragma omp for collapse(2)
                for(size_t i=0; i<blocksize[k]; i++){
                    for(size_t j=0; j<M; j++){
                        blockvectorbt[displs[k] + blocksize[k]*j + i] = bvpretransp[displs[k] + M*i + j];
                    }
                }
            }


            // Shuffling back, from vector to matrix bt:
            #pragma omp for collapse(2)
            for(size_t i = 0; i < blocksize[myrank]; i++){
                for(size_t k=0; k<P; k++){
                    for(size_t j=domsum[k]; j<(domsum[k]+blocksize[k]); j++){
                         bt[i][j] = blockvectorbt[ displs[k]+i*blocksize[k]+j ];
                    }
                }

            }

            // Using fstinv on bt:
            for (size_t i = 0; i < blocksize[myrank]; i++) {
                fstinv_(bt[i], &n, z, &nn);
            }


            /* ------------------------------------------------------------------------------
             * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
             */
            #pragma omp for collapse(2)
            for (size_t i = 0; i < blocksize[myrank]; i++) {  // or switch j and i
                for (size_t j = 0; j < m; j++) {
                    if((diag[i] + diag[j]) == 0) perror("zero diag");
                        bt[i][j] = bt[i][j] / (diag[domsum[myrank]+i] + diag[j]);
                }
            }


            /* ------------------------------------------------------------------------------
             * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
             */


            // Using fst on bt:
            for (size_t i = 0; i < blocksize[myrank]; i++) {
                fst_(bt[i], &n, z, &nn);
            }

            // Shuffling matrix bt (stored row-by-row) into a vector of blocks (subrow-by-subrow):
            #pragma omp for collapse(2)
            for(size_t i = 0; i < blocksize[myrank]; i++){
                for(size_t k=0; k<P; k++){
                    for(size_t j=domsum[k]; j<(domsum[k]+blocksize[k]); j++){
                        blockvectorbt[ displs[k]+i*blocksize[k]+j ] = bt[i][j];
                    }
                }

            }

            // Exchanging data (before transposing):
            #pragma omp master
            {
                MPI_Alltoallv(blockvectorbt, counts, displs, MPI_DOUBLE, bvpretransp, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
            }

            // Transposing, block-wise:
            for(size_t k=0; k<P; k++){
                int M = counts[k]/blocksize[k];          // #rows of original block
                #pragma omp for collapse(2)
                for(size_t i=0; i<blocksize[k]; i++){
                    for(size_t j=0; j<M; j++){
                        blockvectorb[displs[k] + blocksize[k]*j + i] = bvpretransp[displs[k] + M*i + j];
                    }
                }
            }


            // Shuffling back, from vector to matrix b:
            #pragma omp for collapse(2)
            for(size_t i = 0; i < blocksize[myrank]; i++){
                for(size_t k=0; k<P; k++){
                    for(size_t j=domsum[k]; j<(domsum[k]+blocksize[k]); j++){
                         b[i][j] = blockvectorb[ displs[k]+i*blocksize[k]+j ];
                    }
                }

            }

            // Using fstinv on b:
            for (size_t i = 0; i < blocksize[myrank]; i++) {
                fstinv_(b[i], &n, z, &nn);
            }

            /* ------------------------------------------------------------------------------
             * Compute maximal value of solution for convergence analysis in L_\infty
             * norm and the largest difference in absolute value between numerical and exact solution.
             */


            double u_maxall, errmaxall;
            #pragma omp for reduction(max:u_max,errmax) collapse(2)
            for (size_t i = 0; i < blocksize[myrank]; i++) {
                for (size_t j = 0; j < m; j++) {
                    if (u_max < b[i][j])
                        u_max = b[i][j];
                    if(errmax < fabs(b[i][j]-verexact[i][j])){
                        errmax = fabs(b[i][j]-verexact[i][j]);
                    }
                }
            }

            MPI_Reduce(&u_max, &u_maxall, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&errmax, &errmaxall, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


            if (myrank == 0) {
                double timing = MPI_Wtime() - time_start;
                double *T_1 = mk_1D_array(verend-verstart,false);
                #pragma omp master
                {
                    // Serial timings for different values of n, computed with this program for P=1 and no OpenMP:
                    T_1[0] = 0.0001297695562243461609;
                    T_1[1] = 0.0004363982006907463074;
                    T_1[2] = 0.001938818953931331635;
                    T_1[3] = 0.008492053486406803131;
                    T_1[4] = 0.03761203121393918991;
                    T_1[5] = 0.1575700864195823669;
                    T_1[6] = 0.6742947511374950409;
                    T_1[7] = 2.886859129182994366;
                    T_1[8] = 28.86859129182994366;
                    T_1[9] = 52.16260548308491707;
                    T_1[10] = 220.8320795409381390;
                    if (P == 1){
                        double T_1 = timing;
                        cout << "For 1 processor the timing is T1 = " << scientific << setprecision(18) << T_1 << endl;
                    }
                    double timing = MPI_Wtime() - time_start;  // = T_P
                    printf("For loop number %i: u_max = %e and maximum error = %e\n", ver, u_maxall, errmaxall);
                    out_data1 << log(h) << " " << scientific << setprecision(18) << log(errmaxall) << endl;
                    out_data2 << log(h) << " " << scientific << setprecision(18) << log(timing) << endl;
                    double speedup = timing/T_1[ver-verstart];
                    double pareff = speedup/P;  // parallel efficiency
                    out_data3 << n << " " << P << " " << scientific << setprecision(18) << speedup << endl;
                    out_data4 << n << " " << P << " " << scientific << setprecision(18) << pareff << endl;
                }
            }


        } // end #pragma

        // Freeing memory:
        for (size_t i=0; i<blocksize[myrank]; i++){
            free(bt[i]);
            free(b[i]);
        }
        free(displs);
        free(counts);
        free(blocksize);
        free(domsum);
        free(blockvectorb);
        free(blockvectorbt);
        free(bvpretransp);

    } // end verification for-loop


    /* MPI ends --------------------------------------------------------------------------------------------*/
    /* -----------------------------------------------------------------------------------------------------*/
    MPI_Finalize();


    return 0;
} // end function


/* ------------------------------------------------------------------------------
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

double rhs(double x, double y) {
    return 5*PI*PI*sin(PI*x)*sin(2*PI*y);
}

/* ------------------------------------------------------------------------------
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

//MPI_Init(&argc, &argv); //---------------------------------------------------------------------------
//void transpose(real **bt, real **b, size_t m)
//{
//    for (size_t i = domsum[myrank]+1; i <= (domsum[myrank]+blocksize[myrank]); i++) {
//        #pragma omp parallel for
//        for (size_t j = 1; j <= m; j++) {
//            bt[i-1][j-1] = b[j-1][i-1];
//        }
//    }
//}
//MPI_Finalize(); //-----------------------------------------------------------------------------------

/* ------------------------------------------------------------------------------
 * The allocation of a vector of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

double *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (double *)calloc(n, sizeof(double));
    }
    return (double *)malloc(n * sizeof(double));
}

/* ------------------------------------------------------------------------------
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contiguous,
 * 3. pointers are set for each row to the address of first element.
 */

double **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    double **ret = (double **)malloc(n1 * sizeof(double *));

    // 2
    if (zero) {
        ret[0] = (double *)calloc(n1 * n2, sizeof(double));
    }
    else {
        ret[0] = (double *)malloc(n1 * n2 * sizeof(double));
    }

    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}
