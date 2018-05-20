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
    if (argc < 2) {
        printf("Usage:\n");
        printf("  poisson n\n\n");
        printf("Arguments:\n");
        printf("  n: the problem size (must be a power of 2)\n");
    }

    /* Listing the speedups and parallel efficiencies for different values of n and P (P*t=36?) */
  //  ofstream out_data1("Speedups.txt");
  //  ofstream out_data2("ParEffs.txt");

    /*
     *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
     *  conditions are applied on the boundary:
     *  - the number of grid points in each direction is n+1,
     *  - the number of degrees of freedom in each direction is m = n-1,
     *  - the mesh size is constant h = 1/n.
     */

    int n;     // system size is n+1 including the boundary (n is declared in the makefile or as a part of the verification loop)
    int ex = atoi(argv[1]); // makefile
    n = pow(2,ex);

    cout << "Test1: before MPI is initialized" << endl; /////////////////////////////////////


    int m = n-1;
    double h = 1.0 / n;
    int myrank, P, t;


    // Unit test (for n=5)
    double **unitexact = mk_2D_array(m, m, false);
    if(n == 5) {
        unitexact[0][0], unitexact[0][3], unitexact[3][0], unitexact[3][3] = 1/30;
        unitexact[1][1], unitexact[1][2], unitexact[2][1], unitexact[2][2] = 1/15;
        unitexact[0][1], unitexact[0][2], unitexact[3][1], unitexact[3][2], unitexact[1][0], unitexact[2][0], unitexact[1][3], unitexact[2][3] = 7/150;
    }


    /* -----------------------------------------------------------------------------------------------------*/
    /* MPI begins ------------------------------------------------------------------------------------------*/

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (n<P){
	    cout << "P > n-1, some of the processors get zero rows!" << endl;
	    return 0;
	}
	if (myrank==0){
   	    cout << "Test2: MPI-start" << endl; /////////////////////////////////////
    }

    double time_start;

    /*    Divide the domain into blocks for each processor (done by processor 0) */
    if (myrank == 0) {
        time_start = MPI_Wtime();
        if (P > m)
            printf("Not possible to use all processors on this gridded domain. Please decrease P or increase n.");
    }

    //vector <int> blocksize (P);
    //vector <int> domsum (P,0);
	int *blocksize = (int*) calloc(P,sizeof(int));  // length of blocks
	int *domsum = (int*) calloc(P,sizeof(int));    // partial sums of these lengths

    //int *counts,*displs;
    //counts = (int*)malloc(P*sizeof(int));
    //displs = (int*)malloc(P*sizeof(int));
	int *counts = (int*) calloc(P,sizeof(int)); // MPI counts-vector, processor by processor, equal for send and recv
	int *displs = (int*) calloc(P,sizeof(int)); // MPI displacements-vector, equal for send and recv, zeros by default

	//cout << "is malloc ok? my rank is " << myrank << endl;
    blocksize[0] = ceil(m/P);
    domsum[0] = 0;

    int sum = 0;
    int Sum = 0;

    #pragma omp parallel for
    for (size_t i=1; i<P; i++){
        sum += blocksize[i-1];
        blocksize[i] = ceil((m-sum)/(P-i));
        domsum[i] = sum;
    }
    displs[0] = 0;
    counts[0] = blocksize[myrank]*blocksize[0];
    #pragma omp parallel for
    for(size_t l=1; l<P; l++){
        Sum += counts[l-1];
        counts[l] = blocksize[myrank]*blocksize[l];
        displs[l] = Sum;
    }

    if (myrank==0){
        cout << "Test3: After blocksize, domsum etc" << endl; /////////////////////////////////////
    //    int bstot = accumulate(blocksize.begin(), blocksize.end(), 0);
    //    if ( bstot != m)
    //        cout << "The elements of blocksize do not add up to m!" << endl;
    }



    /* ------------------------------------------------------------------------------
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    double *grid = mk_1D_array(m, false);
	#pragma omp parallel for
	for (size_t i = 0; i < m; i++){
        grid[i] = i * h; // As the boundary grid-points aren't used in the code anyway, they won't be calculated
    }

    //     We let processor P-1 check the partition
    if (myrank == P-1) {
        if (domsum[myrank]+blocksize[myrank] < m)
            cout << "Not all of the domain is covered by the blocks! Check domsum and blocksize." << endl;
        else if (domsum[myrank]+blocksize[myrank] > m)
            cout << "The blocks exceed the domain! Check domsum and blocksize." << endl;
    }


    /* ------------------------------------------------------------------------------
     * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
     * defined Chapter 9. page 93 of the Lecture Notes.
     * Note that the indexing starts from zero here, thus i+1.
     */
    double *diag = mk_1D_array(m, false);
	#pragma omp parallel for
    for (size_t i = 0; i <m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }


    /* ------------------------------------------------------------------------------
     * Allocate the matrices b and bt which will be used for storing value of
     * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
     */
    double **b = mk_2D_array(blocksize[myrank], m, false);
    double **bt = mk_2D_array(blocksize[myrank], m, false);
    double **difference = mk_2D_array(blocksize[myrank], m, false);

	if (myrank==0){
	    cout << "Test4: after declaring matrices and vectors" << endl; /////////////////////////////////////
    }

    /* ------------------------------------------------------------------------------
     * Initialize the right hand side data for a given rhs function.
     * Note that the right hand-side is set at nodes corresponding to degrees
     * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
     *
     */
	#pragma omp parallel
	{
	    t = omp_get_num_threads();
	    #pragma omp for collapse(2)
        for (size_t i = 0; i < blocksize[myrank]; i++) {
            for (size_t j = 0; j < m; j++) {
                b[i][j] = h * h * rhs(grid[domsum[myrank]+i], grid[j]); // grid boundary excluded
            }
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
    //bvpretransp = (double*)malloc(m*blocksize[myrank]*sizeof(double));
    //blockvectorbt = (double*)malloc(m*blocksize[myrank]*sizeof(double));
	double *blockvectorb = mk_1D_array(blocksize[myrank]*m,false);
	double *bvpretransp = mk_1D_array(blocksize[myrank]*m,false);
	double *blockvectorbt = mk_1D_array(blocksize[myrank]*m,false);

	if (myrank==0){
	    cout << "Test5: before pragma-block" << endl; /////////////////////////////////////
    }


    double u_max = 0.0;
    double diffmax = 0.0;

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
        //delete [] blockvectorb;

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

        if (myrank==0){
            cout << "Test6: in pragma-block, after step1 in algorithm" << endl; /////////////////////////////////////
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

        if (myrank==0){
            cout << "Test7: in pragma-block, after step2 in algorithm" << endl; /////////////////////////////////////
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
                }   // maybe i-->(i-domsum[k]) in bvb, if b-values come from global matrix not local
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

        if (myrank==0){
            cout << "Test8: in pragma-block, after step3 in algorithm" << endl; /////////////////////////////////////
        }

        /* ------------------------------------------------------------------------------
         * Compute maximal value of solution for convergence analysis in L_\infty
         * norm.
         */
       // double u_max = 0.0;
        double u_maxall;
        #pragma omp for reduction(max:u_max) collapse(2)
        for (size_t i = 0; i < blocksize[myrank]; i++) {
            for (size_t j = 0; j < m; j++) {
                if (u_max < b[i][j])
                    u_max = b[i][j];
                if (n==5){
                    difference[i][j] = fabs( unitexact[i+domsum[myrank]][j] - b[i][j] );
                    if (diffmax < difference[i][j])
                        diffmax = difference[i][j];
                }
            }
        }
        MPI_Reduce(&u_max, &u_maxall, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (myrank == 0) {
            #pragma omp master
            {
                double timing = MPI_Wtime() - time_start;  // = T_P
                printf("For n=%i, P=%i and t=%i we get:\n u_max=%e\n timing=%e\n greatest error=%e\n", n, P, t, u_maxall, timing, diffmax );
                //cout << "For n = " << n << ", P = " << P << " and t = " << t << " we get:" << endl;
                //printf("u_max = %e\n", u_maxall);
                //printf("timing = %e\n", timing);
            }
            cout << "Test9: in pragma-block, after computing and printing max.values and timings" << endl; /////////////////////////////////////
        }

        //free(z);
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
	free(difference);





    /* MPI ends --------------------------------------------------------------------------------------------*/
    /* -----------------------------------------------------------------------------------------------------*/
    MPI_Finalize();

	if (myrank==0){
	    cout << "Test10: after pragma-block and MPI_Finalize" << endl; /////////////////////////////////////
    }

    return 0;

}


/* ------------------------------------------------------------------------------
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

double rhs(double x, double y) {
    return 5*PI*PI*sin(PI*x)*sin(2*PI*y); //1;
}

/* ------------------------------------------------------------------------------
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

/* ------------------------------------------------------------------------------
 * The allocation of a vectore of size n is done with just allocating an array.
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
 *   is contigusous,
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
