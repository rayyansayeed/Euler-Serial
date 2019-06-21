#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Accelerate/Accelerate.h>

/* ADI method to solve linearized Euler equations.
Input args: N, dimension of square grid
            M, # of timesteps */

int main(int argc, char* argv[]){

	// Start runtime clock

	clock_t start = clock();

	/* Defining constants here: gamma, terminal time T,
	# number of timesteps M, grid dimension N, dt=T/M, dx=2/N */

	double gam = 1.4;
	int T = 2;
	int M;
	int N;
	double dt = (double) T/M;
	double dx = (double) 2/N;

	// Set up M and N to be taken as inputs

	if (argc >1){
		N = atoi(argv[1]);
		M = atoi(argv[2]);
	}

    // Initialize matrices used by ADI

    // (N x N) matrices

    	double (*rho)[N] = malloc(N*sizeof(*rho));
    	double (*u)[N] = malloc(N*sizeof(*u));
        double (*v)[N] = malloc(N*sizeof(*v));
        double (*p)[N] = malloc(N*sizeof(*p));
        double (*X)[N] = malloc(N*sizeof(*X));
        double (*Y)[N] = malloc(N*sizeof(*Y));
        double (*I)[N] = malloc(N*sizeof(*I));
        double (*D)[N] = malloc(N*sizeof(*D));
        double (*GD)[N] = malloc(N*sizeof(*GD));
        double (*NGD)[N] = malloc(N*sizeof(*NGD));
        double (*ND)[N] = malloc(N*sizeof(*ND));

    // (4N x N) matrices

		double (*V)[N] = malloc(4*N*sizeof(*V));
        double (*W)[N] = malloc(4*N*sizeof(*W));

		// V_temp and W_temp store V and W temporarily

    	double (*V_temp)[N] = malloc(4*N*sizeof(*V_temp));
        double (*W_temp)[N] = malloc(4*N*sizeof(*W_temp));

	// (4N x 4N) matrices

		double (*MLX)[4*N] = malloc(4*N*sizeof(*MLX));
        double (*MLY)[4*N] = malloc(4*N*sizeof(*MLY));
        double (*MRX)[4*N] = malloc(4*N*sizeof(*MRX));
        double (*MRY)[4*N] = malloc(4*N*sizeof(*MRY));

    // (N x 4N) matrices

    // V_T_temp, W_T_temp temporarily store V transpose, U transpose

        double (*V_T_temp)[4*N] = malloc(N*sizeof(*V_T_temp));
        double (*W_T_temp)[4*N] = malloc(N*sizeof(*W_T_temp));

    // Populate X, Y, rho, p matrices

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){

			X[i][j] = -1 + dx * i;
			Y[i][j] = -1 + dx * j;
			rho[i][j] = 2/gam * exp(-100 * (X[i][j]*X[i][j] + Y[i][j]*Y[i][j]));
			p[i][j] = 2 * exp(-100 * (pow(X[i][j],2) + pow(Y[i][j],2)));
		}
	}

	// Store matrices in V, block form

	for (int i = 0; i < 4 * N; i++){
		for (int j = 0; j < N; j++){

			if (i >= 0 && i < N){
				V[i][j] = rho[i][j];
			}else if (i >= N && i < 2*N){
				V[i][j] = 0;
			}else if (i >= 2*N && i < 3*N){
				V[i][j] = 0;
			}else
				V[i][j] = p[i-3*N][j];
			}
	}

	// Populate identity matrix

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			if (i == j){
				I[i][j] = 1;
			}else{
				I[i][j] = 0;
			}
		}
	}

	// Populate D matrix

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){

			if (i == j-1){
				D[i][j] = -dt/(4*dx);
			}else if (i == j+1){
				D[i][j] = +dt/(4*dx);
			}else if (i == N-1 && j == 0){
				D[i][j] = -dt/(4*dx);
			}else if (i == 0 && j == N-1){
				D[i][j] = +dt/(4*dx);
			}else
				D[i][j] = 0;
		}
	}

	// Populate -D matrix

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			ND[i][j]= -1* D[i][j];
		}
	}

	// Populate product matrix, GD = (gamma*D)

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			GD[i][j]= gam * D[i][j];
		}
	}

	// Populate product matrix, NGD = -GD

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			NGD[i][j]= - GD[i][j];
		}
	}

	// Populate MRX matrix

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
            MRX[i][j] = I[i][j];
            MRX[i][N+j] = 0;
            MRX[i][2*N+j] = 0;
            MRX[i][3*N+j] = 0;
            MRX[N+i][j] = ND[i][j];
            MRX[N+i][N+j] = I[i][j];
            MRX[N+i][2*N+j] = 0;
            MRX[N+i][3*N+j] = NGD[i][j];
            MRX[2*N+i][j]=0;
            MRX[2*N+i][N+j]=0;
            MRX[2*N+i][2*N+j]=I[i][j];
            MRX[2*N+i][3*N+j]=0;
            MRX[3*N+i][j]=0;
            MRX[3*N+i][N+j]=ND[i][j];
            MRX[3*N+i][2*N+j]=0;
            MRX[3*N+i][3*N+j]=I[i][j];
		}
	}

	// Populate MLX matrix

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
            MLX[i][j]=I[i][j];
            MLX[i][N+j]=0;
            MLX[i][2*N+j]=0;
            MLX[i][3*N+j]=0;
            MLX[N+i][j]=D[i][j];
            MLX[N+i][N+j]=I[i][j];
            MLX[N+i][2*N+j]=0;
            MLX[N+i][3*N+j]=GD[i][j];
            MLX[2*N+i][j]=0;
            MLX[2*N+i][N+j]=0;
            MLX[2*N+i][2*N+j]=I[i][j];
            MLX[2*N+i][3*N+j]=0;
            MLX[3*N+i][j]=0;
            MLX[3*N+i][N+j]=D[i][j];
            MLX[3*N+i][2*N+j]=0;
            MLX[3*N+i][3*N+j]=I[i][j];
		}
	}

	// Populate MRY matrix

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
            MRY[i][j]=I[i][j];
            MRY[i][N+j]=0;
            MRY[i][2*N+j]=0;
            MRY[i][3*N+j]=0;
            MRY[N+i][j]=0;
            MRY[N+i][N+j]=I[i][j];
            MRY[N+i][2*N+j]=0;
            MRY[N+i][3*N+j]=0;
            MRY[2*N+i][j]=ND[i][j];
            MRY[2*N+i][N+j]=0;
            MRY[2*N+i][2*N+j]=I[i][j];
            MRY[2*N+i][3*N+j]=NGD[i][j];
            MRY[3*N+i][j]=0;
            MRY[3*N+i][N+j]=0;
            MRY[3*N+i][2*N+j]=ND[i][j];
            MRY[3*N+i][3*N+j]=I[i][j];
		}
	}

	// Populate MLY matrix

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
            MLY[i][j] = I[i][j];
            MLY[i][N+j]=0;
            MLY[i][2*N+j]=0;
            MLY[i][3*N+j]=0;
            MLY[N+i][j]=0;
            MLY[N+i][N+j]=I[i][j];
            MLY[N+i][2*N+j]=0;
            MLY[N+i][3*N+j]=0;
            MLY[2*N+i][j]=D[i][j];
            MLY[2*N+i][N+j]=0;
            MLY[2*N+i][2*N+j]=I[i][j];
            MLY[2*N+i][3*N+j]=GD[i][j];
            MLY[3*N+i][j]=0;
            MLY[3*N+i][N+j]=0;
            MLY[3*N+i][2*N+j]=D[i][j];
            MLY[3*N+i][3*N+j]=I[i][j];
        }
    }

        int key = 0;
		int dim = 4*N;
		int lda = 4*N;
        int nrhs = 4*N;
        int ldb = 4*N;
        int ipiv1[4*N];
        int ipiv2[4*N];

        char trb = 'N';

		// LU factorization for MLY and MLX

		double *A = &(MLY[0][0]);
		double *C = &(MLX[0][0]);
		dgetrf_(&dim,&dim,A,&lda,ipiv1,&key);
		dgetrf_(&dim,&dim,C,&lda,ipiv2,&key);

		// Open target files

		FILE *rhodata = fopen("R.out","w");
		FILE *udata = fopen("U.out","w");
		FILE *vdata = fopen("V.out","w");
		FILE *pdata = fopen("P.out","w");

	double t = 0.0;


	// START MAIN LOOP
	for(int n = 0; n < M+1 ; n++){

		// Forward in x with V matrix, V_T_temp=MRX*V

        for (int i = 0;i < 4*N;i++){
        	for (int j = 0; j < N;j++){
        		for (int k = 0; k < 4*N; k++){
        			V_temp[i][j] += MRX[k][i] * V[k][j];
        		}
        	}
        }
        for (int i = 0; i < 4*N; i++){
        	for (int j = 0; j < N;j++){
        		V[i][j] = V_temp[i][j];
        	}
        }
        for (int i = 0; i < 4*N; i++){
			for(int j = 0; j < N; j++){
				V_temp[i][j] = 0;
			}
		}



		// Backward in y, relating W to V


        for (int i = 0; i < 4*N; i++){
        	for (int j = 0; j < N; j++){
        		 if (i < N){
        		 	W[i][j] = V[j][i];
        		 }else if (i >= N && i < 2*N){
        		 	W[i][j] = V[j+N][i-N];
        		 }else if (i >= 2*N && i < 3*N){
        		 	W[i][j] = V[j+2*N][i-2*N];
        		 }else{
        		 	W[i][j] = V[j+3*N][i-3*N];
        		 }
        	}
        }

        // Transpose W_T_temp to get back updated W

        for (int i = 0; i< 4*N; i++){
        	for (int j=0; j<N; j++){
            	W_T_temp[j][i] = W[i][j];
        	}
    	}

		// Advance with W = MLY\W
		// Segfault occurs in following loop

        // return 0;
		double *B = &(W_T_temp[0][0]);

		// lapack linear solver
		dgetrs_(&trb,&dim,&nrhs,A,&lda,ipiv1,B,&ldb,&key);

		for (int i =0; i<N; ++i){
        	for (int j=0; j<4*N; ++j){
            	W[j][i] = W_T_temp[i][j];
        	}
    	}


		// Forward in y with W matrix, W_T_temp = MRY*W

        for (int i = 0; i < 4*N;i++){
        	for (int j = 0; j < N;j++){
        		for (int k = 0; k < 4*N; k++){
        			W_temp[i][j] += MRY[k][i] * W[k][j];
        		}
        	}
        }
        for (int i = 0; i < 4*N;i++){
        	for (int j = 0; j < N;j++){
        		W[i][j] = W_temp[i][j];
        	}
        }
        for (int i = 0; i < 4*N; i++){
			for(int j = 0; j < N; j++){
				W_temp[i][j] = 0;
			}
		}

		// Backward in x, relating V to W

		for (int i = 0; i < 4*N; i++){
        	for (int j = 0; j < N; j++){
        		 if (i < N){
        		 	V[i][j] = W[j][i];
        		 }else if (i >= N && i < 2*N){
        		 	V[i][j] = W[j+N][i-N];
        		 }else if (i >= 2*N && i < 3*N){
        		 	V[i][j] = W[j+2*N][i-2*N];
        		 }else{
        		 	V[i][j] = W[j+3*N][i-3*N];
        		 }
        	}
        }

        for (int i =0; i<4*N; ++i){
        	for (int j=0; j<N; ++j){
            	V_T_temp[j][i] = V[i][j];
        	}
    	}

        // Advance with V = MLX\V

		double *E = &(V_T_temp[0][0]);
		dgetrs_(&trb,&dim,&nrhs,C,&lda,ipiv2,E,&ldb,&key);
		for (int i =0; i<N; ++i){
        	for (int j=0; j<4*N; ++j){
            	V[j][i] = V_T_temp[i][j];
        	}
    	}

    	for (int i= 0; i < 4*N; i++){
        	for (int j = 0; j< N; j++){
            	V_T_temp[j][i]=0;
            	W_T_temp[j][i]=0;
        	}
    	}

		if (t >= 0.2){

            // Reset t

            t = 0;

			for (int i = 0; i < N; i++){
				for (int j = 0; j < N; j++){

					double rho1 = V[i][j];
					fwrite(&rho1, sizeof(double),1,rhodata);
				}
			}

			for (int i = N; i < 2*N; i++){
				for (int j = 0; j < N; j++){

					double u1 = V[i][j];
					fwrite(&u1, sizeof(double),1,udata);

				}
			}

			for (int i = 2*N; i < 3*N; i++){
				for (int j = 0; j < N; j++){

					double v1 = V[i][j];
					fwrite(&v1, sizeof(double),1,vdata);

				}
			}

			for (int i = 3*N; i < 4*N; i++){
				for (int j = 0; j < N; j++){

					double p1 = V[i][j];
					fwrite(&p1, sizeof(double),1,pdata);

				}
			}

		}
    // advance in time at end of each loop run

    t += dt;
	}

	// END MAIN LOOP

    // Close files used in main loop

    fclose(rhodata);
    fclose(udata);
    fclose(vdata);
    fclose(pdata);

    // End clock, calculate runtime

	clock_t end = clock();
	double runtime = (double)(end - start)/CLOCKS_PER_SEC;

	// Output runtime to console

	printf("Runtime: %1.4f s\n", runtime);
	return 0;
}
