#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <cstdio>

using namespace std;

#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */

#define PI 3.14159265f
#define g 1.0f
#define f 0.1f
#define L 100.0f
#define hx 1.0f
#define hy 1.0f
#define T 10.0f
#define dt 0.01f
#define H0 20000.0f
#define H1 4400.0f
#define H2 2660.0f
#define D 4400000.0f

int nx = L/hx;
int ny = L/hy;

int index(int i, int j){
    return j * nx + i;
}

float value(const float *A, int i, int j){
    return *(A + index(i, j));
}

float dx(const float* A, int i, int j) {
    float r = (i == nx-1) ? value(A, 0, j) : value(A, i+1, j);
    float l = (i == 0) ? value(A, nx-1, j) : value(A, i-1, j);

    return (r - l) / (2 * hx);
}

float dy(const float* A, int i, int j) {
    float u = (j == ny-1) ? value(A, i, 0) : value(A, i, j+1);
    float d = (j == 0) ? value(A, i, ny-1) : value(A, i, j-1);

    return (u - d) / (2 * hy);
}

void calDerivate(float* U, float* V, float* H, float *dU, float *dV, float *dH, int b, int e){
    // float Ux, Uy, Ur, Ud, Uc, Vx, Vy, Vr, Vd, Vc, Hc, Hd, Hr, Hx, Hy;

    for (int j = b; j < e; j++){
            for (int i = 0; i < nx; i++){
            // Uc = value(U, i ,j);
            // Ur = (i == nx-1) ? value(U, 0, j) : value(U, i+1, j);
            // Ud = (j == ny-1) ? value(U, i, 0) : value(U, i, j+1);
            // // Ud = (j == ny-1) ? 0 : value(U, i, j+1);

            // Vc = value(V, i,j);
            // Vr = (i == nx-1) ? value(V, 0, j) : value(V, i+1, j);
            // Vd = (j == ny-1) ? value(V, i, 0) : value(V, i, j+1);
            // // Vd = (j == ny-1) ? 0 : value(V, i, j+1);

            // Hc = value(H, i, j);
            // Hr = (i == nx-1) ? value(H, 0, j) : value(H, i+1, j);
            // Hd = (j == ny-1) ? value(H, i, 0) : value(H, i, j+1);
            // // Hd = (j == ny-1) ? 0 : value(H, i, j+1);

            // Ux = (Ur - Uc)/hx;
            // Uy = (Ud - Uc)/hy;

            // Vx = (Vr - Vc)/hx;
            // Vy = (Vd - Vc)/hy;

            // Hx = (Hr - Hc)/hx;
            // Hy = (Hd - Hc)/hy;

            float Ux = dx(U, i, j);
            float Uy = dy(U, i, j);
            float Vx = dx(V, i, j);
            float Vy = dy(V, i, j);
            float Hx = dx(H, i, j);
            float Hy = dy(H, i, j);

            // *(dU + index(i, j)) = f*Vc - Uc*Ux - Vc*Uy -g*Hx;
            // *(dV + index(i, j)) = -f*Uc - Uc*Vx - Vc*Vy - g*Hy;
            // *(dH + index(i, j)) = -Uc*Hx - Hc*Ux - Vc*Hy - Hc*Vy;

            *(dU + index(i, j)) = - g * Hx;
            *(dV + index(i, j)) = - g * Hy;
            *(dH + index(i, j)) = - (Ux * value(H, i, j) + Hx * value(U, i, j) + Vy * value(H, i, j) + Hy * value(V, i, j)) ;
        }
    }
}

void update(float *U, float *V, float *H, float *dU, float *dV, float *dH, int b, int e) {
	calDerivate(U, V, H, dU, dV, dH, b, e);

	for (int j = b; j < e; j++){
		for (int i = 0; i < nx; i++){
			*(U + index(i,j)) += dt*(value(dU, i, j));
			*(V + index(i,j)) += dt*(value(dV, i, j));
			*(H + index(i,j)) += dt*(value(dH, i, j));
		}
	}
}

void init(float *U, float *V, float *H){
    float Hc, Hd, Hr, Hx, Hy;
    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
            *(H + index(i, j)) = 1.0;
        }
    }
    *(H) = 1.1;

    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
            *(U + index(i, j)) = 0;
            *(V + index(i, j)) = 0;
        }
    }
}

// void writeResult(float *U, float *V, float *H, float t){
//     fstream output;
// 	output.open("outputU.txt", ios::app);
//     output <<"time: " << t << endl;
//     output << setprecision(16);

//     for (int j = 0; j < ny; j++){
//             for (int i = 0; i < nx; i++){
//             output<<value(U, i, j) << " ";
//         }
//         output<<endl;
//     }
//     output.close();

//     output.open("outputV.txt", ios::app);
//     output <<"time: " << t << endl;
//     output << setprecision(16);

//     for (int j = 0; j < ny; j++){
//             for (int i = 0; i < nx; i++){
//             output<<value(V, i, j) << " ";
//         }
//         output<<endl;
//     }
//     output.close();

//     output.open("outputH.txt", ios::app);
//     output <<"time: " << t << endl;
//     output << setprecision(16);

//     for (int j = 0; j < ny; j++){
//             for (int i = 0; i < nx; i++){
//             output<<value(H, i, j) << " ";
//         }
//         output<<endl;
//     }
//     output.close();
// } // ham ghi ket qua ra file

int main(int argc, char** argv) {
	int id;
	int numtasks;
	float *U, *V, *H;
    float *dU, *dV, *dH;
	int i, it;
	int numworkers;
	int dst, src;
	int left, right;
	int avgnum, extra, offset;
	int numele;
	int Ntime;
	MPI_Status status;

	U = (float *) malloc ((nx*ny)*sizeof(float));
    V = (float *) malloc ((nx*ny)*sizeof(float));
    H = (float *) malloc ((nx*ny)*sizeof(float));

    dU = (float *) malloc ((nx*ny)*sizeof(float));
    dV = (float *) malloc ((nx*ny)*sizeof(float));
    dH = (float *) malloc ((nx*ny)*sizeof(float));

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	numworkers = numtasks - 1;

	if (id == MASTER) {
		init(U, V, H);

		avgnum = nx / numworkers;
		extra = nx % numworkers;
		offset = 0;

		for (i = 1; i <= numworkers; i++) {
			numele = (i <= extra) ? avgnum + 1 : avgnum;
			numele *= ny;
			left = (i == 1) ? NONE : i - 1;
			right = (i == numworkers) ? NONE : i + 1;

			dst = i;
			MPI_Send(&offset, 1, MPI_INT, dst, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&numele, 1, MPI_INT, dst, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&left, 1, MPI_INT, dst, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&right, 1, MPI_INT, dst, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&U[offset], numele, MPI_FLOAT, dst, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&V[offset], numele, MPI_FLOAT, dst, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&H[offset], numele, MPI_FLOAT, dst, BEGIN, MPI_COMM_WORLD);

			offset = offset + numele;
		}

		for (i = 1; i <= numworkers; i++) {
			src = i;
			MPI_Recv(&offset, 1, MPI_INT, src, DONE, MPI_COMM_WORLD, &status);
			MPI_Recv(&numele, 1, MPI_INT, src, DONE, MPI_COMM_WORLD, &status);
			MPI_Recv(&U[offset], numele, MPI_FLOAT, src, DONE, MPI_COMM_WORLD, &status);
			MPI_Recv(&V[offset], numele, MPI_FLOAT, src, DONE, MPI_COMM_WORLD, &status);
			MPI_Recv(&H[offset], numele, MPI_FLOAT, src, DONE, MPI_COMM_WORLD, &status);
		}

		// writeResult(U, V, H, 0); // currently write only last result

		MPI_Finalize();
	}

	if (id != MASTER) {
		src = MASTER;

		MPI_Recv(&offset, 1, MPI_INT, src, BEGIN, MPI_COMM_WORLD, &status);
		MPI_Recv(&numele, 1, MPI_INT, src, BEGIN, MPI_COMM_WORLD, &status);
		MPI_Recv(&left, 1, MPI_INT, src, BEGIN, MPI_COMM_WORLD, &status);
		MPI_Recv(&right, 1, MPI_INT, src, BEGIN, MPI_COMM_WORLD, &status);
		MPI_Recv(&U[offset], numele, MPI_FLOAT, src, BEGIN, MPI_COMM_WORLD, &status);
		MPI_Recv(&V[offset], numele, MPI_FLOAT, src, BEGIN, MPI_COMM_WORLD, &status);
		MPI_Recv(&H[offset], numele, MPI_FLOAT, src, BEGIN, MPI_COMM_WORLD, &status);

		//printf("id = %d, offset = %d, numele = %d, left = %d, right = %d\n", id, offset, numele, left, right);
		//display(C, offset, offset + (int)numele/n);

		Ntime = T / dt;
		for (it = 0; it < Ntime; it++) {
			if (left != NONE) {
				MPI_Send(&U[offset], ny, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);
				MPI_Send(&V[offset], ny, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);
				MPI_Send(&H[offset], ny, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);
				MPI_Send(&dU[offset], ny, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);
				MPI_Send(&dV[offset], ny, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);
				MPI_Send(&dH[offset], ny, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);

				src = left;
				MPI_Recv(&U[offset - ny], ny, MPI_FLOAT, src, LTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&V[offset - ny], ny, MPI_FLOAT, src, LTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&H[offset - ny], ny, MPI_FLOAT, src, LTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&dU[offset - ny], ny, MPI_FLOAT, src, LTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&dV[offset - ny], ny, MPI_FLOAT, src, LTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&dH[offset - ny], ny, MPI_FLOAT, src, LTAG, MPI_COMM_WORLD, &status);
			}

			if (right != NONE) {
				MPI_Send(&U[offset + numele - ny], ny, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);
				MPI_Send(&V[offset + numele - ny], ny, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);
				MPI_Send(&H[offset + numele - ny], ny, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);
				MPI_Send(&dU[offset + numele - ny], ny, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);
				MPI_Send(&dV[offset + numele - ny], ny, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);
				MPI_Send(&dH[offset + numele - ny], ny, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);

				src = right;
				MPI_Recv(&U[offset + numele], ny, MPI_FLOAT, src, RTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&V[offset + numele], ny, MPI_FLOAT, src, RTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&H[offset + numele], ny, MPI_FLOAT, src, RTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&dU[offset + numele], ny, MPI_FLOAT, src, RTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&dV[offset + numele], ny, MPI_FLOAT, src, RTAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&dH[offset + numele], ny, MPI_FLOAT, src, RTAG, MPI_COMM_WORLD, &status);
			}

			update(U, V, H, dU, dV, dH, offset/ny, (int)((offset + numele)/ny));
		}

		MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
		MPI_Send(&numele, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
		MPI_Send(&U[offset], numele, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
		MPI_Send(&V[offset], numele, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
		MPI_Send(&H[offset], numele, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);

		MPI_Finalize();
	}

	return 0;
}