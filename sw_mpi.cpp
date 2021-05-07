#include <iostream>
#include <math.h>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <mpi.h>

#define PI 3.14159265f

#define g 1.0f
#define f 0.1f
#define L 100.0f
#define hx 1.0f
#define hy 1.0f
#define T 100.0f
#define dt 0.01f
#define H0 20000.0f
#define H1 4400.0f
#define H2 2660.0f
#define D 4400000.0f

using namespace std;


void init(float *U, float *V, float *H); // ham khoi tao

void writeResult(float *U, float *V, float *H, float t); // ham ghi ket qua ra file

int nx = L/hx;
int ny = L/hy;

//=========================
void truyenTren(float *Uu, float *Vv,float *Hh, float *Ud, float *Vd,float *Hd, int rank, int size){
    MPI_Status stat;
    int src = rank + 1;
    int dst = rank - 1;
    if(rank==0){
        dst = size-1;

    }
    else if(rank == size-1){
        src = 0;
    }
    MPI_Send(Uu, nx, MPI_FLOAT, dst, rank, MPI_COMM_WORLD);
    MPI_Recv(Ud, nx, MPI_FLOAT, src, src, MPI_COMM_WORLD, &stat);
    MPI_Send(Vv, nx, MPI_FLOAT, dst, rank, MPI_COMM_WORLD);
    MPI_Recv(Vd, nx, MPI_FLOAT, src, src, MPI_COMM_WORLD, &stat);
    MPI_Send(Hh, nx, MPI_FLOAT, dst, rank, MPI_COMM_WORLD);
    MPI_Recv(Hd, nx, MPI_FLOAT, src, src, MPI_COMM_WORLD, &stat);
}
//=========================
void truyenDuoi(float *Uu, float *Vv,float *Hh, float *Ut, float *Vt,float *Ht, int rank, int size){
    MPI_Status stat;
    int src = rank - 1;
    int dst = rank + 1;

    if(rank==0){
        src = size - 1;
    }
    else if(rank == size-1){
        dst = 0;
    }
    MPI_Send(Uu + (ny/size-1)*nx, nx, MPI_FLOAT, dst, rank, MPI_COMM_WORLD);
    MPI_Recv(Ut, nx, MPI_FLOAT, src, src, MPI_COMM_WORLD, &stat);
    MPI_Send(Vv + (ny/size-1)*nx, nx, MPI_FLOAT, dst, rank, MPI_COMM_WORLD);
    MPI_Recv(Vt, nx, MPI_FLOAT, src, src, MPI_COMM_WORLD, &stat);
    MPI_Send(Hh + (ny/size-1)*nx, nx, MPI_FLOAT, dst, rank, MPI_COMM_WORLD);
    MPI_Recv(Ht, nx, MPI_FLOAT, src, src, MPI_COMM_WORLD, &stat);
}
//=========================

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

float dy(const float* A, int i, int j, float *At, float*Ad, int nyc) {
    float u = (j == nyc-1) ? *(Ad + i): value(A, i, j+1);
    float d = (j == 0) ? *(At + i) : value(A, i, j-1);

    return (u - d) / (2 * hy);
}

void calDerivate(float* U, float* V, float* H, float *dU, float *dV, float *dH, float *Ut, float *Vt, float *Ht, float *Ud, float *Vd, float *Hd, int nyc){
    // float Ux, Uy, Ur, Ud, Uc, Vx, Vy, Vr, Vd, Vc, Hc, Hd, Hr, Hx, Hy;

    for (int j = 0; j < nyc; j++){
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
            float Uy = dy(U, i, j, Ut, Ud, nyc);
            float Vx = dx(V, i, j);
            float Vy = dy(V, i, j, Vt, Vd, nyc);
            float Hx = dx(H, i, j);
            float Hy = dy(H, i, j, Ht, Hd, nyc);

            // *(dU + index(i, j)) = f*Vc - Uc*Ux - Vc*Uy -g*Hx;
            // *(dV + index(i, j)) = -f*Uc - Uc*Vx - Vc*Vy - g*Hy;
            // *(dH + index(i, j)) = -Uc*Hx - Hc*Ux - Vc*Hy - Hc*Vy;

            *(dU + index(i, j)) = - g * Hx;
            *(dV + index(i, j)) = - g * Hy;
            *(dH + index(i, j)) = - (Ux * value(H, i, j) + Vy * value(H, i, j));
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
    *(H+index(nx/2, ny/2)) = 1.3;

    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
            *(U + index(i, j)) = 0;
            *(V + index(i, j)) = 0;
        }
    }
}

void writeResult(float *U, float *V, float *H, float t){
    std::fstream output;
	// output.open("result/outputU.txt", ios::app);
    // output <<"time: " << t << endl;
    // output << setprecision(16);

    // for (int j = 0; j < ny; j++){
    //         for (int i = 0; i < nx; i++){
    //         output<<value(U, i, j) << " ";
    //     }
    //     output<<endl;
    // }
    // output.close();

    // output.open("result/outputV.txt", ios::app);
    // output <<"time: " << t << endl;
    // output << setprecision(16);

    // for (int j = 0; j < ny; j++){
    //         for (int i = 0; i < nx; i++){
    //         output<<value(V, i, j) << " ";
    //     }
    //     output<<endl;
    // }
    // output.close();

    string h_file = "result_mpi/outputH_" + to_string((int)round(t/dt)) + ".txt";
    output.open(h_file, ios::out | ios::ate);
    output <<"time: " << t << endl;
    output << setprecision(32);

    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
                float val = value(H, i, j);
                // val = max(1.f, min(val, 0.f));
                output<< val << " ";
        }
        output<<endl;
    }
    output.close();
} // ham ghi ket qua ra file


int main(int argc, char **argv){
    float *U, *V, *H;
    float *dU, *dV, *dH;
    float *Uu, *Vv, *Hh;
    float *Ut, *Vt, *Ht;
    float *Ud, *Vd, *Hd;

    float t = 0;
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status stat;

    int nyc = ny/size;
    if (rank == 0){
        U= (float *) malloc ((nx*ny)*sizeof(float));
        V= (float *) malloc ((nx*ny)*sizeof(float));
        H= (float *) malloc ((nx*ny)*sizeof(float));

        init(U, V, H);
    }
    int step = 0;
    Uu= (float *) malloc ((nx*nyc)*sizeof(float));
    Vv= (float *) malloc ((nx*nyc)*sizeof(float));
    Hh= (float *) malloc ((nx*nyc)*sizeof(float));
    dU= (float *) malloc ((nx*nyc)*sizeof(float));
    dV= (float *) malloc ((nx*nyc)*sizeof(float));
    dH= (float *) malloc ((nx*nyc)*sizeof(float));
    
    Ut= (float *) malloc ((nx)*sizeof(float));
    Vt= (float *) malloc ((nx)*sizeof(float));
    Ht= (float *) malloc ((nx)*sizeof(float));

    Ud= (float *) malloc ((nx)*sizeof(float));
    Vd= (float *) malloc ((nx)*sizeof(float));
    Hd= (float *) malloc ((nx)*sizeof(float));


    MPI_Scatter(U, nx*nyc, MPI_FLOAT, Uu, nx*nyc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(V, nx*nyc, MPI_FLOAT, Vv, nx*nyc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(H, nx*nyc, MPI_FLOAT, Hh, nx*nyc, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (rank == 0){
        writeResult(U, V, H, t);
    }

    while (t <= T)
    {
        truyenDuoi(Uu, Vv, Hh, Ut, Vt, Ht, rank, size);
        truyenTren(Uu, Vv, Hh, Ud, Vd, Hd, rank, size);
        MPI_Barrier(MPI_COMM_WORLD);
        calDerivate(Uu, Vv, Hh, dU, dV, dH, Ut, Vt, Ht, Ud, Vd, Hd, nyc);
        for (int j = 0; j < nyc; j++){
            for (int i = 0; i < nx; i++){
                *(Uu + index(i,j)) += dt*(value(dU, i, j));
                *(Vv + index(i,j)) += dt*(value(dV, i, j));
                *(Hh + index(i,j)) += dt*(value(dH, i, j));
            }
        }
        step ++;
        t += dt;

        if (step % 10 == 0){
            MPI_Gather(Uu, nx*nyc, MPI_FLOAT, U, nx*nyc, MPI_FLOAT, 0, MPI_COMM_WORLD);
            MPI_Gather(Vv, nx*nyc, MPI_FLOAT, V, nx*nyc, MPI_FLOAT, 0, MPI_COMM_WORLD);
            MPI_Gather(Hh, nx*nyc, MPI_FLOAT, H, nx*nyc, MPI_FLOAT, 0, MPI_COMM_WORLD);
            
            if (rank == 0){
                writeResult(U, V, H, t);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Gather(Uu, nx*nyc, MPI_FLOAT, U, nx*nyc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(Vv, nx*nyc, MPI_FLOAT, V, nx*nyc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(Hh, nx*nyc, MPI_FLOAT, H, nx*nyc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(rank==0){
        writeResult(U, V, H, t);
    }

    MPI_Finalize();

    return 0;
}
