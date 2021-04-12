#include<iostream>
#include<math.h>
#include <fstream>
#include <iomanip>
#include <cstdio>

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

using namespace std;

void calDerivate(float* U, float* V, float* H, float *dU, float *dV, float *dH); // ham tinh dao ham theo t

void init(float *U, float *V, float *H); // ham khoi tao

void writeResult(float *U, float *V, float *H, float t); // ham ghi ket qua ra file

int nx = L/hx;
int ny = L/hy;

int index(int i, int j){
    return j * nx + i;
}

float value(const float *A, int i, int j){
    return *(A + index(i, j));
}

int main(){
    float *U, *V, *H;
    float *dU, *dV, *dH;
    float t = 0;
    
    U= (float *) malloc ((nx*ny)*sizeof(float));
    V= (float *) malloc ((nx*ny)*sizeof(float));
    H= (float *) malloc ((nx*ny)*sizeof(float));

    dU= (float *) malloc ((nx*ny)*sizeof(float));
    dV= (float *) malloc ((nx*ny)*sizeof(float));
    dH= (float *) malloc ((nx*ny)*sizeof(float));

    remove("outputU.txt");
    remove("outputV.txt");
    remove("outputH.txt");

    init(U, V, H);
    writeResult(U, V, H, t);
    int step = 0;

    while (t <= T)
    {
        calDerivate(U, V, H, dU, dV, dH);
        for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
                *(U + index(i,j)) += dt*(value(dU, i, j));
                *(V + index(i,j)) += dt*(value(dV, i, j));
                *(H + index(i,j)) += dt*(value(dH, i, j));
            }
        }

        t += dt;
        step ++;
        if(step % 10 == 0)
            writeResult(U, V, H, t);
    }
    writeResult(U, V, H, t);

    return 0;
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

void calDerivate(float* U, float* V, float* H, float *dU, float *dV, float *dH){
    // float Ux, Uy, Ur, Ud, Uc, Vx, Vy, Vr, Vd, Vc, Hc, Hd, Hr, Hx, Hy;

    for (int j = 0; j < ny; j++){
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

void writeResult(float *U, float *V, float *H, float t){
    fstream output;
	output.open("outputU.txt", ios::app);
    output <<"time: " << t << endl;
    output << setprecision(16);

    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
            output<<value(U, i, j) << " ";
        }
        output<<endl;
    }
    output.close();

    output.open("outputV.txt", ios::app);
    output <<"time: " << t << endl;
    output << setprecision(16);

    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
            output<<value(V, i, j) << " ";
        }
        output<<endl;
    }
    output.close();

    output.open("outputH.txt", ios::app);
    output <<"time: " << t << endl;
    output << setprecision(16);

    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
            output<<value(H, i, j) << " ";
        }
        output<<endl;
    }
    output.close();
} // ham ghi ket qua ra file