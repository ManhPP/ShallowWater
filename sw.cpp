#include<iostream>
#include<math.h>
#include <fstream>
#include <iomanip>
#include <cstdio>

#define PI 3.14159265f

#define g 9.8f
#define f 0.1f
#define L 10.0f
#define hx 5.0f
#define hy 5.0f
#define T 1.0f
#define dt 0.5f

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

float value(float *A, int i, int j){
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
        writeResult(U, V, H, t);
    }
    writeResult(U, V, H, t);

    return 0;
}

void calDerivate(float* U, float* V, float* H, float *dU, float *dV, float *dH){
    float Ux, Uy, Ur, Ud, Uc, Vx, Vy, Vr, Vd, Vc, Hc, Hd, Hr, Hx, Hy;

    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
            Uc = value(U, i ,j);
            Ur = (i == nx-1) ? value(U, 0, j) : value(U, i+1, j);
            Ud = (j == ny-1) ? 0.0 : value(U, i, j+1);

            Vc = value(V, i,j);
            Vr = (i == nx-1) ? value(V, 0, j) : value(V, i+1, j);
            Vd = (j == ny-1) ? 0.0 : value(V, i, j+1);

            Hc = value(H, i, j);
            Hr = (i == nx-1) ? value(H, 0, j) : value(H, i+1, j);
            Hd = (j == ny-1) ? 0.0 : value(H, i, j+1);

            Ux = (Ur - Uc)/hx;
            Uy = (Ud - Uc)/hy;

            Vx = (Vr - Vc)/hx;
            Vy = (Vd - Vc)/hy;

            Hx = (Hr - Hc)/hx;
            Hy = (Hd - Hc)/hy;

            *(dU + index(i, j)) = f*Vc - Uc*Ux - Vc*Uy -g*Hx;
            *(dV + index(i, j)) = -f*Uc - Uc*Vx - Vc*Vy - g*Hy;
            *(dH + index(i, j)) = -Uc*Hx - Hc*Ux - Vc*Hy - Hc*Vy;
        }
    }
} // ham tinh dao ham theo t

void init(float *U, float *V, float *H){
    float Hc, Hd, Hr, Hx, Hy;
    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
            *(H + index(i, j)) = H0 
                                + H1*tan(9*(j*hy - 6*hx)/(2*D)) 
                                + H2*sin(2*PI*i*hx)/pow(cos((9*j*hy - 6*hx)/D), 2);
        }
    }
    for (int j = 0; j < ny; j++){
            for (int i = 0; i < nx; i++){
            Hc = value(H, i, j);
            Hr = (i == nx-1) ? value(H, 0, j) : value(H, i+1, j);
            Hd = (j == ny-1) ? 0 : value(H, i, j+1);

            Hx = (Hr - Hc)/hx;
            Hy = (Hd - Hc)/hy;

            *(U + index(i, j)) = -Hy/f;
            *(V + index(i, j)) = -Hx/f;
        }
    }
} // ham khoi tao

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
