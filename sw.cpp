#include<iostream>
#include<math.h>

#define PI 3.14159265

#define g 9.8
#define f 0.1
#define L 10
#define hx 1
#define hy 1
#define T 2
#define dt 0.01

#define H0 20000
#define H1 4400
#define H2 2660
#define D 4400000

using namespace std;

void calDerivate(float* U, float* V, float* H, float *dU, float *dV, float *dH); // ham tinh dao ham theo t

void init(float *U, float *V, float *H); // ham khoi tao

void writeResult(float *U, float *V, float *H); // ham ghi ket qua ra file

int nx = L/hx;
int ny = L/hy;

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


    init(U, V, H);
    writeResult(U, V, H);

    while (t <= T)
    {
        calDerivate(U, V, H, dU, dV, dH);

        for (int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                *(U + i*ny + j) += dt*(*(dU + i*ny + j));
                *(V + i*ny + j) += dt*(*(dV + i*ny + j));
                *(H + i*ny + j) += dt*(*(dH + i*ny + j));
            }
        }

        t += dt;
    }

    writeResult(U, V, H);
    return 0;
}

void calDerivate(float* U, float* V, float* H, float *dU, float *dV, float *dH){
    float Ux, Uy, Ur, Ud, Uc, Vx, Vy, Vr, Vd, Vc, Hc, Hd, Hr, Hx, Hy;
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){

            //U
            Uc = *(U + i*ny + j);
            Ur = (i==nx-1) ? *(U + j) : *(U + (i+1)*ny + j);
            Ud = (j==ny-1) ? 0 : *(U + i*ny + j+1);

            Ux = (Ur - Uc)/hx;
            Uy = (Ud - Uc)/hy;

            //V
            Vc = *(V + i*ny + j);
            Vr = (i==nx-1) ? *(V + j) : *(V + (i+1)*ny + j);
            Vd = (j==ny-1) ? 0 : *(V + i*ny + j+1);

            Vx = (Vr - Vc)/hx;
            Vy = (Vd - Vc)/hy;

            //H
            Hc = *(H + i*ny + j);
            Hr = (i==nx-1) ? *(H + j) : *(H + (i+1)*ny + j);
            Hd = (j==ny-1) ? 0 : *(H + i*ny + j+1);

            Hx = (Hr - Hc)/hx;
            Hy = (Hd - Hc)/hy;

            *(dU + i*ny + j) = f*Vc - Uc*Ux - Vc*Uy - g*Hx;
            *(dV + i*ny + j) = -f*Uc - Uc*Vx - Vc*Vy - g*Hy;
            *(dH + i*ny + j) = -Uc*Hx - Hc*Ux - Vc*Hy - Hc*Vy;
        }
    }
}

void init(float *U, float *V, float *H){
    float Hc, Hd, Hr, Hx, Hy;
    
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            *(H + i*ny + j) = H0 + H1 * tan(9*(j*hy - 6*hx)/(2*D)) + H2 * sin(2 * PI * i * hx)/pow((9*j*hy - 6*hx)/D, 2);
        }
    }

    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            Hc = *(H + i*ny + j);
            Hr = (i==nx-1) ? *(H + j) : *(H + (i+1)*ny + j);
            Hd = (j==ny-1) ? 0 : *(H + i*ny + j+1);

            Hx = (Hr - Hc)/hx;
            Hy = (Hd - Hc)/hy;

            *(U + i*ny + j) = -Hy/f;
            *(V + i*ny + j) = -Hx/f;
        }
    }
}

void writeResult(float *U, float *V, float *H){

}