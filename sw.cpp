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

void calDerivateU(float* U, float *dU, float* V, float* H); // ham tinh dao ham cua U theo t
void calDerivateV(float* V, float *dV, float* U, float* H); // ham tinh dao ham cua V theo t
void calDerivateH(float* H, float *dH, float* U, float* V); // ham tinh dao ham cua H theo t

void init(float *U, float *V, float *H); // ham khoi tao

void writeResult(float *U, float *V, float *H); // ham ghi ket qua ra file

int nx = L/hx;
int ny = L/hy;

int main(){
    float *U, *V, *H;
    float *dU, *dV, *dH;
    float t = 0;

    init(U, V, H);
    writeResult(U, V, H);

    while (t <= T)
    {
        calDerivateU(U, dU, V, H);
        calDerivateV(V, dV, U, H);
        calDerivateH(H, dH, U, V);

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

void calDerivateU(float* U, float *dU, float* V, float* H){
    float Ux, Uy, Ur, Ud, Uc, Vc, Hr, Hc, Hx;
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            Uc = *(U + i*ny + j);
            Ur = (i==nx-1) ? *(U + j) : *(U + (i+1)*ny + j);
            Ud = (j==ny-1) ? 0 : *(U + i*ny + j+1);

            Ux = (Ur - Uc)/hx;
            Uy = (Ud - Uc)/hy;

            Vc = *(V + i*ny + j);

            Hr = (i==nx-1) ? *(H + j) : *(H + (i+1)*ny + j);
            Hc = *(H + i*ny + j);

            Hx = (Hr - Hc)/hx;

            *(dU + i*ny + j) = f*Vc - Uc*Ux - Vc*Uy - g*Hx;
        }
    }
}

void calDerivateV(float* V, float *dV, float* U, float* H){
    float Vx, Vy, Vr, Vd, Vc, Uc, Hd, Hc, Hy;
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            Vc = *(V + i*ny + j);
            Vr = (i==nx-1) ? *(V + j) : *(V + (i+1)*ny + j);
            Vd = (j==ny-1) ? 0 : *(V + i*ny + j+1);

            Vx = (Vr - Vc)/hx;
            Vy = (Vd - Vc)/hy;

            Uc = *(U + i*ny + j);
            
            Hd = (j==ny-1) ? 0 : *(H + i*ny + j+1);
            Hc = *(H + i*ny + j);

            Hy = (Hd - Hc)/hy;

            *(dV + i*ny + j) = -f*Uc - Uc*Vx - Vc*Vy -g*Hy;
        }
    }
}

void calDerivateH(float* H, float *dH, float* U, float* V){
    float Hc, Hd, Hr, Hx, Hy, Uc, Ur, Ux, Vc, Vd, Vy;

    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            Hc = *(H + i*ny + j);
            Hr = (i==nx-1) ? *(H + j) : *(H + (i+1)*ny + j);
            Hd = (j==ny-1) ? 0 : *(H + i*ny + j+1);

            Hx = (Hr - Hc)/hx;
            Hy = (Hd - Hc)/hy;

            Uc = *(U + i*ny + j);
            Ur = (i==nx-1) ? *(U + j) : *(U + (i+1)*ny + j);

            Ux = (Ur - Uc)/hx;

            Vc = *(V + i*ny + j);
            Vd = (j==ny-1) ? 0 : *(V + i*ny + j+1);

            Vy = (Vd - Vc)/hy;

            *(dH + i*ny + j) = -Uc*Hx -Hc*Ux -Vc*Hy -Hc*Vy;
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