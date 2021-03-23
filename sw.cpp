#include<iostream>

#define g = 9.8
#define f = 0.1
#define L = 10
#define hx = 1
#define hy = 1

void calDerivateU(float* U, float *dU); // ham tinh dao ham cua U theo t
void calDerivateV(float* V, float *dV); // ham tinh dao ham cua V theo t
void calDerivateH(float* H, float *dH); // ham tinh dao ham cua H theo t

void init(float *C); // ham khoi tao

void writeResult(float *C); // ham ghi ket qua ra file

int main(){
    float* U, V, H;
    float* dU, dV, dH;


    
}