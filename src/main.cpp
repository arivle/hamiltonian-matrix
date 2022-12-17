#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <stdlib.h>
#include <complex>

using namespace std;

vector<complex<double>> kinpropX(vector<complex<double>>& in, int b, int nx, float dt){
    // vector<complex<double>> ps0_out(in.size());
    complex<double> ct (cos(-dt * 1),0);
    complex<double> st (0, (-1 * sin(-dt*1)));
    int nrx = nx/2;
    int sx = 0;
    int sx1 = 0;

    for (int i = 0; i < nx; i++)
    {
        if (b==0)
        {
            for (int j=0; j<nrx; j++)
            {
                sx = 2*j*b+(nx*i);
                sx1 = sx+1;

                in[sx] = ct * in[sx] - st * in[sx1];
                in[sx1] = -st * in[sx] + ct * in[sx1];
            }
        } 
        else if (b==1)
        {
            int nm2 = nx *i;
            int nm3 = nx - 1;

            in[nm2] = ct * in[nm2] - st * in[nm3];
            in[nm3] = -st * in[nm2] + ct * in[nm3];

            for (int j = 0; j < nrx-1; j++)
            {
                sx = 2 * j + b + (nx*i);
                sx1 = sx+1;

                in[sx] = ct * in[sx] - st*in[sx1];
                in[sx1] = -st * in[sx] + ct * in[sx1];
            }
            
        }  
    }
    return in;


}

vector<complex<double>> kinpropY(vector<complex<double>>& in, int b, int nx, int ny, float dt){
    // vector<complex<double>> ps0_out(in.size());
    complex<double> ct (cos(-dt * 1),0);
    complex<double> st (0, (-1 * sin(-dt*1)));
    int nrx = nx/2;
    int nry = ny/2;
    int sy = 0;
    int sy1 = 0;
    
    if (b==0)
    {
        for (int i = 0; i < nrx; i++)
        {
            for (int j = 0; j < nry; j++)
            {
                sy = (2*j*ny) + (2*i);
                sy1 = sy+ny;

                in[sy] = ct * in[sy] - st * in[sy1];
                in[sy1] = -st * in[sy] + ct * in[sy1];
            }
            
        }
        
    }
    else if (b==1)
    {
        for (int i = 0; i < nrx; i++)
        {
            int nm4 = nry-1;
            int nm5 = 2*i +1;
            int nm6 = (2*i + (2 * nm4 * ny)+ b +ny);

            in[nm5] = ct * in[nm5] - st * in[nm6];
            in[nm6] = -st * in[nm5] + ct * in[nm6];

            for (int j = 0; j < nm4; j++)
            {
                sy = (2*j*ny) + b + ny + (2*i);
                sy1 = sy + ny;

                in[sy] = ct * in[sy] - st * in[sy1];
                in[sy1] = -st * in[sy] + ct * in[sy1];
            }
            

        }
        
    }
    
    
    return in;


}

int main(int argc, char const *argv[])
{
    /*variable declaration*/
    int Lx = 1;
    int Ly = 1;
    /*for shorten my computation time. please change the nx and ny*/
    int nx = 4;
    int ny = 4;
    // int nm = nx/2;
    int nmat = nx*nx;
    float dx = Lx/nx;
    float dy = Ly/ny;
    float dt = 0.02;
    float h = 0.01;
    /*just make my computational time shorter. please change this nstep as you need*/
    int nstep = 10; 
    float ta = (nstep/2)*dt;
    float tau = 0.0;
    int jj = 0;
    int nta = round(nstep/2);
    vector<complex<double> > corrF(nstep);
   
    /*vector of initial wave. the vector contains random integer of -1 to 1*/
    vector<complex<double> > psi0(nmat);
    complex<double> I = complex<double>(-1.0, 1.0);
    for (int i =0; i<nmat; i++)
    {
        // x[i] = ((double)rand()/RAND_MAX)-1;
        psi0[i] = ((double)rand()/RAND_MAX)-1 + I*(((double)rand()/RAND_MAX)-1);
    }

    /*conjugate the vector complex psi0*/
    vector<complex<double> > psi0_conj(nmat);
    for (int i =0; i<nmat; i++)
    {
        psi0_conj[i]=conj(psi0[i]);
    }
    
    //this loop just for showing the array contents
    cout<<"psi0: ";
    for (int i=0; i<nmat; i++){
        // cout<<"ini x "<<x[i]<<" ";
        cout<<psi0[i]<< " ";
    }

    cout<<"psi0_conj: ";
    for (int i=0; i<nmat; i++){
        // cout<<"ini x "<<x[i]<<" ";
        cout<<psi0_conj[i]<< " ";
    }

    /*calculating the kinetics*/
    // while (tau<ta)
    // {
    //     for (int kk = 0; kk < nx; kk++)
    //     {
    //         /* code */
    //     }
        
    // }
    return 0;
}
