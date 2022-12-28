#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <stdlib.h>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iterator>
#include<cmath>

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

complex<double> dot_product(vector<complex<double>> &a, vector<complex<double>> &b, int nmat){
    complex<double> product = 0.0;
    for (int i = 0; i < nmat; i++)
    {
        product = product + a[i] * b[i];

    }
    return product;
    

}

complex<double> normalize_psi(vector<complex<double>> &a, int num, int nmat){
    complex<double> sum0n = 0.0;
    complex<double> sum0(nmat);

    for (int i = 0; i < nmat; i++)
    {
        sum0n = a[i] * a[i];
        if (num==0)
        {
            sum0[i] = sqrt(4*sum0n[i]);
        }
        else if (num==1)
        {
            sum0[i] = sqrt(2*sum0n[i]);
        }
   
    }
    return sum0;    
}

vector<complex<double>> fermi(int ff, int c, vector<complex<double>>& in, vector<complex<double>>& fermiDirac, vector<complex<double>>& fermiDirac_t1, complex<double> psi0N, int nx, int ny, float dt, int nm){
    int T = 300;
    double Kb = 0.695;
    int pot = 0;
    int mu = 0;
    int nrx = nx/2;
    int nry = ny/2;
    int sx1 = 0;
    int sy1 = 0;
    int sx2 = 0;
    int sy2 = 0;

    double TEn = T*Kb;
    complex<double> potN = pot / psi0N;
    complex<double> bN = psi0N / TEn;
    complex<double> aN = 1/psi0N;

    complex<double> ct (cos(-dt * 1),0);
    complex<double> st (0, (-1 * sin(-dt*1)));

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (i==j)
            {
                fermiDirac[i][j] = ct * sqrt((bN*potN) - (mu/TEn));
            }
            
        }
        
    }

    // filling the inverse matrix of fermi dirac: 2x2 block diagonal parth with even index
    for (int i = 0; i < nrx; i++)
    {
        for (int j = 0; j < nry; j++)
        {
            if (i == j)
            {
                sx1 = 2*i;
                sy1 = 2*j;
                fermiDirac[sx1][sy1+1] = -st * sqrt((bN*potN) - (mu/TEn));
                fermiDirac[sx1+1][sy1] = -st * sqrt((bN*potN) - (mu/TEn));
            }
            
        }
        
    }

    // filling the invers matrix of fermi dirac: 2x2 block diagonal with odd index
    for (int j = 0; j < ny; i++)
    {
        sx2 = 2*0+0;
        sy2 = 2*j+0;
        fermiDirac_t1[sx2][sy2] = ct*fermiDirac[sx2][sy2] - st*fermiDirac[nm][sy2];
        fermiDirac_t1[nm][sy2] = -st*fermiDirac[sx2][sy2] + ct*fermiDirac[nm][sy2];
        sx2 = 0;
        sy2 = 0;

    }

    for (int i = 0; i < nrx-1; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            sx2 = 2*i+1;
            sy2 = j;
            fermiDirac_t1[sx2][sy2] = ct*fermiDirac[sx2][sy2] - st*fermiDirac[sx2+1][sy2];
            fermiDirac_t1[sx2+1][sy] = -ct*fermiDirac[sx2][sy2] + ct*fermiDirac[sx2+1][sy2];
        }
        
    }
    
    
    

    

}

void write_txt(string filename, vector<complex<double>> vect, int nmat){
    ofstream myFile(filename);
    // myFile.open(filename, ios::out | ios::binary);

    // myFile << colname<<"\n";

    for (int i = 0; i < nmat; i++)
    {
        myFile<<setprecision(16)<<vect[i]<<"\n";
        // myFile.write(vect[i], sizeof(float));
    }
    myFile.close();   

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
    double dx = Lx/nx;
    double dy = Ly/ny;
    float dt = 0.02;
    float h = 0.01;
    /*just make my computational time shorter. please change this nstep as you need*/
    int nstep = 10; 
    int nm = nstep-1;
    float ta = (nstep/2)*dt;
    float tau = 0.0;
    int jj = 0;
    int nta = round(nstep/2);
    vector<complex<double> > corrF(nstep);
   
    /*vector of initial wave. the vector contains random integer of -1 to 1*/
    vector<complex<double> > psi0(nmat);

    /*vector of normalization psi*/
    vector<complex<double> > psi0N(nmat);

    /*vector of zeros for fermi function*/
    vector<complex<double> > fermiDirac(nmat);
    vector<complex<double> > fermiDirac_t1(nmat);

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

    psi0N = normalize_psi(psi0, 0, nmat);
    psi0N = normalize_psi(psi0, 1, nmat);
    
    //this loop just for showing the array contents
    // cout<<"psi0: ";
    // for (int i=0; i<nmat; i++){
    //     // cout<<"ini x "<<x[i]<<" ";
    //     cout<<psi0[i]<< " ";
    // }

    // cout<<"psi0_conj: ";
    // for (int i=0; i<nmat; i++){
    //     // cout<<"ini x "<<x[i]<<" ";
    //     cout<<psi0_conj[i]<< " ";
    // }

    /*calculating the kinetics*/
    while (tau<ta)
    {
        psi0 = kinpropX(psi0,0,nx,dt);
        psi0 = kinpropX(psi0,1,nx,dt);
        psi0 = kinpropY(psi0, 0,nx, ny, dt);
        psi0 = kinpropY(psi0, 1,nx, ny, dt);

        complex<double> cf0 = dot_product(psi0_conj, psi0, nmat);
        complex<double> cf1 = cf0 * dx * dy;
        corrF[jj + nta] = cf1;

        jj+=1;
        tau = dt * jj;

    }
    complex<double> matrix_fermi = fermi(1, 0, corrF, fermiDirac, fermiDirac_t1, psi0N, nx, ny, dt, nm);
    write_txt("corrf.txt", corrF, nmat);
    
    return 0;
}
