#include <iostream>
#include <string>
#include <vector>
#include <random>

using namespace std;


int main(int argc, char const *argv[])
{
    /*variable declaration*/
    int Lx = 1;
    int Ly = 1;
    /*for shorten my computation. please change the nx and ny*/
    int nx = 4;
    int ny = 4;
    float nm = nx/2;
    int nmat = nx*nx;
    float dx = Lx/nx;
    float dy = Ly/ny;
    float dt = 0.02;
    float h = 0.01;
    /*just make my computational time shorter. please change this nstep as you need*/
    int nstep = 10; 
    float ta = (nstep/2)*dt;

    float linex[nstep];
    for (int i=0; i<nstep; i++){
        linex[i] = -ta+(i*dt);
    }
    /*this loop just for showing the array contents 
    
    for (int i=0; i<nstep; i++){
        cout<<linex[i]<<" ";
    }

    */
    
    /*vector of initial wave. the vector contains random integer of -1 to 1*/
    vector<float> psio(nmat);
    for (int i =0; i<nmat; i++){
        psio[i] = random () % 3 + (-1);
    }
    
    /*this loop just for showing the array contents
    for (int i=0; i<nmat; i++){
        cout<<psio[i]<<" ";
    }
    */
    int jj = 0;
    float tau = 0.0;
    int corrF[nstep];
    int nta=round(nstep/2);

    /*calculating the kinetics*/
    while (tau<ta)
    {
        for (int kk = 0; kk < nx; kk++)
        {
            /* code */
        }
        
    }
    
    






    return 0;
}
