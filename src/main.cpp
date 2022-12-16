#include <iostream>
#include <string>

using namespace std;


int main(int argc, char const *argv[])
{
    /*deklarasi variabel yang digunakan*/
    int Lx = 1;
    int Ly = 1;
    int nx = 16;
    int ny = 16;
    float nm = nx/2;
    int nmat = nx*nx;
    float dx = Lx/nx;
    float dy = Ly/ny;
    float dt = 0.02;
    float h = 0.01;
    int nstep = 8192;
    float ta = (nstep/2)*dt;
    return 0;
}
