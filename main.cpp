#include <iostream>
#include "particlesystem.h"
#include <cstring>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[])
{
    string filename = "config.cfg";
    if(argc == 2){
        filename = (argv[1]);
    }
    ParticleSystem ps(filename);
    ps.fittingB();
    ps.fittingAB();
    ps.fittingA();
    return 0;
}
