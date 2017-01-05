#include <iostream>
#include "particlesystem.h"
#include <cstring>

using namespace std;
//1483139279

int main(int argc, char *argv[])
{
    ParticleSystem ps(3, 3, 3);
    //ps.checkInOn();
    for (int i = 1; i < argc; ++i){
        for (int j = 0; j < strlen(argv[i]); ++j){
            if (argv[i][j] == '='){
                if (strncmp(argv[i], "-seed", j) == 0){
                    ps.seed = atoi(argv[i] + j + 1);
               }
                else{
                    if (strncmp(argv[i], "-thr", j) == 0){
                        ps.threads = atoi(argv[i] + j + 1);
                        ps.OMP_ON = ps.threads > 1;
                    }
                }
                break;
            }
        }
    }
    if (ps.seed == -1){
        ps.seed = time(NULL);
    }
    srand(ps.seed);
    cout << "seed = " << ps.seed << "\n";
    cout << "threads = " << ps.threads << "\n";
    cout << "openmp on = " << ps.OMP_ON << "\n";
    ps.print_table_params();
    if (!SKIP_FITTING_B) ps.fittingB();
    if (!SKIP_FITTING_AB) ps.fittingAB();
    ps.fittingA();
    cerr << "seed = " << ps.seed << "\n";
    return 0;
}
