#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#include <vector>
#include <cmath>
#include <string>
#include "particle.h"

using std::vector;
using std::string;


class ParticleSystem
{
private:
    int printErrorFlag, printConfigurationFlag;
    double tableLatticeConstant, tableCohesiveEnergy, tableB,
            tableC11, tableC12, tableC44, tableEsol, tableEin, tableEon, cohEnergyPurityAtom;
    double bb_a0[2], bb_a1[2], bb_ksi[2], bb_p[2], bb_q[2], bb_r0[2];
    double ab_a0[2], ab_a1[2], ab_ksi[2], ab_p[2], ab_q[2], ab_r0[2];
    double aa_a0[2], aa_a1[2], aa_ksi[2], aa_p[2], aa_q[2], aa_r0[2];
    int OMP_ON;
    int seed, threads;
    string nameAtom, namePurityAtom;
    vector <Particle> particles;
    int xNumber, yNumber, zNumber;
    int isBulk;
    double latticeConstant, cutoff;
    double r0[3], A1[3], A0[3], p[3], q[3], ksi[3];
    double V0, C11, C12, C44, B, Esol, cohEnergy, energyIn, energyOn;

    double euclideanNorm(Vector3d);
    double distance(int i, int j, int, double);
    double distance(Vector3d i, Vector3d j);
    bool get_out(Vector3d &);
    bool between(double, double, double);
    void generate_params(int type, vector<double> &);
    void load_params(int type, vector<double> &);
    void print_params(vector<double> &);
    void setV0();
    void setC11_C12();
    void setC44();
    void setB();
    void setEsol();
    double error_b();
    double error_ab();
    double error_a();
    void setEnergyIn();
    void setEnergyOn();
    bool stop_rule(vector<vector<double>>&, vector<double>&);
    int type(int i, int j);
    void printConfiguration();
public:
    void fittingB();
    void fittingAB();
    void fittingA();
    void print_table_params();
    void set_isBulk(int flag){isBulk = flag;}
    void print_params(int);
    double NMA(double (ParticleSystem::*error)(), int type, double alpha=1.0, double beta=0.5, double gamma=2.0);
    void set_r0(int type, double alpha) {r0[type] = latticeConstant * (1.0 + 0) / sqrt(2.0);}
    void set_cutoff(double alpha) {cutoff = 1.5 * latticeConstant * (1.0 + 0);}
    ParticleSystem(string &);
    double cohesiveEnergy(int type = 0, double alpha = 0.0);
    double cohesiveEnergy_omp(int type = 0, double alpha = 0.0);
    void printV0();
    void printC11_C12();
    void printC44();
    void printB();
    void printEsol();
};

#endif // PARTICLESYSTEM_H
