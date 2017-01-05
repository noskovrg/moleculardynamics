#include "particlesystem.h"
#include <iostream>
#include <algorithm>
#include <omp.h>

using std::cerr;
using std::cout;
using std::endl;
using std::pair;
const double eps_b = 3e-2;
const double eps_ab = 3e-4; // -2
const double eps_a = 3e-2;

double ParticleSystem::euclideanNorm(Vector3d r)
{
    return sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
}


ParticleSystem::ParticleSystem(int _xNumber, int _yNumber, int _zNumber)
{
    OMP_ON = 0;
    seed = -1;
    threads = 1;
    isBulk = 1;
    if (SKIP_FITTING_B){
        A0[0] = 0.12051;
        A1[0] = 0.0320955;
        ksi[0] = 1.21461;
        p[0] = 10.6915;
        q[0] = 2.94819;
        r0[0] = 2.8704;
        latticeConstant = 4.05936;
        cohEnergy = -2.90166;
    }
    if (SKIP_FITTING_AB){
        A0[1] = 0.0773566;
        A1[1] = -0.00913184;
        ksi[1] = 0.851237;
        p[1] = 15.574;
        q[1] = 2.41335;
        r0[1] = 2.88853;
        Esol = -0.578009;
    }
    xNumber = _xNumber; yNumber = _yNumber; zNumber = _zNumber;
    if (isBulk){
        for (int itX = 0; itX < xNumber; ++itX){
            for (int itY = 0; itY < yNumber; ++itY){
                for (int itZ = 0; itZ < zNumber;  ++itZ){
                    Particle tmpParticle0(itX, itY, itZ, 0);
                    particles.push_back(tmpParticle0);
                    //if (itX == xNumber || itY == yNumber || itZ == zNumber) continue;
                    Particle tmpParticle1(itX + 0.5, itY + 0.5, itZ, 0);
                    Particle tmpParticle2(itX + 0.5, itY, itZ + 0.5, 0);
                    Particle tmpParticle3(itX, itY + 0.5, itZ + 0.5, 0);
                    particles.push_back(tmpParticle1);
                    particles.push_back(tmpParticle2);
                    particles.push_back(tmpParticle3);
                }
            }
        }
    }
    else{
        for (int itX = 0; itX <= xNumber; ++itX){
            for (int itY = 0; itY <= yNumber; ++itY){
                for (int itZ = 0; itZ <= zNumber;  ++itZ){
                    Particle tmpParticle0(itX, itY, itZ, 0);
                    particles.push_back(tmpParticle0);
                    if (itX == xNumber || itY == yNumber || itZ == zNumber) continue;
                    Particle tmpParticle1(itX + 0.5, itY + 0.5, itZ, 0);
                    Particle tmpParticle2(itX + 0.5, itY, itZ + 0.5, 0);
                    Particle tmpParticle3(itX, itY + 0.5, itZ + 0.5, 0);
                    particles.push_back(tmpParticle1);
                    particles.push_back(tmpParticle2);
                    particles.push_back(tmpParticle3);
                }
            }
        }
    }
    cerr << "number of particles: " << particles.size() << endl;
}


double ParticleSystem::cohesiveEnergy_omp(int typeTransform, double alpha)
{
    cerr << "ecoh omp in\n";
    cutoff = 1.7 * latticeConstant;
    double energy = 0.0;
    vector<double> repulsiveEnergy(threads, 0.0);
    vector<double> bandEnergy(threads, 0.0);
    vector<double> tmpEnergy(threads, 0.0);
    #pragma omp parallel for num_threads(threads)

    for (int i = 0; i < particles.size(); ++i){
        int thr = omp_get_thread_num();
        repulsiveEnergy[thr] = bandEnergy[thr] = 0.0;
        Vector3d pos_i = particles[i].getPosition(latticeConstant, typeTransform, alpha);
        for (int j = 0; j < particles.size(); ++j){
           if (i == j) continue;
           for (int dx : {-1, 0, 1}){
                for (int dy : {-1, 0, 1}){
                    if (isBulk == 1){
                        for (int dz : {-1, 0, 1}){
                            Vector3d pos_j = particles[j].getPosition(latticeConstant, typeTransform, alpha,
                                                                      dx * xNumber, dy * yNumber, dz * zNumber);
                            double dist = distance(pos_i, pos_j);
                            if (dist < cutoff) {
                                repulsiveEnergy[thr] += (A1[type(i, j)] * (dist - r0[type(i, j)]) + A0[type(i, j)]) * exp(-p[type(i, j)] * (dist / r0[type(i, j)] - 1.0));
                                bandEnergy[thr] += ksi[type(i, j)] * ksi[type(i, j)] * exp(-2.0 * q[type(i, j)] * ( dist / r0[type(i, j)] - 1.0));
                            }
                        }
                    }
                    else{
                        for (int dz : {0}){
                            Vector3d pos_j = particles[j].getPosition(latticeConstant, typeTransform, alpha,
                                                                      dx * xNumber, dy * yNumber, dz * zNumber);
                            double dist = distance(pos_i, pos_j);
                            if (dist < cutoff) {
                                repulsiveEnergy[thr] += (A1[type(i, j)] * (dist - r0[type(i, j)]) + A0[type(i, j)]) * exp(-p[type(i, j)] * (dist / r0[type(i, j)] - 1.0));
                                bandEnergy[thr] += ksi[type(i, j)] * ksi[type(i, j)] * exp(-2.0 * q[type(i, j)] * ( dist / r0[type(i, j)] - 1.0));
                            }
                        }
                    }
                }
            }
        }
        bandEnergy[thr] = -sqrt(bandEnergy[thr]);
        tmpEnergy[thr] += repulsiveEnergy[thr] + bandEnergy[thr];
       // cerr << repulsiveEnergy << "\t" << bandEnergy << "\n";
    }
    for (int i = 0; i < threads; ++i){
        energy += tmpEnergy[i];
    }
   // cerr << "total energy = " << energy << endl;
   // cohEnergy =  energy / particles.size();
    return energy / particles.size();
}



double ParticleSystem::cohesiveEnergy(int typeTransform, double alpha)
{
    //set_r0(alpha);
    //set_cutoff(alpha);
    cutoff = 1.7 * latticeConstant;
    double energy = 0.0;
    double repulsiveEnergy = 0.0;
    double bandEnergy = 0.0;

    for (int i = 0; i < particles.size(); ++i){
       // cerr << "i = " << i << "\n";
        repulsiveEnergy = bandEnergy = 0.0;
        Vector3d pos_i = particles[i].getPosition(latticeConstant, typeTransform, alpha);
        for (int j = 0; j < particles.size(); ++j){
           if (i == j) continue;
           for (int dx : {-1, 0, 1}){
                for (int dy : {-1, 0, 1}){
                    if (isBulk == 1){
                        for (int dz : {-1, 0, 1}){
                            Vector3d pos_j = particles[j].getPosition(latticeConstant, typeTransform, alpha,
                                                                      dx * xNumber, dy * yNumber, dz * zNumber);
                            double dist = distance(pos_i, pos_j);
                            if (dist < cutoff) {
                                repulsiveEnergy += (A1[type(i, j)] * (dist - r0[type(i, j)]) + A0[type(i, j)]) * exp(-p[type(i, j)] * (dist / r0[type(i, j)] - 1.0));
                                bandEnergy += ksi[type(i, j)] * ksi[type(i, j)] * exp(-2.0 * q[type(i, j)] * ( dist / r0[type(i, j)] - 1.0));
                            }
                        }
                    }
                    else{
                        for (int dz : {0}){
                            Vector3d pos_j = particles[j].getPosition(latticeConstant, typeTransform, alpha,
                                                                      dx * xNumber, dy * yNumber, dz * zNumber);
                            double dist = distance(pos_i, pos_j);
                            if (dist < cutoff) {
                                repulsiveEnergy += (A1[type(i, j)] * (dist - r0[type(i, j)]) + A0[type(i, j)]) * exp(-p[type(i, j)] * (dist / r0[type(i, j)] - 1.0));
                                bandEnergy += ksi[type(i, j)] * ksi[type(i, j)] * exp(-2.0 * q[type(i, j)] * ( dist / r0[type(i, j)] - 1.0));
                            }
                        }
                    }
                }
            }
        }

        bandEnergy = -sqrt(bandEnergy);
        energy += repulsiveEnergy + bandEnergy;
    }
   // cerr << "total energy = " << energy << endl;
   // cohEnergy =  energy / particles.size();
    return energy / particles.size();
}

double ParticleSystem::distance(int i, int j, int typeTransform, double alpha)
{
    Vector3d pos_i = particles[i].getPosition(latticeConstant, typeTransform, alpha);
    Vector3d pos_j = particles[j].getPosition(latticeConstant, typeTransform, alpha);
    return euclideanNorm(pos_i - pos_j);
}

double ParticleSystem::distance(Vector3d i, Vector3d j)
{
    return euclideanNorm(i - j);
}

bool ParticleSystem::get_out(Vector3d & a)
{
    return !((a.x >= 0 && a.x < xNumber * latticeConstant) &&
           (a.y >= 0 && a.y < yNumber * latticeConstant) &&
             (a.z >= 0 && a.z < zNumber * latticeConstant));
}

bool ParticleSystem::between(double x, double l, double r)
{
    return x >= l && x < r;
}

int comp(pair<double, int> a, pair<double, int> b){
    return a.first < b.first;
}

double ParticleSystem::NMA(double (ParticleSystem::*error)(), int _type, double alpha, double beta, double gamma)
{
    //cerr << "NMA start " << PRINT_FLAG << "\n";
    const int DIM = 6;
    vector<vector<double> > simplex(DIM + 1, vector<double> (DIM));
    vector<pair<double, int>> f(DIM + 1);
    int x_h, x_g, x_l;
    double f_r, f_h, f_g, f_l, f_e, f_s, f_c;
    vector<double> x_c(DIM, 0);
    vector<double> x_r(DIM, 0);
    vector<double> x_e(DIM, 0);
    vector<double> x_s(DIM, 0);
    double err = 1e10;
    double globalError = 1e10;
    double eps;
    int cntIt = 0;
    if (_type == 0){
        eps = eps_b;
    }
    else if (_type == 1){
        eps = eps_ab;
    }
    else{
        eps = eps_a;
    }
    // step 1
    int cntSteps = 0;
    while (true){
        if(err < eps){
            cerr << "error: " << err << "\n";
            cerr << "OPTIMIZED\n";
            print_params(_type);
            break;
        }
        for (int i = 0; i < DIM + 1; ++i){
            generate_params(_type, simplex[i]);
            //print_params(simplex[i]);
            cntSteps = 0;
            //cntIt++;
        }
        while (true){
            cntSteps++;
            if (cntSteps > 1000){
                break;
            }
           // cerr << "while start\n";
            for (int i = 0; i < DIM + 1; ++i){
                load_params(_type, simplex[i]);
                f[i].first = (this->*error)();
                f[i].second = i;
            }
            // step 2
            sort(f.begin(), f.end(), comp);
            x_h = f[5].second; load_params(_type, simplex[x_h]); f_h = (this->*error)();
            x_g = f[4].second; load_params(_type, simplex[x_g]); f_g = (this->*error)();
            x_l = f[0].second; load_params(_type, simplex[x_l]); f_l = (this->*error)();
            // step 3
            x_c.assign(DIM, 0);
            for (int i = 0; i < simplex.size(); ++i){
                if (i == x_h) continue;
                for (int dim = 0; dim < DIM; ++dim){
                    x_c[dim] += simplex[i][dim];
                }
            }
            for (int dim = 0; dim < DIM; ++dim){
                x_c[dim] = x_c[dim] / (double)DIM;
            }
            globalError = (globalError < f_l ? globalError : f_l);
           // cerr << "FLAG " << PRINT_FLAG <<  _type << "\n";
            // step 9
            if (PRINT_FLAG){
            	cout << cntIt++ << " " << globalError << "\n";
               /* if (_type == 0) cerr << "error: " << f_l << "\n";
                if (_type == 1) cerr << "error: " << f_l << "\t " << Esol << "\n";
                if (_type == 2){
                    //if (f_l < 5e-3) {
                        cerr << "error: " << f_l << "\t " << energyIn << "\t" << energyOn << "\n";
                        //print_params(_type);
                    //}
                }*/

            }
            // print_params(params[x_l]);
            if (f_l < eps || stop_rule(simplex, x_c)){
                load_params(_type, x_c);
                f_c = (this->*error)();
                if (f_c < f_l){
                    err = f_c;
                }
                else{
                    err = f_l;
                    load_params(_type, simplex[x_l]);
                    if (_type == 1) setEsol();
                    if (_type == 2){
                        setEnergyIn();
                        setEnergyOn();
                    }
                }
                break;
            }
            // step 4
            for (int dim = 0; dim < DIM; ++dim){
                x_r[dim] = (1.0 + alpha) * x_c[dim] - alpha * simplex[x_h][dim];
            }
            load_params(_type, x_r);
            f_r = (this->*error)();
            // step 5
            if (f_r < f_l){
                for (int dim = 0; dim < DIM; ++dim){
                    x_e[dim] = (1.0 - gamma) * x_c[dim] + gamma * x_r[dim];
                }
                load_params(_type, x_e);
                f_e = (this->*error)();
                if (f_e < f_l){
                    simplex[x_h] = x_e;
                    continue;
                }
                else{
                    simplex[x_h] = x_r;
                    continue;
                }
            }
            else if (f_l < f_r && f_r < f_g){
                simplex[x_h] = x_r;
                continue;
            }
            else if (f_g < f_r && f_r < f_h){
                simplex[x_h] = x_r;
            }
            // step 6
            for (int dim = 0; dim < DIM; ++dim){
                x_s[dim] = (1.0 - beta) * x_c[dim] + beta * simplex[x_h][dim];
            }
            load_params(_type, x_s);
            f_s = (this->*error)();
            // step 7
            if (f_s < f_h){
                simplex[x_h] = x_s;
                continue;
            }
            // step 8
            else{
                for (int i = 0; i < simplex.size(); ++i){
                    if (i == x_l) continue;
                    for (int dim = 0; dim < DIM; ++dim){
                        simplex[i][dim] = simplex[x_l][dim] + (simplex[i][dim] - simplex[x_l][dim]) / 2.0;
                    }
                }
            }
        }
    }

    cout << "\n\n";
}

void ParticleSystem::load_params(int type, vector<double> & params)
{
    A0[type] = params[0];
    A1[type] = params[1];
    ksi[type] = params[2];
    p[type] = params[3];
    q[type] = params[4];
    r0[type] = params[5];
    if (type == 0) {
        latticeConstant = r0[type] * sqrt(2);
    }
}

double ParticleSystem::error_b()
{
    cohEnergy = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy());
    setB();
    setC11_C12();
    setC44();
   /* return ((-2.96 - cohEnergy) / -2.96) * ((-2.96 - cohEnergy) / -2.96) +
            ((4.085 - latticeConstant) / 4.085) * ((4.085 - latticeConstant) / 4.085)
            + ((1.08 - B) / 1.08) * ((1.08 - B) / 1.08)
            + ((1.32 - C11) / 1.32) * ((1.32 - C11) / 1.32)
            + ((0.97 - C12) / 0.97) * ((0.97 - C12) / 0.97)
            + ((0.51 - C44) / 0.51) * ((0.51 - C44) / 0.51);*/
    return sqrt((((-2.96 - cohEnergy)) * ((-2.96 - cohEnergy))/pow(-2.96, 2) +
                ((4.085 - latticeConstant) ) * ((4.085 - latticeConstant))/pow(4.085, 2)
                + ((1.08 - B)) * ((1.08 - B))/pow(1.08, 2)
                + ((1.32 - C11)) * ((1.32 - C11))/pow(1.32, 2)
                + ((0.97 - C12)) * ((0.97 - C12))/pow(0.97, 2)
            + ((0.51 - C44)) * ((0.51 - C44))/pow(0.51, 2))/6.0);
}

double ParticleSystem::error_ab()
{
    setEsol();
    double err = sqrt((0.539 - Esol) * (0.539 - Esol)/pow(0.539, 2));
    return err;
}

double ParticleSystem::error_a()
{
    setEnergyIn();
    setEnergyOn();
    double err = sqrt(((0.06 - energyIn) * (0.06 - energyIn)/pow(0.06, 2) +
             (-0.56 - energyOn) * (-0.56 - energyOn)/pow(-0.56, 2))/2.0);
    return err;
}

/*
double ParticleSystem::getEnergyIn()
{
    double res = 0;
    isBulk = 0;
    double energy_surf = cohesiveEnergy() * particles.size();
    particles[94].setType(1);
    double energy_adatom_surf = cohesiveEnergy() * particles.size();
    particles[95].setType(1);
    double energy_dimer_surf = cohesiveEnergy() * particles.size();
    res = (energy_dimer_surf - energy_surf) - 2 * (energy_adatom_surf - energy_surf);
    particles[94].setType(0);
    particles[95].setType(0);
    isBulk = 1;
    energyIn = res;
    return res;
}

double ParticleSystem::getEnergyOn()
{
    double res = 0;
    isBulk = 0;
    double energy_surf = cohesiveEnergy() * particles.size();
    Particle tmp0(0.0, 0.0, 3.0, 1);
    Particle tmp1(0.5, 0.5, 3.0, 1);
    particles.push_back(tmp0);
    double energy_adatom_surf = cohesiveEnergy() * particles.size();
    particles.push_back(tmp1);
    double energy_dimer_surf = cohesiveEnergy() * particles.size();
    particles.pop_back();
    particles.pop_back();
    res = (energy_dimer_surf - energy_surf) - 2 * (energy_adatom_surf - energy_surf);
    isBulk = 1;
    energyOn = res;
    return res;
}
*/


void ParticleSystem::setEnergyIn()
{
    double res = 0;
//    double energy_bulk = cohesiveEnergy() * particles.size();
//    particles[94].setType(1);
//    double energy_bulk_adatom = cohesiveEnergy() * particles.size();
//    particles[95].setType(1);
//    double energy_bulk_dimer = cohesiveEnergy() * particles.size();
//    particles[94].setType(0);
//    particles[95].setType(0)

    isBulk = 0;
    double energy_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    //double energy_surf = (energy_bulk - cohesiveEnergy() * particles.size()) / (18.0 * latticeConstant * latticeConstant);
    particles[94].setType(1);
    double energy_adatom_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    //double energy_adatom_surf = (energy_bulk_adatom - cohesiveEnergy() * particles.size()) / (18.0 * latticeConstant * latticeConstant);
    particles[95].setType(1);
    double energy_dimer_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    //double energy_dimer_surf = (energy_bulk_dimer - cohesiveEnergy() * particles.size()) / (18.0 * latticeConstant * latticeConstant);
    res = (energy_dimer_surf - energy_surf) - 2 * (energy_adatom_surf - energy_surf);
    particles[94].setType(0);
    particles[95].setType(0);
    isBulk = 1;
    energyIn = res;
}

void ParticleSystem::setEnergyOn()
{
    double res = 0;
    Particle tmp0(0.0, 0.0, 3.0, 1);
    Particle tmp1(0.5, 0.5, 3.0, 1);
//    double energy_bulk = cohesiveEnergy() * particles.size();
//    particles.push_back(tmp0);
//    double energy_bulk_adatom = cohesiveEnergy() * particles.size();
//    particles.push_back(tmp1);
//    double energy_bulk_dimer = cohesiveEnergy() * particles.size();
//    particles.pop_back();
//    particles.pop_back();

    isBulk = 0;
    double energy_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    //double energy_surf = (energy_bulk - cohesiveEnergy() * particles.size()) / (18.0 * latticeConstant * latticeConstant);
    particles.push_back(tmp0);
    double energy_adatom_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    //double energy_adatom_surf = (energy_bulk_adatom - cohesiveEnergy() * particles.size()) / (18.0 * latticeConstant * latticeConstant);
    particles.push_back(tmp1);
//    cout << particles.size() << "\n\n";
//    for (int i = 0; i < particles.size(); ++i){
//        Vector3d pos = particles[i].getPosition(4.085, 0, 0);
//        cout << particles[i].getType() << " " << pos.x << " " << pos.y << " " << pos.z << "\n";
//    }
    double energy_dimer_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    //double energy_dimer_surf = (energy_bulk_dimer - cohesiveEnergy() * particles.size()) / (18.0 * latticeConstant * latticeConstant);
    particles.pop_back();
    particles.pop_back();
    res = (energy_dimer_surf - energy_surf) - 2 * (energy_adatom_surf - energy_surf);
    isBulk = 1;
    energyOn = res;
}



void ParticleSystem::print_params(vector<double> & params)
{
    cout << "A0 = " << params[0] << "\nA1 = " << params[1] << "\nksi = " << params[2] <<
            "\np = " << params[3] << "\nq = " << params[4] << "\nr0 = " << params[5] << "\n";
}

void ParticleSystem::print_params(int type)
{
    //set_r0(type, 0);
    cerr << "A0 = " << A0[type] << "\nA1 = " << A1[type] << "\nksi = " << ksi[type] <<
            "\np = " << p[type] << "\nq = " << q[type] << "\nr0 = " << r0[type] << "\n";
}

void ParticleSystem::setV0()
{
    V0 = latticeConstant * latticeConstant * latticeConstant * xNumber * yNumber * zNumber;
}

void ParticleSystem::printV0()
{
    setV0();
    cerr << "V0 = " << V0 << endl;
}

void ParticleSystem::setC11_C12()
{
    C11 = 0;
    C12 = 0;
    {
    double energy_i1 = (OMP_ON ? cohesiveEnergy_omp(11, 0.0): cohesiveEnergy()) * particles.size();
    double energy_i0 = (OMP_ON ? cohesiveEnergy_omp(11, -0.01): cohesiveEnergy(11, -0.01)) * particles.size();
    double energy_i2 = (OMP_ON ? cohesiveEnergy_omp(11, 0.01): cohesiveEnergy(11, 0.01)) * particles.size();
   // cerr << energy_i0 << "\t" << energy_i1 << "\t" << energy_i2 << endl;
    C11 += 1. / V0 * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001 * 0.8018993929636421;
    C12 += 1. / V0 * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001 * 0.8018993929636421;
    }
    {
    double energy_i1 = (OMP_ON ? cohesiveEnergy_omp(12, 0.0): cohesiveEnergy(12, 0.0)) * particles.size();
    double energy_i0 = (OMP_ON ? cohesiveEnergy_omp(12, -0.01): cohesiveEnergy(12, -0.01)) * particles.size();
    double energy_i2 = (OMP_ON ? cohesiveEnergy_omp(12, 0.01): cohesiveEnergy(12, 0.01)) * particles.size();
   // cerr << energy_i0 << "\t" << energy_i1 << "\t" << energy_i2 << endl;
    C11 += 1. / V0 * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001 * 0.8018993929636421;
    C12 -= 1. / V0 * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001 * 0.8018993929636421;
    }
    C11 /= 2;
    C12 /= 2;
}

void ParticleSystem::printC11_C12()
{
    setC11_C12();
    cerr << "C11 = " << C11 << endl;
    cerr << "C12 = " << C12 << endl;
}


void ParticleSystem::setC44()
{
   // latticeConstant = 4.085;
    double energy_i1 = (OMP_ON ? cohesiveEnergy_omp(44, 0.0): cohesiveEnergy(44, 0.0)) * particles.size();
    double energy_i0 = (OMP_ON ? cohesiveEnergy_omp(44, -0.01): cohesiveEnergy(44, -0.01)) * particles.size();
    double energy_i2 = (OMP_ON ? cohesiveEnergy_omp(44, 0.01): cohesiveEnergy(44, 0.01)) * particles.size();
    C44 = 1.0 / (2.0 * V0) * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001 * 0.8018993929636421;
}

void ParticleSystem::printC44()
{
    setC44();
    cerr << "C44 = " << C44 << endl;
}

void ParticleSystem::setB()
{
   // set_r0();
    double energy_i1 = (OMP_ON ? cohesiveEnergy_omp(0, 0.0): cohesiveEnergy(0, 0.0)) * particles.size();
    double energy_i0 = (OMP_ON ? cohesiveEnergy_omp(0, -0.01): cohesiveEnergy(0, -0.01)) * particles.size();
    double energy_i2 = (OMP_ON ? cohesiveEnergy_omp(0, 0.01): cohesiveEnergy(0, 0.01)) * particles.size();
    //cerr << energy_i0 << "\t" << energy_i1 << "\t" << energy_i2 << "\t" << V0 << endl;
    B = 2.0 / (9.0 * V0) * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001 * 0.8018993929636421;
}


void ParticleSystem::setEsol()
{
    double cohEnergy_b = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy());
    //cerr << "coh b = " << cohEnergy_b << "\n";
    double cohEnergy_a = -4.435;
    double energy_b = cohEnergy_b * particles.size();
    particles[0].setType(1);
    double energy_ab = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
   // cerr << "ener " << cohesiveEnergy() << "\n";
    Esol = energy_ab - energy_b - cohEnergy_a + cohEnergy_b;
    particles[0].setType(0);
}

bool ParticleSystem::stop_rule(vector<vector<double>> & simplex, vector<double> & c)
{
    for (int i = 0; i < simplex.size(); ++i){
        double dist = 0;
        for (int dim = 0; dim < c.size(); ++dim){
            dist += (c[dim] - simplex[i][dim]) * (c[dim] - simplex[i][dim]);
        }
        dist = sqrt(dist);
        if (dist > 1e-3){
            return false;
        }
    }
   // cerr << "stop rule\n";
    return true;
}

int ParticleSystem::type(int i, int j)
{
    return particles[i].getType() + particles[j].getType();
}

void ParticleSystem::fittingB()
{
    cerr << "\n\nValues by fitted params for BB type interconection of elements\n";
    NMA(&ParticleSystem::error_b, 0);
    cerr << "lattice constant = " << latticeConstant << "\n";
    cohEnergy = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy());
    cerr <<  "Cohesive energy: " << cohEnergy << endl;
    //printV0();
    setV0();
    printB();
    printC11_C12();
    printC44();
}

void ParticleSystem::fittingAB()
{
    cerr << "\n\nValues by fitted params for AB type interconection of elements\n";
   // particles[53].setType(1);
    NMA(&ParticleSystem::error_ab, 1);
    //setEsol();
    cerr << "Esol = " << Esol << "\n";
   // particles[53].setType(0);
}

void ParticleSystem::fittingA()
{
    cerr << "\n\nValues by fitted params for AA type interconection of elements\n";
    NMA(&ParticleSystem::error_a, 2);
    cerr << "Energy dim in = " << energyIn << "\n";
    cerr << "Energy dim out = " << energyOn << "\n";
}


void ParticleSystem::print_table_params()
{
    cerr << "Values by table params\n";
    setTableParams(0);
    cohEnergy = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy());
    print_params(0);
    cerr << "lattice constant = " << latticeConstant  << "\n";
    cerr << "Cohesive energy: " << cohEnergy << endl;
    //printV0();
    setV0();
    printB();
    printC11_C12();
    printC44();
}

void ParticleSystem::printEsol()
{
    setEsol();
    cerr << "Esol = " << Esol << endl;
}

void ParticleSystem::printB()
{
    setB();
    cerr << "B = " << B << endl;
}

double ParticleSystem::setTableParams(int type)
{
     A1[type] = 0.0; A0[type] = 0.1028; ksi[type] = 1.178; p[type] = 10.928; q[type] = 3.139; latticeConstant = 4.085; r0[type] = latticeConstant / sqrt(2); // Ag
}


void ParticleSystem::generate_params(int _type, vector<double> & params)
{
    if (_type == 0){
        params[0] = random() / (double)RAND_MAX * 0.1028 * 2.0 / 3.0 + 0.1028 * 2.0 / 3.0;
        params[1] = random() / (double)RAND_MAX * 0.2 - 0.1;
        params[2] = random() / (double)RAND_MAX * 1.178 * 2.0 / 3.0 + 1.178 * 2.0 / 3.0;
        params[3] = random() / (double)RAND_MAX * 10.928 * 2.0 / 3.0 + 10.928 * 2.0 / 3.0;
        params[4] = random() / (double)RAND_MAX * 3.139 * 2.0 / 3.0 + 3.139 * 2.0 / 3.0;
        params[5] = random() / (double)RAND_MAX * 2.88853 * 2.0 / 3.0 + 2.88853 * 2.0 / 3.0;
    }
    else if (_type == 1){
        params[0] = random() / (double)RAND_MAX * 0.1028 * 2.0 / 3.0 + 0.1028 * 2.0 / 3.0;
        params[1] =random() / (double)RAND_MAX * 0.2 - 0.1;
        params[2] = random() / (double)RAND_MAX * 1.178 * 2.0 / 3.0 + 1.178 * 2.0 / 3.0;
        params[3] = random() / (double)RAND_MAX * 10.928 * 2.0 / 3.0 + 10.928 * 2.0 / 3.0;
        params[4] = random() / (double)RAND_MAX * 3.139 * 2.0 / 3.0 + 3.139 * 2.0 / 3.0;
        params[5] = random() / (double)RAND_MAX * 2.88853 * 2.0 / 3.0 + 2.88853 * 2.0 / 3.0;

    }
    else{
        params[0] = random() / (double)RAND_MAX * 2 - 1;
        params[1] = random() / (double)RAND_MAX * 2 - 1;
        params[2] = random() / (double)RAND_MAX * 2 - 1;
        params[3] = random() / (double)RAND_MAX * 20 - 10;
        params[4] = random() / (double)RAND_MAX * 2 - 1;
        params[5] = random() / (double)RAND_MAX * 5 - 2.5;

        params[0] = random() / (double)RAND_MAX * 1.06684 * 2.0 / 3.0 + 1.06684 * 2.0 / 3.0;
        params[1] = random() / (double)RAND_MAX * -0.897184 * 2.0 / 3.0 + -0.897184 * 2.0 / 3.0;
        params[2] = random() / (double)RAND_MAX * -1.8706 * 2.0 / 3.0 + -1.8706 * 2.0 / 3.0;
        params[3] = random() / (double)RAND_MAX * 1.63617 * 2.0 / 3.0 + 1.63617 * 2.0 / 3.0;
        params[4] = random() / (double)RAND_MAX * -0.0985134 * 2.0 / 3.0 + -0.0985134 * 2.0 / 3.0;
        params[5] = random() / (double)RAND_MAX * 2.89092 * 2.0 / 3.0 + 2.89092 * 2.0 / 3.0;

    }
}


void ParticleSystem::checkInOn()
{
   // A1[0] = 0.0; A0[0] = 0.1028; ksi[0] = 1.178; p[0] = 10.928; q[0] = 3.139; latticeConstant = 4.085; r0[0] = latticeConstant / sqrt(2); // Ag
    A0[0] = 0.0854;
    A1[0] = 0.0;
    ksi[0] = 1.2243;
    p[0] = 10.939;
    q[0] = 2.2799;
    r0[0] = 2.5563;
    latticeConstant = 3.615;
    A0[1] = -0.0487;
    A1[1] = -0.7922;
    ksi[1] = 0.7356;
    p[1] = 8.1825;
    q[1] = 3.344;
    r0[1] = 2.4049;
    A0[2] = 0.1385;
    A1[2] = -0.3583;
    ksi[2] = 1.5247;
    p[2] = 7.6788;
    q[2] = 2.139;
    r0[2] = 2.3780;
    cohEnergy = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy());
    setV0();
    setC11_C12();
    setB();
    setC44();
    setEsol();
    setEnergyIn();
    setEnergyOn();
    cerr << "A0: " <<  A0[0] << " " << A0[1] << " " << A0[2] << "\n";
    cerr << "A1: " <<  A1[0] << " " << A1[1] << " " << A1[2] << "\n";
    cerr << "xi: " <<  ksi[0] << " " << ksi[1] << " " << ksi[2] << "\n";
    cerr << "p: " <<  p[0] << " " << p[1] << " " << p[2] << "\n";
    cerr << "q: " <<  q[0] << " " << q[1] << " " << q[2] << "\n";
    cerr << "r0: " <<  r0[0] << " " << r0[1] << " " << r0[2] << "\n";
    cerr << "a: " << latticeConstant << "\n";
    cerr << "coh = " << cohEnergy << "\n";
    cerr << "B = " << B << "\n";
    cerr << "C11 = " << C11 << "\n";
    cerr << "C12 = " << C12 << "\n";
    cerr << "C44  = " << C44 << "\n";
    cerr << "Esol = " << Esol << "\n";
    cerr << "energy in " << energyIn << "\n";
    cerr << "energy on " << energyOn << "\n";

}
