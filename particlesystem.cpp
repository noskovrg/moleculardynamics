#include "particlesystem.h"
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <cstring>
#include <cstdlib>
#include <string>
#include <fstream>

using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::string;
using std::ifstream;
using std::getline;

const double eps_b = 3e-2;
const double eps_ab = 3e-4; // -2
const double eps_a = 3e-2;

double ParticleSystem::euclideanNorm(Vector3d r)
{
    return sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
}

ParticleSystem::ParticleSystem(string & filename){
    threads = 1;
    isBulk = 1;
    xNumber = 3; yNumber = 3; zNumber = 3;
    printErrorFlag = printConfigurationFlag = 0;
    string thr_str = "Number of threads: ";
    string seed_str = "Seed: ";
    string nameAtom_str = "Name of atom: ";
    string latticeConstant_str = "latticeConstant: ";
    string cohEnergy_str = "cohesive energy: ";
    string b_str = "B: ";
    string c11_str = "C11: ";
    string c12_str = "C12: ";
    string c44_str = "C44: ";
    string namePurityAtom_str = "Name of impurity atom: ";
    string eSol_str = "e_sol: ";
    string energyIn_str = "e_in_dim: ";
    string energyOn_str = "e_on_dim: ";
    string x_str = "count of cells(x): ";
    string y_str = "count of cells(y): ";
    string z_str = "count of cells(z): ";
    string printError_str = "Print error: ";
    string printConfiguration_str = "Print configuration: ";
    string cohEnergyPurityAtom_str = "cohesive energy of purity atom: ";
    ifstream inFile(filename.c_str());
    while(inFile){
        string line;
        getline(inFile, line);
        if (line[0] == '#') continue;
        if (line.compare(0, seed_str.size(), seed_str) == 0){
            seed = atoi(line.c_str() + seed_str.size());
        }
        else if (line.compare(0, thr_str.size(), thr_str) == 0){
            threads = atoi(line.c_str() + thr_str.size());
        }
        else if (line.compare(0, nameAtom_str.size(), nameAtom_str) == 0){
            nameAtom = line.substr(nameAtom_str.size(), line.size() - nameAtom_str.size());
        }
        else if (line.compare(0, namePurityAtom_str.size(), namePurityAtom_str) == 0){
            namePurityAtom = line.substr(namePurityAtom_str.size(), line.size() - namePurityAtom_str.size());
        }
        else if (line.compare(0, latticeConstant_str.size(), latticeConstant_str) == 0){
            tableLatticeConstant = atof(line.c_str() + latticeConstant_str.size());
        }
        else if (line.compare(0, cohEnergy_str.size(), cohEnergy_str) == 0){
            tableCohesiveEnergy = atof(line.c_str() + cohEnergy_str.size());
        }
        else if (line.compare(0, b_str.size(), b_str) == 0){
            tableB = atof(line.c_str() + b_str.size());
        }
        else if (line.compare(0, c11_str.size(), c11_str) == 0){
            tableC11 = atof(line.c_str() + c11_str.size());
        }
        else if (line.compare(0, c12_str.size(), c12_str) == 0){
            tableC12 = atof(line.c_str() + c12_str.size());
        }
        else if (line.compare(0, c44_str.size(), c44_str) == 0){
            tableC44 = atof(line.c_str() + c44_str.size());
        }
        else if (line.compare(0, eSol_str.size(), eSol_str) == 0){
            tableEsol = atof(line.c_str() + eSol_str.size());
        }
        else if (line.compare(0, energyIn_str.size(), energyIn_str) == 0){
            tableEin = atof(line.c_str() + energyIn_str.size());
        }
        else if (line.compare(0, energyOn_str.size(), energyOn_str) == 0){
            tableEon = atof(line.c_str() + energyOn_str.size());
        }
        else if (line.compare(0, cohEnergyPurityAtom_str.size(), cohEnergyPurityAtom_str) == 0){
            cohEnergyPurityAtom = atof(line.c_str() + cohEnergyPurityAtom_str.size());
        }
        else if (line.compare(0, x_str.size(), x_str) == 0){
            xNumber = atoi(line.c_str() + x_str.size());
        }
        else if (line.compare(0, y_str.size(), y_str) == 0){
            yNumber = atoi(line.c_str() + y_str.size());
        }
        else if (line.compare(0, z_str.size(), z_str) == 0){
            zNumber = atoi(line.c_str() + z_str.size());
        }
        else if (line.compare(0, printError_str.size(), printError_str) == 0){
            printErrorFlag = atoi(line.c_str() + printError_str.size());
        }
        else if (line.compare(0, printConfiguration_str.size(), printConfiguration_str) == 0){
            printConfigurationFlag = atoi(line.c_str() + printConfiguration_str.size());
        }


        else if (line.size() >= 3 && (line.compare(3, 5, "B-B: ") == 0 ||
                                      line.compare(2, 5, "B-B: ") == 0 ||
                                      line.compare(4, 5, "B-B: ") == 0)){
            if (line.compare(0, 3, "A0 ") == 0){
                char * end;
                bb_a0[0] = strtod(line.c_str() + 8, &end);
                bb_a0[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 3, "A1 ") == 0){
                char * end;
                bb_a1[0] = strtod(line.c_str() + 8, &end);
                bb_a1[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 4, "ksi ") == 0){
                char * end;
                bb_ksi[0] = strtod(line.c_str() + 9, &end);
                bb_ksi[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 2, "p ") == 0){
                char * end;
                bb_p[0] = strtod(line.c_str() + 7, &end);
                bb_p[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 2, "q ") == 0){
                char * end;
                bb_q[0] = strtod(line.c_str() + 7, &end);
                bb_q[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 3, "r0 ") == 0){
                char * end;
                bb_r0[0] = strtod(line.c_str() + 8, &end);
                bb_r0[1] = strtod(end, NULL);
            }
        }
        else if (line.size() >= 3 && (line.compare(3, 5, "A-B: ") == 0 ||
                                      line.compare(2, 5, "A-B: ") == 0 ||
                                      line.compare(4, 5, "A-B: ") == 0)){
            if (line.compare(0, 3, "A0 ") == 0){
                char * end;
                ab_a0[0] = strtod(line.c_str() + 8, &end);
                ab_a0[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 3, "A1 ") == 0){
                char * end;
                ab_a1[0] = strtod(line.c_str() + 8, &end);
                ab_a1[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 4, "ksi ") == 0){
                char * end;
                ab_ksi[0] = strtod(line.c_str() + 9, &end);
                ab_ksi[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 2, "p ") == 0){
                char * end;
                ab_p[0] = strtod(line.c_str() + 7, &end);
                ab_p[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 2, "q ") == 0){
                char * end;
                ab_q[0] = strtod(line.c_str() + 7, &end);
                ab_q[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 3, "r0 ") == 0){
                char * end;
                ab_r0[0] = strtod(line.c_str() + 8, &end);
                ab_r0[1] = strtod(end, NULL);
            }
        }
        else if (line.size() >= 3 && (line.compare(3, 5, "A-A: ") == 0 ||
                                      line.compare(2, 5, "A-A: ") == 0 ||
                                      line.compare(4, 5, "A-A: ") == 0)){
            if (line.compare(0, 3, "A0 ") == 0){
                char * end;
                aa_a0[0] = strtod(line.c_str() + 8, &end);
                aa_a0[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 3, "A1 ") == 0){
                char * end;
                aa_a1[0] = strtod(line.c_str() + 8, &end);
                aa_a1[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 4, "ksi ") == 0){
                char * end;
                aa_ksi[0] = strtod(line.c_str() + 9, &end);
                aa_ksi[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 2, "p ") == 0){
                char * end;
                aa_p[0] = strtod(line.c_str() + 7, &end);
                aa_p[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 2, "q ") == 0){
                char * end;
                aa_q[0] = strtod(line.c_str() + 7, &end);
                aa_q[1] = strtod(end, NULL);
            }
            else if (line.compare(0, 3, "r0 ") == 0){
                char * end;
                aa_r0[0] = strtod(line.c_str() + 8, &end);
                aa_r0[1] = strtod(end, NULL);
            }
        }

    }
    OMP_ON = threads > 1;
    if (seed == -1){
        seed = time(NULL);
    }
    srand(seed);
    for (int itX = 0; itX < xNumber; ++itX){
        for (int itY = 0; itY < yNumber; ++itY){
            for (int itZ = 0; itZ < zNumber;  ++itZ){
                Particle tmpParticle0(itX, itY, itZ, 0);
                particles.push_back(tmpParticle0);
                Particle tmpParticle1(itX + 0.5, itY + 0.5, itZ, 0);
                Particle tmpParticle2(itX + 0.5, itY, itZ + 0.5, 0);
                Particle tmpParticle3(itX, itY + 0.5, itZ + 0.5, 0);
                particles.push_back(tmpParticle1);
                particles.push_back(tmpParticle2);
                particles.push_back(tmpParticle3);
            }
        }
    }
    if (printConfigurationFlag) printConfiguration();
}

void ParticleSystem::printConfiguration()
{
    cout << "Configuration of program system:\n______________________________________________________\n";
    cout << "number of particles: " << particles.size() << endl;
    cout << "Number of threads: " << threads << "\n";
    cout << "Seed: " << seed << "\n";
    cout << "Print error: " << printErrorFlag << "\n";
    cout << "Print configuration: " << printConfigurationFlag << "\n";
    cout << "count of cells(x): " << xNumber << "\n";
    cout << "count of cells(y): " << yNumber << "\n";
    cout << "count of cells(z): " << zNumber << "\n";
    cout << "Name of atom: " << nameAtom << "\n";
    cout << "latticeConstant: " << tableLatticeConstant << "\n";
    cout << "cohesive energy: " << tableCohesiveEnergy << "\n";
    cout << "B: " << tableB << "\n";
    cout << "C11: " << tableC11 << "\n";
    cout << "C12: " << tableC12 << "\n";
    cout << "c44: " << tableC44 << "\n";
    cout << "Name of purity atom: " << namePurityAtom << "\n";
    cout << "cohesive energy of purity atom: " << cohEnergyPurityAtom << "\n";
    cout << "e_sol: " << tableEsol << "\n";
    cout << "e_in_dim: " << tableEin << "\n";
    cout << "e_on_dim: " << tableEon << "\n";


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
    }
    for (int i = 0; i < threads; ++i){
        energy += tmpEnergy[i];
    }
    return energy / particles.size();
}



double ParticleSystem::cohesiveEnergy(int typeTransform, double alpha)
{
    cutoff = 1.7 * latticeConstant;
    double energy = 0.0;
    double repulsiveEnergy = 0.0;
    double bandEnergy = 0.0;
    for (int i = 0; i < particles.size(); ++i){
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

int comp(pair<double, int> a, pair<double, int> b){
    return a.first < b.first;
}

double ParticleSystem::NMA(double (ParticleSystem::*error)(), int _type, double alpha, double beta, double gamma)
{
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
            cout << "error: " << err << "\n";
            cout << "OPTIMIZED\n";
            print_params(_type);
            break;
        }
        for (int i = 0; i < DIM + 1; ++i){
            generate_params(_type, simplex[i]);
            cntSteps = 0;
        }
        while (true){
            cntSteps++;
            if (cntSteps > 1000){
                break;
            }
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
            // step 9
            if (printErrorFlag){
            	cout << cntIt++ << " " << globalError << "\n";
            }
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
}

void ParticleSystem::load_params(int type, vector<double> & params)
{
    A0[type] = params[0];
    A1[type] = params[1];
    ksi[type] = fabs(params[2]);
    p[type] = fabs(params[3]);
    q[type] = fabs(params[4]);
    r0[type] = fabs(params[5]);
    if (type == 0) {
        latticeConstant = r0[type] * sqrt(2);
    }
    setV0();
}

double ParticleSystem::error_b()
{
    cohEnergy = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy());
    setB();
    setC11_C12();
    setC44();
    return sqrt((pow((tableCohesiveEnergy - cohEnergy) / tableCohesiveEnergy, 2) +
                pow((tableLatticeConstant - latticeConstant) / tableLatticeConstant, 2)
                + pow((tableB - B) / tableB, 2)
                + pow((tableC11 - C11) / tableC11, 2)
                + pow((tableC12 - C12) / tableC12, 2)
                + pow((tableC44 - C44) / tableC44, 2)) / 6.0);
}

double ParticleSystem::error_ab()
{
    setEsol();
    double err = sqrt(pow((tableEsol - Esol) / tableEsol, 2));
    return err;
}

double ParticleSystem::error_a()
{
    setEnergyIn();
    setEnergyOn();
    double err = sqrt((pow((tableEin - energyIn) / tableEin, 2) +
             pow((tableEon - energyOn) / tableEon, 2)) / 2.0);
    return err;
}



void ParticleSystem::setEnergyIn()
{
    int first = -1;
    int second = -1;
    for (int i = 0; i < particles.size(); ++i){
        Vector3d pos = particles[i].getPositionLC();
        if (pos.z == zNumber - 0.5 && pos.x == 0.5 && pos.y == 0){
            first = i;
        }
        if (pos.z == zNumber - 0.5 && pos.x == 0 && pos.y == 0.5){
            second = i;
        }
    }
    double res = 0;
    isBulk = 0;
    double energy_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    particles[first].setType(1);
    double energy_adatom_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    particles[second].setType(1);
    double energy_dimer_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    res = (energy_dimer_surf - energy_surf) - 2 * (energy_adatom_surf - energy_surf);
    particles[first].setType(0);
    particles[second].setType(0);
    isBulk = 1;
    energyIn = res;
}

void ParticleSystem::setEnergyOn()
{
    double res = 0;
    Particle tmp0(0.0, 0.0, zNumber, 1);
    Particle tmp1(0.5, 0.5, zNumber, 1);
    isBulk = 0;
    double energy_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    particles.push_back(tmp0);
    double energy_adatom_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    particles.push_back(tmp1);
    double energy_dimer_surf = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
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
    cout << "A0 = " << A0[type] << "\nA1 = " << A1[type] << "\nksi = " << ksi[type] <<
            "\np = " << p[type] << "\nq = " << q[type] << "\nr0 = " << r0[type] << "\n";
   // cout << A0[type] << ", " << A1[type] << ", " << p[type] << ", " << q[type] << ", " <<
     //       ksi[type] << ", " << r0[type] << "\n";
}

void ParticleSystem::setV0()
{
    V0 = latticeConstant * latticeConstant * latticeConstant * xNumber * yNumber * zNumber;
}

void ParticleSystem::printV0()
{
    setV0();
    cout << "V0 = " << V0 << endl;
}

void ParticleSystem::setC11_C12()
{
    C11 = 0;
    C12 = 0;
    {
    double energy_i1 = (OMP_ON ? cohesiveEnergy_omp(11, 0.0): cohesiveEnergy()) * particles.size();
    double energy_i0 = (OMP_ON ? cohesiveEnergy_omp(11, -0.01): cohesiveEnergy(11, -0.01)) * particles.size();
    double energy_i2 = (OMP_ON ? cohesiveEnergy_omp(11, 0.01): cohesiveEnergy(11, 0.01)) * particles.size();
    C11 += 1. / V0 * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001;
    C12 += 1. / V0 * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001 ;
    }
    {
    double energy_i1 = (OMP_ON ? cohesiveEnergy_omp(12, 0.0): cohesiveEnergy(12, 0.0)) * particles.size();
    double energy_i0 = (OMP_ON ? cohesiveEnergy_omp(12, -0.01): cohesiveEnergy(12, -0.01)) * particles.size();
    double energy_i2 = (OMP_ON ? cohesiveEnergy_omp(12, 0.01): cohesiveEnergy(12, 0.01)) * particles.size();
    C11 += 1. / V0 * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001;
    C12 -= 1. / V0 * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001;
    }
    C11 /= 2;
    C12 /= 2;
}

void ParticleSystem::printC11_C12()
{
    setC11_C12();
    cout << "C11 = " << C11 << endl;
    cout << "C12 = " << C12 << endl;
}


void ParticleSystem::setC44(){
    double energy_i1 = (OMP_ON ? cohesiveEnergy_omp(44, 0.0): cohesiveEnergy(44, 0.0)) * particles.size();
    double energy_i0 = (OMP_ON ? cohesiveEnergy_omp(44, -0.01): cohesiveEnergy(44, -0.01)) * particles.size();
    double energy_i2 = (OMP_ON ? cohesiveEnergy_omp(44, 0.01): cohesiveEnergy(44, 0.01)) * particles.size();
    C44 = 1.0 / (2.0 * V0) * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001;
}

void ParticleSystem::printC44(){
    setC44();
    cout << "C44 = " << C44 << endl;
}

void ParticleSystem::setB(){
    double energy_i1 = (OMP_ON ? cohesiveEnergy_omp(0, 0.0): cohesiveEnergy(0, 0.0)) * particles.size();
    double energy_i0 = (OMP_ON ? cohesiveEnergy_omp(0, -0.01): cohesiveEnergy(0, -0.01)) * particles.size();
    double energy_i2 = (OMP_ON ? cohesiveEnergy_omp(0, 0.01): cohesiveEnergy(0, 0.01)) * particles.size();
    B = 2.0 / (9.0 * V0) * (energy_i0 - 2 * energy_i1 + energy_i2) / 0.0001;
}


void ParticleSystem::setEsol()
{
    double cohEnergy_b = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy());
    double energy_b = cohEnergy_b * particles.size();
    particles[0].setType(1);
    double energy_ab = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy()) * particles.size();
    Esol = energy_ab - energy_b - cohEnergyPurityAtom + cohEnergy_b;
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
    return true;
}

int ParticleSystem::type(int i, int j)
{
    return particles[i].getType() + particles[j].getType();
}

void ParticleSystem::fittingB()
{
    if (printConfigurationFlag) cout << "\n";
    cout << "Values by fitted params for " + nameAtom + "-" + nameAtom + " type interconection of elements\n";
    NMA(&ParticleSystem::error_b, 0);
    cout << "lattice constant = " << latticeConstant << "\n";
    cohEnergy = (OMP_ON ? cohesiveEnergy_omp(): cohesiveEnergy());
    cout <<  "Cohesive energy: " << cohEnergy << endl;
    printB();
    printC11_C12();
    printC44();
}

void ParticleSystem::fittingAB()
{
    cout << "\n\nValues by fitted params for " + namePurityAtom + "-" + nameAtom + " type interconection of elements\n";
    NMA(&ParticleSystem::error_ab, 1);
    cout << "Esol = " << Esol << "\n";
}

void ParticleSystem::fittingA()
{
    cout << "\n\nValues by fitted params for " + namePurityAtom + "-" + namePurityAtom + " type interconection of elements\n";
    NMA(&ParticleSystem::error_a, 2);
    cout << "Energy dim in = " << energyIn << "\n";
    cout << "Energy dim out = " << energyOn << "\n";
}

void ParticleSystem::printEsol()
{
    setEsol();
    cout << "Esol = " << Esol << endl;
}

void ParticleSystem::printB()
{
    setB();
    cout << "B = " << B << endl;
}



void ParticleSystem::generate_params(int _type, vector<double> & params)
{
    if (_type == 0){
        params[0] = random() / (double)RAND_MAX * (bb_a0[1] - bb_a0[0]) + bb_a0[0];
        params[1] = random() / (double)RAND_MAX * (bb_a1[1] - bb_a1[0]) + bb_a1[0];
        params[2] = random() / (double)RAND_MAX * (bb_ksi[1] - bb_ksi[0]) + bb_ksi[0];
        params[3] = random() / (double)RAND_MAX * (bb_p[1] - bb_p[0]) + bb_p[0];
        params[4] = random() / (double)RAND_MAX * (bb_q[1] - bb_q[0]) + bb_q[0];
        params[5] = random() / (double)RAND_MAX * (bb_r0[1] - bb_r0[0]) + bb_r0[0];
    }
    else if (_type == 1){
        params[0] = random() / (double)RAND_MAX * (ab_a0[1] - ab_a0[0]) + ab_a0[0];
        params[1] = random() / (double)RAND_MAX * (ab_a1[1] - ab_a1[0]) + ab_a1[0];
        params[2] = random() / (double)RAND_MAX * (ab_ksi[1] - ab_ksi[0]) + ab_ksi[0];
        params[3] = random() / (double)RAND_MAX * (ab_p[1] - ab_p[0]) + ab_p[0];
        params[4] = random() / (double)RAND_MAX * (ab_q[1] - ab_q[0]) + ab_q[0];
        params[5] = random() / (double)RAND_MAX * (ab_r0[1] - ab_r0[0]) + ab_r0[0];
    }
    else{
        params[0] = random() / (double)RAND_MAX * (aa_a0[1] - aa_a0[0]) + aa_a0[0];
        params[1] = random() / (double)RAND_MAX * (aa_a1[1] - aa_a1[0]) + aa_a1[0];
        params[2] = random() / (double)RAND_MAX * (aa_ksi[1] - aa_ksi[0]) + aa_ksi[0];
        params[3] = random() / (double)RAND_MAX * (aa_p[1] - aa_p[0]) + aa_p[0];
        params[4] = random() / (double)RAND_MAX * (aa_q[1] - aa_q[0]) + aa_q[0];
        params[5] = random() / (double)RAND_MAX * (aa_r0[1] - aa_r0[0]) + aa_r0[0];
    }
}

