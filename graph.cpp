#include <iostream>
#include <random>
#include <iomanip>
#include <chrono>
#include <cstdlib>
#include <fstream>
 
using namespace std;
 
int main(int argc,char *argv[])
{
    double a00 = 4.03469;
    double A00 = 0.104586; 
    double A10 = 0.010153;
    double p0 = 10.9801;
    double q0 = 3.1301;
    double ksi0 = 1.17288; 
    double r00 = 2.85296;
 
   /* double A01 = strtod(start, &end);
    double A11 = strtod(start, &end);
    double p1 = strtod(start, &end);
    double q1 = strtod(start, &end);
    double ksi1 = strtod(start, &end);
    double r01 = strtod(start, &end);
  
    double A02 = strtod(start, &end);
    double A12 = strtod(start, &end);
    double p2 = strtod(start, &end);
    double q2 = strtod(start, &end);
    double ksi2 = strtod(start, &end);
    double r02 = strtod(start, &end);
   
    double a01 = r01*sqrt(2), a02 = r02*sqrt(2);*/
 
    ofstream BB,AB,AA;
    BB.open("GRAPHICBB");
   // AB.open("GRAPHICAB");
   // AA.open("GRAPHICAA");
    for(double rij = 1.8; rij < 1.5*a00; rij+=0.1){
        double Efull = 0;
        Efull = (2.0 * (A10*(rij - r00) + A00)*exp(-abs(p0)*(rij/r00 - 1))
                         - sqrt(2.0 * ksi0*ksi0*exp(-2*abs(q0)*(rij/r00 - 1.))));
        BB << rij << "\t" << Efull << endl;

	cerr << 2.0 * (A10*(rij - r00) + A00)*exp(-abs(p0)*(rij/r00 - 1)) << " " << A10*(rij - r00) + A00 << " "
	 << exp(-abs(p0)*(rij/r00 - 1)) << "\n";



	// << - sqrt(2.0 * ksi0*ksi0*exp(-2*abs(q0)*(rij/r00 - 1.)))  << "\n";
    }
  /*  for(double rij = 1.8; rij < 1.5*a01; rij+=0.1){
        double Efull = 0;
        Efull = (A11*(rij - r01) + A01)*exp(-abs(p1)*(rij/r01 - 1)) + (A11*(rij - r01) + A01)*exp(-abs(p1)*(rij/r01 - 1))
                         - sqrt(ksi1*ksi1*exp(-2*abs(q1)*(rij/r01 - 1.)) + ksi1*ksi1*exp(-2*abs(q1)*(rij/r01 - 1.)));
        AB << rij << "\t" << Efull << endl;
    }
    for(double rij = 1.8; rij < 4.5*a02; rij+=0.1){
        double Efull = 0;
        Efull = (A12*(rij - r02) + A02)*exp(-abs(p2)*(rij/r02 - 1)) + (A12*(rij - r02) + A02)*exp(-abs(p2)*(rij/r02 - 1))
                         - sqrt(ksi2*ksi2*exp(-2*abs(q2)*(rij/r02 - 1.)) + ksi2*ksi2*exp(-2*abs(q2)*(rij/r02 - 1.)));
        AA << rij << "\t" << Efull << endl;
    }*/
 
    return 0;
}
