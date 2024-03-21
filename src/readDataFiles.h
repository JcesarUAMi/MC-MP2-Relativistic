#ifndef READ_H
#define READ_H

#include "Basis.h"
#include "Geometry.h"
#include "Coefficients.h"
#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

struct basisCoeff {
  vector<vector<double>> bC;
};

struct finalBasisData { 
  vector<double> basisExponents;
  vector<double> basisCoefficients;
  vector<basisCoeff> lastBasisCoeff;
  vector<int> shellSize;
  vector<int> functionType;
  vector<int> angularMoment;
  array<double, 3> nucleusCoord;
  int atomType;
}; 

struct densMatrix {
  vector<vector<double>> dE;
};

class readDataInput {
  public:
    int numberOfMolOrb;
    int numberOfOccOrb;
    int numberOfVirOrb;
    int numberOfPrimFunc;
    vector<double> moi;
    vector<vector<double>> moisCoeffs;
    vector<densMatrix> dEi;
    vector<finalBasisData> lastBasis;
    double generateMoiSpherical(double, double, double);
    double SphericalEvaluation(double, double, double, double, int, int);
    void densityMatrix();
    void coefficientsFinal(fstream&, fstream&, fstream&);
    void nucleusBasis();
    void getNumOccOrb();
    vector<double> Occ;
    vector<double> Vir;
  private:
        double sc[7] = {sqrt(3),
                sqrt(5/2),
                sqrt(15),
                sqrt(3/2),
                sqrt(35),
                sqrt(35/2),
                sqrt(5)};
    Basis basis;
    Molecule molec;
    Coefficients coeff;

};

#endif


