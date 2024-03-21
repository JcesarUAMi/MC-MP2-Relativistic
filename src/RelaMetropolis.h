#ifndef RELAMETROPOLIS_H
#define RELAMETROPOLIS_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <array>
#include <cstdlib>
#include "Random.h"
#include "Spinors.h"
using namespace std;

struct RelaSteps {
	array<double, 3> e1;
  array<double, 3> e2;
  double w12;
  double r12;
};

struct RelaWalkPairs {
	vector<RelaSteps> rSteps;
};

class RelaMetropolis {
  private:
    RandomWalk walker;
		relBasis rB;
  public:
    void motion(int);
    void startWalking(double&, double&, double&, double&, double&, double&, double&, double&);
    unsigned int fSteps;
    unsigned int sSteps;
    int numberWalkers;
    int Nbsize;
    unsigned long long int reSteps;
    unsigned long long int nSteps;
    vector<RelaWalkPairs> wp;
    double Eg;
    double moveLength;
    void moveScheme();
    void burn(int);
    double rescaleMoveLengthBurn(int, double);
    double rescaleMoveLengthStep(int, double);
		void Data(fstream&, fstream&);
};

#endif
