#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <array>
#include <cstdlib>
#include "Random.h"
#include "WeightFunction.h"
using namespace std;

struct Steps {
	array<double, 3> e1;
  array<double, 3> e2;
  double w12;
  double r12;
};

struct walkPairs {
	vector<Steps> rSteps;
  double moveLength;
};

class Metropolis {
  private:
    RandomWalk walker;
    Weights wg;
  public:
    void startWalking(double&, double&, double&, double&, double&, double&, double&, double&);
    void motion();
    unsigned int fSteps;
    unsigned int sSteps;
		int numberWalkers = 16;
    unsigned long long int reSteps;
    unsigned long long int nSteps;
    vector<walkPairs> wp;
    void Data(fstream&);
    double Eg;
    double moveLength;
    void moveScheme();
    void burn();
    double rescaleMoveLengthBurn(int, double);
    double rescaleMoveLengthStep(int, double);
};

#endif
