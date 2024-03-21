#ifndef WEIGHT_H
#define WEIGHT_H

#include "Geometry.h"
#include "Random.h"
#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <cstdlib>
using namespace std;

struct weightsValues {
  array<double, 2> coeffWeight;
  array<double, 2> expWeight;
	array<double, 2> normWeight;
  int atomType;
};

struct finalWeightsData {
  array<double, 2> coeffWeight;
  array<double, 2> expWeight;
  array<double, 3> nucleusCoord;
	array<double, 2> normWeight;
  int atomType;
};


class Weights {
  private:
    Molecule molec;
    RandomWalk walker;
  public:
		double normalize(double);
    void Data(fstream&);
    double constant();
		const double pi25 = pow(M_PI,2.5);
    double weightIntegral(double, double, double);
    double weightFunction(double, double, double, double, double, double, double);
    vector<weightsValues> wV;
    vector<finalWeightsData> wD;
};

#endif



