#ifndef GRID_H
#define GRID_H

#include "Basis.h"
#include "Geometry.h"
#include "Coefficients.h"
#include "readDataFiles.h"
#include "IntegralEvaluation.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>
using namespace std;

class Grid {
  private:
    readDataInput readData;
  public:
    void detMinMax ();
    void Data(fstream&, fstream&, fstream&);
    string size = "Medium";
    void defNumSteps ();
    array<double, 3> sizeStep;
    array<vector<double>, 3> pointsGrid;
    array<double, 6> minmax;
    int totPoints;
    double factorTot;
    void evaluateDen();
    vector<double> den;

};

#endif


