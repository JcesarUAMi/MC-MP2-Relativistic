#ifndef SPINORS_H
#define SPINORS_H


#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <algorithm>
#include "WeightFunction.h"
using namespace std;

struct geomInfo {
  string atomType;
  int atomicNum;
  double charge;
  vector<int> rpsb;
  array<double, 3> nucCoor;
};

class relGeometry {
  public:
    int atmNum;
    int nE;
		int nDifAtoms;
		vector<int> difAtoms;
    vector<geomInfo> geoRel;
    void geomLecture (fstream&);
};

struct cgto {
  int n;
  int k;
  vector<double> ak;
  vector<double> RL1;
  vector<double> RS1;
  vector<double> RL2;
  vector<double> RS2;
};

struct basisData {
  int l;
  vector<cgto> cg;
};

struct cgcValues {
  vector<vector<double>> cgc;
};

struct coeffValues {
  vector<double> moRe;
  vector<double> moIm;
};

struct spiMos {
  vector<coeffValues> moS;
  vector<vector<double>> ov;
};

struct valaoStruct {
  array<vector<double>, 2> val;
};

struct valmoValues {
  array<vector<double>, 2> valuesMos;
};

struct Spinors {
  array<valmoValues, 4> Occvalmo;
  array<valmoValues, 4> Virvalmo;
};

struct basisInfo {
  int maxBid;
  int maxl;
  int rdim;
  int occOrb, virOrb;
  vector<basisData> bD;
  array<cgcValues, 4> cG;
  array<vector<double>, 2> sph;
  vector<double> sqV;
  array<spiMos, 2> mo;
  vector<double> orbEne;
  array<valaoStruct, 4> valao;
  array<Spinors, 4> moSpi;
};

class relBasis {
  public:
    basisInfo bI;
    relGeometry rG;
		Weights wg;
    void readBasisRel (fstream&, fstream&, fstream&, fstream&, fstream&);
    void readOrbitalCoeff(fstream&);
    void readOrbitalEne(fstream&);
		void generateSpinors(int);
    void mofunc(double, double, double, int);
    void make_sph(int, double, double, double, double, double);
    double complexReal(double, double, double, double);
    double complexIm(double, double, double, double);
    void source();
		void weightFunctionData(fstream&);
		double constant();
		double weightFunction(double, double, double, double, double, double, double);
    int nWalkers;
    int nSteps;
    int Nbsize;
};

#endif
