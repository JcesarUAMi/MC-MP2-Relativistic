#ifndef RELAINTEGRATION_H
#define RELAINTEGRATION_H 

#include "RelaMetropolis.h"
#include "WeightFunction.h"
#include "Spinors.h"
#include "Random.h"
#include "IntegralEvaluation.h"
#include <iostream>
#include <array>
#include <vector>
#include <iomanip>
#include <cstdlib>
using namespace std;

struct values {
  vector<double> factor;
};


class RelaIntegration {
  public:
    gaussKronrod gk;
    array<values, 2> occs;
    array<values, 2> virsA;
    array<values, 2> virsB;
    array<values, 2> vir13;
    array<values, 2> vir24;
    array<values, 2> vir23;
    array<values, 2> vir14;
    array<values, 2> occ13;
    array<values, 2> occ24;
		double lambda;
    double count;
		vector<double> taus;
		vector<double> wgTaus;
		vector<vector<double>> tauOccExpos;
		vector<vector<double>> tauVirExpos;
		array<expoValues, 21> RelaTauExpos;
    double FinalMP2Ener, FinalSCSEMP2;
		void gaussKronrodEval();
    void RelaMP2Energy(int, int, fstream&, fstream&, fstream&, fstream&, fstream&, fstream&);
    void MetropolisData(fstream&, fstream&, fstream&, fstream&, fstream&);
		void gaussKronrodRela();
		void tauIntegration();
		void oIJ(double&, double&, int, int, int);
		void oIJTau(array<values, 2>& occ13, double, int, int);
		void vABTau(array<values, 2>& vir13, double, int, int);
		void oIJK(double&, double&, int, int, int);
		void vABK(double&, double&, int, int, int);
		double wg = 0.5;
  private:
    RelaMetropolis MRC;
		relBasis rB;
		RandomWalk walker;
};




#endif
