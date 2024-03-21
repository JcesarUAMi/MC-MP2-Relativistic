#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
using namespace std;
#include "Coefficients.h"
#include "Geometry.h"
#include "Basis.h"
#include "readDataFiles.h"
#include "WeightFunction.h"
#include "Metropolis.h"
#include "RelaMetropolis.h"
#include "Random.h"
#include "IntegralEvaluation.h"
#include "RelaIntEval.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif



int main (int argc, char *argv[]) {

// For non-relativistica calculations
  fstream xyz;
  fstream movecs;
  fstream basis;
  Integration ie;
//  Grid gd;
//  gd.Data(basis, xyz, movecs);
  ie.Data(basis, xyz, movecs);
  ie.MP2Energy(xyz);

  int mynode, totalnodes, initialized, finalized;
// For Relativistic calculations
  fstream info;
  fstream coeffs;
  fstream ene;
  fstream MC;
  RelaIntegration MCR;
  

/*
  auto begin = chrono::high_resolution_clock::now();

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  MCR.RelaMP2Energy(mynode, totalnodes, info, basis, xyz, coeffs, ene, MC);

  MPI_Finalize();

  auto end = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

  cout << "Time measured: " <<  setprecision(6) << elapsed.count() * 1e-9 << endl;
*/

  return 0;
}

