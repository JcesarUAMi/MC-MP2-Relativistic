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
#include "mpi_info.h"
#include "input_words.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main (int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  MPI_info mpi_info;

  if (argc != 2) {
    if (mpi_info.sys_master) {
      printf("Usage: mcmpN.x <input>\n");
    }
    exit(EXIT_FAILURE);
  } else {
    if (mpi_info.sys_master) {
      printf("MC-GFn program developed by the Hirata lab\n");
    }
  }
  mpi_info.print();

  IOPs iops;
  iops.read(mpi_info, argv[1]);
  iops.print(mpi_info, argv[1]);

  Integration ie;
  

  return 0;
}

