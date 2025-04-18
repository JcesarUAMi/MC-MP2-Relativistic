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

#include "mpi.h"

int main (int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  MPI_info mpi_info;

  if (argc != 2) {
    if (mpi_info.sys_master) {
      printf("Usage: MC-Rela <input>\n");
    }
  } else {
    if (mpi_info.sys_master) {
      printf("Program developed by Julio Cruz\n");
    }
  }
  
  mpi_info.print();

  IOPs iops;

  iops.read(mpi_info, argv[1]);
  iops.print(mpi_info, argv[1]);

  Integration ie;
  Molecule molec;

  molec.readCoordinates(mpi_info, iops.sopns[KEYS::GEOM]);

  if (iops.iopns[KEYS::JOBTYPE] == JOBTYPE::RELATIVISTIC) {
  
    cout << "Relativistic Type Job" << endl;


  } else if (iops.iopns[KEYS::JOBTYPE] == JOBTYPE::NON_RELATIVISTIC) {
      
    cout << "NON Relativistic Type Jpb" << endl;
  }



  MPI_Finalize();

  return 0;
}

