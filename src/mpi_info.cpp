#ifndef MPI_INFO_H
#define MPI_INFO_H

#include "mpi.h"

#include <cstdio>
#include "mpi_info.h"

MPI_info::MPI_info() {
  /*
   * MPI_info constructor
   */
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

  sys_master = 0 == taskid;
  

}

void MPI_info::print() {
  /*
   * Print number of mpi tasks if sys_master
   */
  if (sys_master) {
    printf("MPI IS RUNNING WITH %i CPUS\n", numtasks);
  }
}

#endif

