#ifndef MPI_H
#define MPI_H

#include <string>
#include <vector>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

class MPI_info {
 public:
  MPI_info();
  void print();

  int numtasks;
  int taskid;
  bool sys_master;
};

#endif 
