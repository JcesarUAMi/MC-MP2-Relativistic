#ifndef MPI_H
#define MPI_H

class MPI_info {
 public:
  MPI_info();
  void print();

  int numtasks;
  int taskid;
  bool sys_master;
};

#endif 
