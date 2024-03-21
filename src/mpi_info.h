#include <string>
#include <vector>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

class MPI_info {
 public:
  MPI_info();
  void print();

  static void barrier();
  static void broadcast_int(int*, size_t);
  static void broadcast_char(char*, size_t);
  static void broadcast_double(double*, size_t);
  static void broadcast_vector_double(std::vector<double>&);
  static void broadcast_string(std::string&);
  static void comm_rank(int*);
  static void comm_size(int*);

  int numtasks;
  int taskid;
  bool sys_master;
};

#endif 
