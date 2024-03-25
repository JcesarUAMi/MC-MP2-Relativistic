#ifndef INPUT_H_
#define INPUT_H_

#include <array>
#include <string>
#include <vector>
#include <map>

#include "mpi_info.h"

namespace KEYS {
  enum KeyVal {
    JOBNAME = 0,
    JOBTYPE,
    SPHERICAL,
    FREEZE_CORE,
    MC_TRIAL,
    ELECTRON_PAIRS,
    MC_DELX,  // 0-4
    GEOM,
    BASIS,
    MC_BASIS,
    NBLOCK,
    MOVECS,  // 5-9
    DEBUG,
  };
}

namespace JOBTYPE {
  enum JOBTYPE {
    RELATIVISTIC,
    NON_RELATIVISTIC,
  };
  const std::vector<std::string> jobtype_strings = {
    "RELATIVISTIC",
    "NON_RELATIVISTIC",
  };
}

class IOPs {
 public:
  IOPs();

  void print(const MPI_info&, const std::string&);
  void read(const MPI_info&, const std::string&);

  std::array<int, 100> iopns;
  std::array<double, 100> dopns;
  std::array<bool, 100> bopns;
  std::array<std::string, 100> sopns;
};

#endif
