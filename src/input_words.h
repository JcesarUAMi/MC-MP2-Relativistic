#ifndef INPUT_H_
#define INPUT_H_

#include <array>
#include <string>
#include <vector>
#include <map>

namespace KEYS {
  enum KEY_TYPE {
    STRING,
    INT,
    DOUBLE,
    BOOL,
    OTHER
  };

  enum KEYS {
#define FORMAT(X, Y) X,
#include "keys.h"
#undef FORMAT
  };
  const std::map<KEYS, KEY_TYPE> KEY_TYPE_TABLE = {
#define FORMAT(X, Y) {X, Y},
#include "keys.h"
#undef FORMAT
  };
  const std::vector<std::string> key_strings = {
#define FORMAT(X, Y) #X,
#include "keys.h"
#undef FORMAT
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
