#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


#include "input_words.h"
#include "mpi_info.h"
#include "mpi.h"

IOPs::IOPs() {
  /*
   * IOPs constructor
   * initializes variables to mostly rational values
   *
   * TODO
   *  -read values directly instead of setting them
   */
  bopns.fill(false);
  dopns.fill(0.0);
  iopns.fill(0);
  bopns[KEYS::SPHERICAL] = true;
  bopns[KEYS::FREEZE_CORE] = true;
  dopns[KEYS::MC_DELX] = 0.1;
  iopns[KEYS::ELECTRON_PAIRS] = 16;
  iopns[KEYS::MC_TRIAL] = 1024;
  iopns[KEYS::NBLOCK] = 1;
  iopns[KEYS::JOBTYPE] = JOBTYPE::NON_RELATIVISTIC;

  sopns[KEYS::GEOM] = "geom.xyz";
  sopns[KEYS::BASIS] = "basis.dat";
  sopns[KEYS::MC_BASIS] = "mc_basis.dat";
  sopns[KEYS::MOVECS] = "nwchem.movecs";
  iopns[KEYS::MOVECS] = 0;  // default is binary
}

void IOPs::read(const MPI_info& mpi_info, const std::string& file) {
  /*
   * reads and stores options mcin file provided as command line argument
   *
   * Arguments
   *
   *  -needs input validation
   *  -make key setter by type
   *  -write functions to convert strings to enums
   *  -clean up keys
   */
  KEYS::KeyVal keyval;

  bool keySet;
  std::string str;
  std::string key;

  const std::vector<std::string> key_vals = {
      "JOBNAME", "JOBTYPE", "SPHERICAL", "FREEZE_CORE", "MC_TRIAL",  // 0-4 
      "ELECTRON_PAIRS", "MC_DELX", "GEOM", "BASIS", "MC_BASIS",      // 5-9
      "NBLOCK", "MOVECS", "DEBUG"};                                  // 10-12

  if (mpi_info.sys_master) {
    std::ifstream input(file.c_str());
    if (!input.is_open()) {
      std::cerr << "FILE \"" << file << "\" DOES NOT EXIST" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    keySet = false;
    while (getline(input, str)) {
      std::istringstream ss(str);
      while (ss >> key) {
        if (keySet == false) {  // if key not set, determine key value from key_vals arrays
          std::transform(key.begin(), key.end(), key.begin(), ::toupper);
          auto it = std::find(key_vals.begin(), key_vals.end(), key);
          if (it != key_vals.end()) {
            keyval = static_cast<KEYS::KeyVal>(std::distance(key_vals.begin(), it));
          } else {
            std::cerr << "Key " << key << " not reconginzed" << std::endl;
            exit(0);
          }
          keySet = true;
        } else {
            switch (keyval) {
              case KEYS::JOBTYPE:
              std::transform(key.begin(), key.end(), key.begin(), ::toupper);
              {
                auto it = std::find(key_vals.begin(), key_vals.end(), key);
                if (it != key_vals.end()) {
                  iopns[keyval] = std::distance(key_vals.begin(), it);
                } else {
                  std::cerr << "Job type " << key << " not reconginzed" << std::endl;
                  exit(0);
                }
              }
              keySet = false;
              break;
            case KEYS::JOBNAME:
              sopns[keyval] = key;
              keySet = false;
              break;
            case KEYS::SPHERICAL:
              bopns[keyval] = stoi(key, nullptr);
              keySet = false;
              break;
            case KEYS::MC_TRIAL:
              iopns[keyval] = stoi(key, nullptr);
              keySet = false;
              break;
            case KEYS::ELECTRON_PAIRS:
              iopns[keyval] = stoi(key, nullptr);
              keySet = false;
              break;
            case KEYS::MC_DELX:
              dopns[keyval] = stod(key, nullptr);
              keySet = false;
              break;
            case KEYS::GEOM:
              sopns[keyval] = key;
              keySet = false;
              break;
            case KEYS::BASIS:
              sopns[keyval] = key;
              keySet = false;
              break;
            case KEYS::MC_BASIS:
              sopns[keyval] = key;
              keySet = false;
              break;
            case KEYS::NBLOCK:
              iopns[keyval] = stoi(key, nullptr);
              if (iopns[keyval] > 20) {
                std::cerr << "NBlock must be less than 20" << std::endl;
                exit(EXIT_FAILURE);
              } else if (iopns[keyval] <= 0) {
                iopns[keyval] = 1;
              }
              keySet = false;
              break;
            case KEYS::MOVECS:
              if (key == "ASCII") {
                iopns[keyval] = 1;
              } else if (key == "END") {
                keySet = false;
              } else {
                sopns[keyval] = key;
              }
              break;
            default:
              std::cerr << "KEY \"" << key << "\" NOT RECONGNIZED" << std::endl;
              exit(EXIT_FAILURE);
          }
        }
      }
    }
    input.close();
  }



  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(iopns.data(), iopns.size(), MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(dopns.data(), dopns.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(bopns.data(), bopns.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

}

void IOPs::print(const MPI_info& mpi_info, const std::string& file) {

  if (mpi_info.sys_master) {
    std::cout << std::endl;
    std::cout << "Input read from " << file << std::endl;
    std::cout << "JOBNAME: " << sopns[KEYS::JOBNAME] << std::endl;
    std::cout << " JOBTYPE: " << JOBTYPE::jobtype_strings[iopns[KEYS::JOBTYPE]] << std::endl;
    std::cout << " MC_TRIAL: " << iopns[KEYS::MC_TRIAL] << std::endl;
    std::cout << " Electron Pairs: " << iopns[KEYS::ELECTRON_PAIRS] << std::endl;
    std::cout << " FREEZE_CORE: " << bopns[KEYS::FREEZE_CORE] << std::endl;
    std::cout << " SPHERICAL: " << bopns[KEYS::SPHERICAL] << std::endl;

    std::cout << std::endl;
  }
}

