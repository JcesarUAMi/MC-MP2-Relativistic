#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <string>
#include <array>
#include <vector>
using namespace std;

#include "mpi_info.h"

struct Geometry {
  array<double, 3> nCoor;
  int nAtomic;
};

class Molecule {
  private:
    string atomType;
  public:
    int nNuc;
    int nDifAtoms;
    int atomicNumber;
    void getNumberOfDifAtoms();
    void readCoordinates(MPI_info&, string&);
    void getAtomicNumber (string);
    vector<Geometry> geom;
    vector<int> numDifAtoms;

};


#endif
