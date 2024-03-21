#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include "makeGrid.h"
using namespace std;

void Grid::Data(fstream& basisFile, fstream& xyzFile, fstream& movecsFile) {

  readData.coefficientsFinal(basisFile, xyzFile, movecsFile);
  readData.nucleusBasis();
  readData.getNumOccOrb();
//  readData.densityMatrix();
  evaluateDen();
}

void Grid::detMinMax () {

  double minx, miny, minz;
  double maxx, maxy, maxz;
  double Ax, Ay, Az;

  minx = readData.lastBasis[0].nucleusCoord[0];
  miny = readData.lastBasis[0].nucleusCoord[1];
  minz = readData.lastBasis[0].nucleusCoord[2];
  
  maxx = readData.lastBasis[0].nucleusCoord[0];
  maxy = readData.lastBasis[0].nucleusCoord[1];
  maxz = readData.lastBasis[0].nucleusCoord[2];

  for (int i=1; i<readData.lastBasis.size(); i++) {
    Ax = readData.lastBasis[i].nucleusCoord[0];
    Ay = readData.lastBasis[i].nucleusCoord[1];
    Az = readData.lastBasis[i].nucleusCoord[2];

    if ( Ax <= minx)
      minx = Ax;
    if ( Ay <= miny)
      miny = Ay;
    if ( Az <= minz)
      minz = Az;

    if ( Ax >= maxx)
      maxx = Ax;
    if ( Ay >= maxy)
      maxy = Ay;
    if ( Az >= maxz)
      maxz = Az;
  }

  minmax[0] = minx - 4.0;
  minmax[1] = miny - 4.0;
  minmax[2] = minz - 4.0;
  minmax[3] = maxx + 4.0;
  minmax[4] = maxy + 4.0;
  minmax[5] = maxz + 4.0;

}

void Grid::defNumSteps () {
  double distx, disty, distz;
  double pointsx, pointsy, pointsz;
  int numSteps;
  distx = abs(minmax[3] - minmax[0]);
  disty = abs(minmax[4] - minmax[1]);
  distz = abs(minmax[5] - minmax[2]);

  if (size == "Large")
    numSteps = 200;
  else if (size == "Medium")
    numSteps = 100;
  else if (size == "Short")
    numSteps = 50;
  else
    numSteps = 50;

  pointsx = distx * numSteps;
  pointsy = disty * numSteps;
  pointsz = distz * numSteps;

  sizeStep[0] = distx/numSteps;
  sizeStep[1] = disty/numSteps;
  sizeStep[2] = distz/numSteps;

  factorTot = sizeStep[0] * sizeStep[1] * sizeStep[2];

  pointsGrid[0].resize(numSteps+1);
  pointsGrid[1].resize(numSteps+1);
  pointsGrid[2].resize(numSteps+1);

  totPoints = pow(numSteps+1,3);
  
  den.resize(totPoints);

  pointsGrid[0][0] = minmax[0];
  pointsGrid[1][0] = minmax[1];
  pointsGrid[2][0] = minmax[2];

  for (int i=1; i<pointsGrid[0].size(); i++)
    pointsGrid[0][i] = minmax[0] + sizeStep[0] * i;
  
  for (int i=1; i<pointsGrid[1].size(); i++)
    pointsGrid[1][i] = minmax[1] + sizeStep[1] * i;
  
  for (int i=1; i<pointsGrid[2].size(); i++)
    pointsGrid[2][i] = minmax[2] + sizeStep[2] * i;
/*
  for (int i=0; i<pointsGrid[0].size(); i++)
    for (int j=0; j<pointsGrid[1].size(); j++)
      for (int k=0; k<pointsGrid[2].size(); k++)
        cout << i << " " << j << " " << k << " " << pointsGrid[0][i] << " " << pointsGrid[1][j] << " " << pointsGrid[2][k] << endl;
*/
}

void Grid::evaluateDen() {
  
  int mu;
  double x, y, z;
  double totalDen=0.0; 
  detMinMax();
  defNumSteps();
  mu = 0;
  for (int i=0; i<pointsGrid[0].size(); i++)
    for (int j=0; j<pointsGrid[1].size(); j++)
      for (int k=0; k<pointsGrid[2].size(); k++) {
        x = pointsGrid[0][i];
        y = pointsGrid[1][j];
        z = pointsGrid[2][k]; 
        den[mu] = readData.generateMoiSpherical(x, y, z);
        totalDen += den[mu];
        mu++;
      }
  cout << factorTot*totalDen << endl;
}






