#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include "WeightFunction.h"
using namespace std;

double Weights::normalize(double exponent) {
	double factor;
	factor = pow(2.0*exponent/M_PI, 0.75);
	return factor;
}

void Weights::Data(fstream& xyzFile) {

//  molec.readCoordinates(xyzFile);
  molec.getNumberOfDifAtoms();
  int j, atom;
  size_t l = molec.nNuc;
  size_t t = molec.numDifAtoms.size();
  wV.resize(t);
  wD.resize(l);

  for (int i=0; i<t; i++) {
    cout << "For the atom with atomic number: " << molec.numDifAtoms[i] << endl;
    wV[i].atomType = molec.numDifAtoms[i];
    cout << "Introduce the first coefficient of the atom: \n";
    cin >> wV[i].coeffWeight[0];
    cout << "Introduce the second coefficient of the atom: \n";
    cin >> wV[i].coeffWeight[1];
    cout << "Introduce the first exponent of the atom: \n";
    cin >> wV[i].expWeight[0];
    cout << "Introduce the second exponent of the atom: \n";
    cin >> wV[i].expWeight[1];
    cout << "\n\n";
		wV[i].normWeight[0] = normalize(wV[i].expWeight[0]);
		wV[i].normWeight[1] = normalize(wV[i].expWeight[1]);
  }

  for (int i=0; i<l; i++) {
    atom = molec.geom[i].nAtomic;
    j = 0;
    while (j<molec.numDifAtoms.size()) {
      if (atom == wV[j].atomType) {
        wD[i].atomType = atom;
        copy(begin(molec.geom[i].nCoor), end(molec.geom[i].nCoor), begin(wD[i].nucleusCoord));
        copy(begin(wV[j].coeffWeight), end(wV[j].coeffWeight), begin(wD[i].coeffWeight));
        copy(begin(wV[j].expWeight), end(wV[j].expWeight), begin(wD[i].expWeight));
        copy(begin(wV[j].normWeight), end(wV[j].normWeight), begin(wD[i].normWeight));
        j = t;
      }
      else
        j++;
    }
  }
/*
  for (int i=0; i<wD.size(); i++) {
    cout << wD[i].atomType << endl;
    cout << wD[i].nucleusCoord[0] << " " << wD[i].nucleusCoord[1] << " " << wD[i].nucleusCoord[2] << endl;
    for (int j=0; j<wD[i].coeffWeight.size(); j++)
      cout << " coeff: " << j << " " << wD[i].coeffWeight[j] << " " << " expo: " << wD[i].expWeight[j] << endl;
  }
*/

}

double Weights::constant () {
  double result;
  double result1;
  double Ax, Bx, Ay, By, Az, Bz;
  double Rab;
  size_t t = wD.size();

  result = 0.0;
  for (int i=0; i<t; i++) 
    for (int j=0; j<t; j++) {
      Ax = wD[i].nucleusCoord[0];
      Ay = wD[i].nucleusCoord[1];
      Az = wD[i].nucleusCoord[2];

      Bx = wD[j].nucleusCoord[0];
      By = wD[j].nucleusCoord[1];
      Bz = wD[j].nucleusCoord[2];

      Rab = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);
      Rab = sqrt(Rab);
			
			result1 = 0.0;
      for (int u=0; u<wD[i].coeffWeight.size(); u++)
      	for (int v=0; v<wD[j].coeffWeight.size(); v++) 
        	result1 += wD[i].normWeight[u]*wD[j].normWeight[v]*wD[i].coeffWeight[u]*wD[j].coeffWeight[v]*weightIntegral(wD[i].expWeight[u], wD[j].expWeight[v], Rab);

    result += result1;
  }

  return result;
}

double Weights::weightIntegral (double alpha, double beta, double Rab) {

  double eta, result, factor;
  double expos;
  double preExp;
  double ERF;
  
  expos = alpha*beta;
	factor = 2.0*pi25/(expos*sqrt(alpha+beta)); 
	if (Rab == 0)
		result = factor;
	else {
  	eta = expos/(alpha + beta);
  	eta = sqrt(eta);
  	ERF = erf(eta * Rab);
  	result = 0.5*sqrt(M_PI)*factor * ERF / (eta * Rab);
  }

  return result;
}

double Weights::weightFunction (double x1, double y1, double z1, double x2, double y2, double z2, double r12) {
  
  double Ax, Ay, Az, Bx, By, Bz;
  double result, coef1, expo1, coef2, expo2; 
  double X1a, Y1a, Z1a, X2b, Y2b, Z2b, R1a, R2b, Rab;
	double gf1, gf2;
  size_t t = wD.size();
  
	gf1 = 0.0;
	gf2 = 0.0;

//  for (int i=0; i<t; i++) 
  	for (int j=0; j<t; j++) {
  		Ax = wD[j].nucleusCoord[0];
   		Ay = wD[j].nucleusCoord[1];
   		Az = wD[j].nucleusCoord[2];

   		Bx = wD[j].nucleusCoord[0];
   		By = wD[j].nucleusCoord[1];
   		Bz = wD[j].nucleusCoord[2];

   		X1a = x1 - Ax;
   		Y1a = y1 - Ay;
   		Z1a = z1 - Az;
   		R1a = X1a*X1a + Y1a*Y1a + Z1a*Z1a;
     
   		X2b = x2 - Bx;
   		Y2b = y2 - By;
   		Z2b = z2 - Bz;
   		R2b = X2b*X2b + Y2b*Y2b + Z2b*Z2b;

 //   	for (int u=0; u<wD[i].coeffWeight.size(); u++) 
    		for (int v=0; v<wD[j].coeffWeight.size(); v++) {
					coef2 = wD[j].coeffWeight[v] * wD[j].normWeight[v];
					expo2 = wD[j].expWeight[v];
					coef1 = wD[j].coeffWeight[v] * wD[j].normWeight[v];
					expo1 = wD[j].expWeight[v];
    			gf1 += coef1 * exp(-expo1*R1a);
 					gf2 += coef2 * exp(-expo2*R2b);
 			}
	}
 	
	result = gf1 * gf2 / r12; 

//	cout << "voy: " << r12 << "   " << result << endl;

  return result;
}
