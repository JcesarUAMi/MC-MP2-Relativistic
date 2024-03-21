#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <random>
#include <cstdlib>
#include "RelaMetropolis.h"
using namespace std;

void RelaMetropolis::Data(fstream& xyzFile, fstream& MCFile)  {

  rB.rG.geomLecture(xyzFile);
  rB.weightFunctionData(MCFile);
  Eg = rB.constant();
//  cout << "CONSTANT: " << setprecision(6) << Eg << endl;
  numberWalkers = rB.nWalkers;
  nSteps = rB.nSteps;
  Nbsize = rB.Nbsize;
  wp.resize(numberWalkers);
  for (int i=0; i<numberWalkers; i++)
    wp[i].rSteps.resize(nSteps);

}


void RelaMetropolis::startWalking(double &x01, double &y01, double &z01, double &x02, double &y02, double &z02, double &wgt, double &r12) {
  
  int atom1, atom2;
  double difx, dify, difz, other;

  walker.random();  
	other = rB.wg.wD.size();
  atom1 = walker.randomNumber(other);
  atom2 = walker.randomNumber(other);
 
  x01 = rB.wg.wD[atom1].nucleusCoord[0];
  y01 = rB.wg.wD[atom1].nucleusCoord[1];
  z01 = rB.wg.wD[atom1].nucleusCoord[2];
  
  x02 = rB.wg.wD[atom2].nucleusCoord[0];
  y02 = rB.wg.wD[atom2].nucleusCoord[1];
  z02 = rB.wg.wD[atom2].nucleusCoord[2];

	walker.initialStep(x01, y01, z01, x02, y02, z02);
	
	difx = x01 - x02;
  dify = y01 - y02;
  difz = z01 - z02;
  
  r12 = sqrt(difx*difx + dify*dify + difz*difz);
 
  wgt = rB.weightFunction(x01, y01, z01, x02, y02, z02, r12);

//  cout << "que pues:   " << wgt << endl;

	wgt /= Eg;


}

void RelaMetropolis::motion(int mynode) {
	
	walker.random();

	for (int i=0; i<numberWalkers; i++) {
  	startWalking(wp[i].rSteps[0].e1[0], wp[i].rSteps[0].e1[1], wp[i].rSteps[0].e1[2], wp[i].rSteps[0].e2[0], wp[i].rSteps[0].e2[1], wp[i].rSteps[0].e2[2], wp[i].rSteps[0].w12, wp[i].rSteps[0].r12);  
//		cout << i << "   " << wp[i].rSteps[0].e1[0] << "  " << wp[i].rSteps[0].e1[1] << "   " << wp[i].rSteps[0].e1[2] << "   " << wp[i].rSteps[0].e2[0] << "  " << wp[i].rSteps[0].e2[1] << "   " << wp[i].rSteps[0].e2[2] << wp[i].rSteps[0].w12 << "   " << wp[i].rSteps[0].r12 << endl; 
	}

  burn(mynode); 

  moveScheme();

}

void RelaMetropolis::moveScheme() {
  
  double s1[3], s2[3];
  double r12, wgt;
  double ratio;
  double rval;
  int j;

  for (int n=0; n<numberWalkers; n++) {
  	reSteps = 0;
  	sSteps = 0;
  	fSteps = 0;
		j = 0;
		do { 
			r12 = 0;
      for (int i=0; i<3; i++) {
 //       s1[i] = wp[n].rSteps[j].e1[i] + walker.uniform(-moveLength, moveLength);
  //      s2[i] = wp[n].rSteps[j].e2[i] + walker.uniform(-moveLength, moveLength);
        s1[i] = wp[n].rSteps[j].e1[i] + (walker.uniform() - 0.5) * moveLength;
        s2[i] = wp[n].rSteps[j].e2[i] + (walker.uniform() - 0.5) * moveLength;
        r12 += (s1[i]-s2[i])*(s1[i]-s2[i]);
      }
      r12 = sqrt(r12);
			
      wgt = rB.weightFunction(s1[0], s1[1], s1[2], s2[0], s2[1], s2[2], r12) / Eg;

      ratio = wgt/wp[n].rSteps[j].w12;
      rval = walker.uniform();
			if (rval < 1E-3)
				rval = 1E-3;
      if (ratio > rval) {
        for (int u=0; u<3; u++) {
          wp[n].rSteps[j+1].e1[u] = s1[u];
          wp[n].rSteps[j+1].e2[u] = s2[u];
        }
      	wp[n].rSteps[j+1].w12 = wgt;
      	wp[n].rSteps[j+1].r12 = r12;
      	sSteps++;
				j++;
      } else
        fSteps++;

    if (reSteps == 1000) {
      moveLength = rescaleMoveLengthStep(fSteps, moveLength);
      sSteps = 0;
      fSteps = 0;
      reSteps = 0;
    }

    reSteps++;
		} while (j < nSteps);

	}
 
}


void RelaMetropolis::burn(int mynode) {
   
  double s1[3], s2[3];
  double r12, wgt;
  double ratio;
  double rval;
  int j;

  moveLength = 0.1;
  
  reSteps = 0;
  sSteps = 0;
  fSteps = 0;

  j = 0;
  do {
		for (int n=0; n<numberWalkers; n++) {
    	r12 = 0;
    	for (int i=0; i<3; i++) {
   //   	s1[i] = wp[n].rSteps[0].e1[i] + walker.uniform(-moveLength, moveLength);
   //   	s2[i] = wp[n].rSteps[0].e2[i] + walker.uniform(-moveLength, moveLength);
      	s1[i] = wp[n].rSteps[0].e1[i] + (walker.uniform() - 0.5) * moveLength;
      	s2[i] = wp[n].rSteps[0].e2[i] + (walker.uniform() - 0.5) * moveLength;
      	r12 += (s1[i]-s2[i])*(s1[i]-s2[i]);
  		}
    		r12 = sqrt(r12);
			
    		wgt = rB.weightFunction(s1[0], s1[1], s1[2], s2[0], s2[1], s2[2], r12) / Eg;
  
	  		ratio = wgt/wp[n].rSteps[0].w12;

				rval = walker.uniform();

				if (rval < 1E-3) 
					rval = 1E-3;
    		if (ratio > rval) {
      		for (int u=0; u<3; u++) {
        		wp[n].rSteps[0].e1[u] = s1[u];
        		wp[n].rSteps[0].e2[u] = s2[u];
      		}
      	wp[n].rSteps[0].w12 = wgt;
      	wp[n].rSteps[0].r12 = r12;
      	sSteps++; 
    		} else 
      		fSteps++;
			}

    if (reSteps == 1000) {
      moveLength = rescaleMoveLengthBurn(fSteps, moveLength);
      sSteps = 0;
      fSteps = 0;
      reSteps = 0;
    }
    
    reSteps++;
    j++;
  } while (j<100000); 
//  	cout << "is burning in the node  " << mynode << "   "   << j << " " << moveLength << endl;

}

double RelaMetropolis::rescaleMoveLengthBurn(int failedMoves, double moveLength) {

  double ratio = double((failedMoves*1.0)/(numberWalkers*1000));

//	cout << setprecision(16) << ratio << "  " << moveLength << endl;
  if (ratio < 0.5) 
    ratio = max(1.0/(2.0*ratio), 1.1);
  else
    ratio = max(0.9, 1.0/(2.0*ratio));

  return moveLength*ratio;

}

double RelaMetropolis::rescaleMoveLengthStep(int failedMoves, double moveLength) {

  double ratio = double((failedMoves*1.0)/1000);

//	cout << setprecision(6) << ratio << "  " << moveLength << endl;
  if (ratio < 0.5)
    ratio = max(1.0/(2.0*ratio), 1.1);
  else
    ratio = max(0.9, 1.0/(2.0*ratio));

  return moveLength*ratio;

}


